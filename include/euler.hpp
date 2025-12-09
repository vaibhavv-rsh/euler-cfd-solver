#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

enum SCHEME
{
    CONSTANT = 0,
    MUSCL
};

enum LIMITER
{
    NONE = 0,
    MINMOD,
    VANLEER
};

enum FACE
{
    WEST = 0,//i-1/2
    EAST//i+1/2
};

struct CaseParameters
{
    int number_of_points;
    double gamma;
    double domain_length;
    double end_time;
    double cfl;
    double dx;
    double time;
    int time_step;
};

class Euler
{
public:
Euler(const std::string& params_file_name);
~Euler()=default;
void print_solution();

private:
void init_solution();
void update_solution();
void piecewise_constant_reconstruction();
void MUSCL_scheme();
std::pair<std::array<double, 3>, std::array<double, 3>> compute_flux(int x_index,int face);
void Rusanov_Riemann_solver(std::array<double,3>flux_left,std::array<double,3>flux_right,int x_index,int face);
void integrate_solution_in_time();
void solve();

SCHEME string_to_scheme_enum(const std::string& scheme){
    if(scheme=="CONSTANT") return SCHEME::CONSTANT;
    if(scheme=="MUSCL") return SCHEME::MUSCL;
    throw std::invalid_argument("Unknown scheme: " + scheme);
}
LIMITER string_to_limiter_enum(const std::string& limiter){
    if(limiter=="NONE") return LIMITER::NONE;
    if(limiter=="MINMOD") return LIMITER::MINMOD;
    if(limiter=="VANLEER") return LIMITER::VANLEER;
    throw std::invalid_argument("Unknown scheme: " + limiter);
}

private:
CaseParameters case_parameters;
SCHEME numerical_scheme;
LIMITER limiter;

std::vector <double> x;//centroids of finite volume cells
std::vector <std::array<double, 3>> U,U_old;//conserved variable vector U=[rho,rho*vel,E]
std::vector <std::array<std::array<double, 3>, 2>>U_faces;//U at faces of cell
std::vector <std::array<std::array<double, 3>, 2>> F_faces;//flux F at faces of cell
double wave_speed_max,dt;

};