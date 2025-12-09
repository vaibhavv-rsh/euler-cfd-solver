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
    MUSCL,
    WENO
};

enum LIMITER
{
    NONE = 0,
    MINMOD,
    VANLEER
};

enum FACE
{
    WEST = 0,  // i-1/2
    EAST = 1,  // i+1/2
    SOUTH = 2, // j-1/2
    NORTH = 3  // j+1/2
};

struct CaseParameters
{
    int number_of_points_x;
    int number_of_points_y;
    double gamma;
    double domain_length_x;
    double domain_length_y;
    double end_time;
    double cfl;
    double dx;
    double dy;
    double time;
    int time_step;
};

class Euler
{
public:
    Euler(const std::string &params_file_name);
    ~Euler() = default;

    void solve();
    void print_solution();

private:
    void set_parameters(const std::string &params_file_name);
    void check_solution(std::string message = "");
    void init_solution();
    void update_solution();
    void piecewise_constant_reconstruction();
    void MUSCL_scheme();
    std::pair<std::array<double, 4>, std::array<double, 4>> compute_flux(int i, int j, FACE face);
    void Rusanov_Riemann_solver(std::array<double, 4> flux_left, std::array<double, 4> flux_right, int i, int j, FACE face);
    void integrate_solution_in_time();
    // void solve();

    SCHEME string_to_scheme_enum(const std::string &scheme)
    {
        if (scheme == "CONSTANT")
            return SCHEME::CONSTANT;
        if (scheme == "MUSCL")
            return SCHEME::MUSCL;
        throw std::invalid_argument("Unknown scheme: " + scheme);
    }
    LIMITER string_to_limiter_enum(const std::string &limiter)
    {
        if (limiter == "NONE")
            return LIMITER::NONE;
        if (limiter == "MINMOD")
            return LIMITER::MINMOD;
        if (limiter == "VANLEER")
            return LIMITER::VANLEER;
        throw std::invalid_argument("Unknown scheme: " + limiter);
    }

private:
    CaseParameters case_parameters;
    SCHEME numerical_scheme;
    LIMITER limiter;

    std::vector<double> x, y;                                               // centroids of finite volume cells
    std::vector<std::vector<std::array<double, 4>>> U, U_old;               // conserved variable vector U=[rho,rho*u,rho*v,E]
    std::vector<std::vector<std::array<std::array<double, 4>, 4>>> U_faces; // U at faces of cell (4 faces: west, east, south, north)
    std::vector<std::vector<std::array<std::array<double, 4>, 4>>> F_faces; // flux F at faces of cell (4 faces: west, east, south, north)
    double wave_speed_max, dt;
    std::string param_file;
};