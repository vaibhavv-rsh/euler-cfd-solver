#include "euler_2d.hpp"
#include "param_loader.hpp"
#include "utils.hpp"

#include <omp.h>

Euler::Euler(const std::string &params_file_name)
{
    this->set_parameters(params_file_name);
    this->init_solution();
}
void Euler::set_parameters(const std::string &params_file_name)
{
    param_loader::ParamLoader &param_loader = param_loader::ParamLoader::getInstance(params_file_name);
    this->param_file = params_file_name;
    this->case_parameters.number_of_points_x = param_loader.get_params<int>("case_parameters", "number_of_points_x");
    this->case_parameters.number_of_points_y = param_loader.get_params<int>("case_parameters", "number_of_points_y");
    this->case_parameters.gamma = param_loader.get_params<double>("case_parameters", "gamma");
    this->case_parameters.domain_length_x = param_loader.get_params<double>("case_parameters", "domain_length_x");
    this->case_parameters.domain_length_y = param_loader.get_params<double>("case_parameters", "domain_length_y");
    this->case_parameters.end_time = param_loader.get_params<double>("case_parameters", "end_time");
    this->case_parameters.cfl = param_loader.get_params<double>("case_parameters", "cfl");
    this->case_parameters.dx = this->case_parameters.domain_length_x / (case_parameters.number_of_points_x - 1);
    this->case_parameters.dy = this->case_parameters.domain_length_y / (case_parameters.number_of_points_y - 1);
    this->case_parameters.time = param_loader.get_params<double>("case_parameters", "time");
    this->case_parameters.time_step = param_loader.get_params<int>("case_parameters", "time_step");

    this->numerical_scheme = this->string_to_scheme_enum(param_loader.get_params<std::string>("case_parameters", "numerical_scheme"));
    this->limiter = this->string_to_limiter_enum(param_loader.get_params<std::string>("case_parameters", "limiter"));
}
void Euler::init_solution()
{
    // Create mesh
    this->x = std::vector<double>(case_parameters.number_of_points_x);
    this->y = std::vector<double>(case_parameters.number_of_points_y);
    this->U = std::vector<std::vector<std::array<double, 4>>>(case_parameters.number_of_points_x,
                                                              std::vector<std::array<double, 4>>(case_parameters.number_of_points_y));
    this->U_old = std::vector<std::vector<std::array<double, 4>>>(case_parameters.number_of_points_x,
                                                                  std::vector<std::array<double, 4>>(case_parameters.number_of_points_y));
    this->U_faces = std::vector<std::vector<std::array<std::array<double, 4>, 4>>>(case_parameters.number_of_points_x,
                                                                                   std::vector<std::array<std::array<double, 4>, 4>>(case_parameters.number_of_points_y));
    this->F_faces = std::vector<std::vector<std::array<std::array<double, 4>, 4>>>(case_parameters.number_of_points_x,
                                                                                   std::vector<std::array<std::array<double, 4>, 4>>(case_parameters.number_of_points_y));

    for (int i = 0; i < this->case_parameters.number_of_points_x; i++)
    {
        this->x[i] = i * case_parameters.dx;
    }
    for (int j = 0; j < this->case_parameters.number_of_points_y; j++)
    {
        this->y[j] = j * case_parameters.dy;
    }

    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double p = 0.0;

    for (int i = 0; i < this->case_parameters.number_of_points_x; i++)
    {
        for (int j = 0; j < this->case_parameters.number_of_points_y; j++)
        {
            if (this->x[i] < 0.5)
            {
                rho = 1.0;
                u = 0.0;
                v = 0.0;
                p = 1.0;
            }
            else
            {
                rho = 0.125;
                u = 0.0;
                v = 0.0;
                p = 0.1;
            }
            if (i == 0)
                u = 10.0;
            double E = p / (case_parameters.gamma - 1) + 0.5 * rho * (std::pow(u, 2) + std::pow(v, 2));

            this->U[i][j][0] = rho;
            this->U[i][j][1] = rho * u;
            this->U[i][j][2] = rho * v;
            this->U[i][j][3] = E;

            this->U_old[i][j][0] = rho;
            this->U_old[i][j][1] = rho * u;
            this->U_old[i][j][2] = rho * v;
            this->U_old[i][j][3] = E;
        }
    }
}

void Euler::MUSCL_scheme()
{
    for (int variable = 0; variable < 4; variable++)
    {
#pragma omp parallel for
        for (int j = 0; j < case_parameters.number_of_points_y; j++)
        {
            this->U_faces[0][j][FACE::WEST][variable] = this->U[0][j][variable];
            this->U_faces[0][j][FACE::EAST][variable] = this->U[0][j][variable];
        }
#pragma omp parallel for
        for (int j = 0; j < case_parameters.number_of_points_y; j++)
        {
            this->U_faces[case_parameters.number_of_points_x - 1][j][FACE::WEST][variable] =
                U[case_parameters.number_of_points_x - 1][j][variable];
            this->U_faces[case_parameters.number_of_points_x - 1][j][FACE::EAST][variable] =
                U[case_parameters.number_of_points_x - 1][j][variable];
        }
#pragma omp parallel for
        for (int i = 0; i < case_parameters.number_of_points_x; i++)
        {
            this->U_faces[i][0][FACE::SOUTH][variable] = U[i][0][variable];
            this->U_faces[i][0][FACE::NORTH][variable] = U[i][0][variable];
        }
#pragma omp parallel for
        for (int i = 0; i < case_parameters.number_of_points_x; i++)
        {
            this->U_faces[i][case_parameters.number_of_points_y - 1][FACE::SOUTH][variable] =
                U[i][case_parameters.number_of_points_y - 1][variable];
            this->U_faces[i][case_parameters.number_of_points_y - 1][FACE::NORTH][variable] =
                U[i][case_parameters.number_of_points_y - 1][variable];
        }
    }

#pragma omp parallel for collapse(3)
    for (int i = 1; i < case_parameters.number_of_points_x - 1; i++)
    {
        for (int j = 0; j < case_parameters.number_of_points_y; j++)
        {
            for (int variable = 0; variable < 4; variable++)
            {
                double du_i_plus_half = this->U[i + 1][j][variable] - this->U[i][j][variable];
                double du_i_minus_half = this->U[i][j][variable] - this->U[i - 1][j][variable];

                double r_left = du_i_minus_half / (du_i_plus_half + 1e-8);
                double r_right = du_i_plus_half / (du_i_minus_half + 1e-8);

                double psi_left = 1.0;
                double psi_right = 1.0;

                if (this->limiter == LIMITER::MINMOD)
                {
                    psi_left = std::max(0.0, std::min(1.0, r_left));
                    psi_right = std::max(0.0, std::min(1.0, r_right));
                }
                else if (this->limiter == LIMITER::VANLEER)
                {
                    psi_left = (r_left + std::fabs(r_left)) / (1.0 + std::fabs(r_left));
                    psi_right = (r_right + std::fabs(r_right)) / (1.0 + std::fabs(r_right));
                }

                this->U_faces[i][j][FACE::WEST][variable] = this->U[i][j][variable] - 0.5 * psi_left * du_i_plus_half;
                this->U_faces[i][j][FACE::EAST][variable] = this->U[i][j][variable] + 0.5 * psi_right * du_i_minus_half;
            }
        }
    }

    for (int i = 0; i < case_parameters.number_of_points_x; i++)
    {
        for (int j = 1; j < case_parameters.number_of_points_y - 1; j++)
        {
            for (int variable = 0; variable < 4; variable++)
            {
                double dv_j_plus_half = this->U[i][j + 1][variable] - this->U[i][j][variable];
                double dv_j_minus_half = this->U[i][j][variable] - this->U[i][j - 1][variable];

                double r_left = dv_j_minus_half / (dv_j_plus_half + 1e-8);
                double r_right = dv_j_plus_half / (dv_j_minus_half + 1e-8);

                double psi_left = 1.0;
                double psi_right = 1.0;

                if (this->limiter == LIMITER::MINMOD)
                {
                    psi_left = std::max(0.0, std::min(1.0, r_left));
                    psi_right = std::max(0.0, std::min(1.0, r_right));
                }
                else if (this->limiter == LIMITER::VANLEER)
                {
                    psi_left = (r_left + std::fabs(r_left)) / (1.0 + std::fabs(r_left));
                    psi_right = (r_right + std::fabs(r_right)) / (1.0 + std::fabs(r_right));
                }
                this->U_faces[i][j][FACE::SOUTH][variable] = this->U[i][j][variable] - 0.5 * psi_left * dv_j_plus_half;
                this->U_faces[i][j][FACE::NORTH][variable] = this->U[i][j][variable] + 0.5 * psi_right * dv_j_minus_half;
            }
        }
    }
}

void Euler::update_solution()
{
    this->U_old = this->U;

    double speed_max_local = 0.0;
#pragma omp parallel for collapse(2) reduction(max : speed_max_local)
    for (int i = 0; i < this->case_parameters.number_of_points_x; i++)
    {
        for (int j = 0; j < this->case_parameters.number_of_points_y; j++)
        {
            double rho = this->U[i][j][0];
            if (rho < 1e-6)
            {
                std::cout << "rho is zero at cell (" << i << ", " << j << ")" << std::endl;
                exit(1);
            }

            double u = this->U[i][j][1] / rho;
            double v = this->U[i][j][2] / rho;
            double E = this->U[i][j][3];
            double p = (case_parameters.gamma - 1.0) * (E - 0.5 * rho * (std::pow(u, 2) + std::pow(v, 2)));

            if (p < 0)
            {
                std::cout << "Negative pressure at cell (" << i << ", " << j << ")" << std::endl;
                exit(1);
            }

            double speed_of_sound = std::sqrt(case_parameters.gamma * p / rho);
            double speed = std::max(std::abs(u), std::abs(v)) + speed_of_sound;
            speed_max_local = std::max(speed_max_local, speed);
        }
    }

    this->dt = (case_parameters.cfl * case_parameters.dx) / speed_max_local;
}

void Euler::Rusanov_Riemann_solver(std::array<double, 4> flux_left, std::array<double, 4> flux_right, int i, int j, FACE face)
{
    int index_offset = face;
    for (int variable = 0; variable < 4; variable++)
    {
        if (face == FACE::WEST || face == FACE::EAST)
        {
            const double &q_left = this->U_faces[i - 1 + index_offset][j][face][variable];
            const double &q_right = U_faces[i + index_offset][j][face][variable];
            const double &f_left = flux_left[variable];
            const double &f_right = flux_right[variable];
            if (std::isnan(q_left) || std::isnan(q_right) || std::isnan(f_left) || std::isnan(f_right))
            {
                check_solution("q_left,q_right,f_left,f_right");
                exit(1);
            }
            this->F_faces[i][j][face][variable] = 0.5 * (f_left + f_right) - this->wave_speed_max * (q_right - q_left);
            if (std::isnan(this->F_faces[i][j][face][variable]))
            {
                std::cout << "q_left: " << q_left << " q_right: " << q_right << " f_left: " << f_left << " f_right: " << f_right << std::endl;
                std::cout << "wave_speed_max: " << this->wave_speed_max << std::endl;
                exit(1);
            }
        }
        else
        {
            int index_offset = face - 2;
            const double &q_left = this->U_faces[i][j - 1 + index_offset][face][variable];
            const double &q_right = U_faces[i][j + index_offset][face][variable];
            const double &f_left = flux_left[variable];
            const double &f_right = flux_right[variable];
            this->F_faces[i][j][face][variable] = 0.5 * (f_left + f_right) - this->wave_speed_max * (q_right - q_left);
        }
    }
}

std::pair<std::array<double, 4>, std::array<double, 4>> Euler::compute_flux(int i, int j, FACE face)
{
    std::array<double, 4> flux_left = {0.0, 0.0, 0.0, 0.0};
    std::array<double, 4> flux_right = {0.0, 0.0, 0.0, 0.0};

    double rho_left = 0.0;
    double u_left = 0.0;
    double v_left = 0.0;
    double E_left = 0.0;
    double p_left = 0.0;
    double rho_right = 0.0;
    double u_right = 0.0;
    double v_right = 0.0;
    double E_right = 0.0;
    double p_right = 0.0;
    if (face == FACE::WEST || face == FACE::EAST)
    {
        rho_left = this->U_faces[i - 1 + face][j][FACE::EAST][0];
        u_left = this->U_faces[i - 1 + face][j][FACE::EAST][1] / rho_left;
        v_left = this->U_faces[i - 1 + face][j][FACE::EAST][2] / rho_left;
        E_left = this->U_faces[i - 1 + face][j][FACE::EAST][3];
        p_left = (case_parameters.gamma - 1.0) * (E_left - 0.5 * rho_left * (std::pow(u_left, 2) + std::pow(v_left, 2)));
        rho_right = this->U_faces[i + face][j][FACE::WEST][0];
        u_right = this->U_faces[i + face][j][FACE::WEST][1] / rho_right;
        v_right = this->U_faces[i + face][j][FACE::WEST][2] / rho_right;
        E_right = this->U_faces[i + face][j][FACE::WEST][3];
        p_right = (case_parameters.gamma - 1.0) * (E_right - 0.5 * rho_right * (std::pow(u_right, 2) + std::pow(v_right, 2)));
    }
    if (face == FACE::SOUTH || face == FACE::NORTH)
    {
        int index_offset = face - 2;
        rho_left = this->U_faces[i][j - 1 + index_offset][FACE::NORTH][0];
        u_left = this->U_faces[i][j - 1 + index_offset][FACE::NORTH][1] / rho_left;
        v_left = this->U_faces[i][j - 1 + index_offset][FACE::NORTH][2] / rho_left;
        E_left = this->U_faces[i][j - 1 + index_offset][FACE::NORTH][3];
        p_left = (case_parameters.gamma - 1.0) * (E_left - 0.5 * rho_left * (std::pow(u_left, 2) + std::pow(v_left, 2)));
        rho_right = this->U_faces[i][j + index_offset][FACE::SOUTH][0];
        u_right = this->U_faces[i][j + index_offset][FACE::SOUTH][1] / rho_right;
        v_right = this->U_faces[i][j + index_offset][FACE::SOUTH][2] / rho_right;
        E_right = this->U_faces[i][j + index_offset][FACE::SOUTH][3];
        p_right = (case_parameters.gamma - 1.0) * (E_right - 0.5 * rho_right * (std::pow(u_right, 2) + std::pow(v_right, 2)));
    }
    double a_left = std::sqrt(case_parameters.gamma * p_left / rho_left);
    double a_right = std::sqrt(case_parameters.gamma * p_right / rho_right);
    if (face == FACE::WEST || face == FACE::EAST)
    {
        flux_left[0] = rho_left * u_left;
        flux_left[1] = rho_left * std::pow(u_left, 2) + p_left;
        flux_left[2] = rho_left * u_left * v_left;
        flux_left[3] = (E_left + p_left) * u_left;

        flux_right[0] = rho_right * u_right;
        flux_right[1] = rho_right * std::pow(u_right, 2) + p_right;
        flux_right[2] = rho_right * u_right * v_right;
        flux_right[3] = (E_right + p_right) * u_right;
        this->wave_speed_max = std::max(std::fabs(u_left) + a_left, std::fabs(u_right) + a_right);
    }
    else
    {
        flux_left[0] = rho_left * v_left;
        flux_left[1] = rho_left * u_left * v_left;
        flux_left[2] = rho_left * std::pow(v_left, 2) + p_left;
        flux_left[3] = (E_left + p_left) * v_left;

        flux_right[0] = rho_right * v_right;
        flux_right[1] = rho_right * u_right * v_right;
        flux_right[2] = rho_right * std::pow(v_right, 2) + p_right;
        flux_right[3] = (E_right + p_right) * v_right;
        this->wave_speed_max = std::max(std::fabs(v_left) + a_left, std::fabs(v_right) + a_right);
    }
    return std::make_pair(flux_left, flux_right);
}

void Euler::integrate_solution_in_time()
{
#pragma omp parallel for collapse(2)
    for (int i = 0; i < this->case_parameters.number_of_points_x; i++)
    {
        for (int j = 0; j < this->case_parameters.number_of_points_y; j++)
        {

            double dFx[4] = {
                this->F_faces[i][j][FACE::EAST][0] - this->F_faces[i][j][FACE::WEST][0],
                this->F_faces[i][j][FACE::EAST][1] - this->F_faces[i][j][FACE::WEST][1],
                this->F_faces[i][j][FACE::EAST][2] - this->F_faces[i][j][FACE::WEST][2],
                this->F_faces[i][j][FACE::EAST][3] - this->F_faces[i][j][FACE::WEST][3]};

            double dFy[4] = {
                this->F_faces[i][j][FACE::NORTH][0] - this->F_faces[i][j][FACE::SOUTH][0],
                this->F_faces[i][j][FACE::NORTH][1] - this->F_faces[i][j][FACE::SOUTH][1],
                this->F_faces[i][j][FACE::NORTH][2] - this->F_faces[i][j][FACE::SOUTH][2],
                this->F_faces[i][j][FACE::NORTH][3] - this->F_faces[i][j][FACE::SOUTH][3]};

            for (int variable = 0; variable < 4; variable++)
            {
                this->U[i][j][variable] = this->U_old[i][j][variable] -
                                          (this->dt / case_parameters.dx) * dFx[variable] -
                                          (this->dt / case_parameters.dy) * dFy[variable];
            }
        }
    }
}
void Euler::check_solution(std::string message)
{
    for (int i = 0; i < case_parameters.number_of_points_x; i++)
    {
        for (int j = 0; j < case_parameters.number_of_points_y; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                if (std::isnan(this->U[i][j][k]))
                {
                    std::cout << message << "\n";
                    std::cout << "Warning: NaN detected at cell (" << i << ", " << j << ")" << std::endl;
                    std::cout << "Time: " << this->case_parameters.time << ", dt: " << this->dt << std::endl;
                    std::cout << "Solution values: " << std::endl;
                    std::cout << "rho: " << this->U[i][j][0] << ", rho*u: " << this->U[i][j][1]
                              << ", rho*v: " << this->U[i][j][2] << ", E: " << this->U[i][j][3] << std::endl;
                    std::cout << "Exiting simulation..." << std::endl;
                    exit(1);
                }
            }
        }
    }
}
void Euler::solve()
{

    // Save initial solution
    std::stringstream ss;
    ss << "solution_"
       << std::setw(6) << std::setfill('0') << this->case_parameters.number_of_points_x << "x"
       << std::setw(6) << std::setfill('0') << this->case_parameters.number_of_points_y << "_"
       << std::setw(6) << std::setfill('0') << static_cast<int>(this->case_parameters.time) << ".csv";

    std::string filename = "/home/vaibhavvikas/euler_cfd/continuous_output/" + ss.str();
    while (this->case_parameters.time < case_parameters.end_time)
    {
        if (interrupt_requested)
        {
            std::cout << "Exiting early due to interrupt. Saving partial results...\n";
            this->print_solution();
            break;
        }
        this->update_solution();
        if (this->numerical_scheme == SCHEME::CONSTANT)
        {
            this->U_faces[0][0][FACE::SOUTH] = this->U[this->case_parameters.number_of_points_x - 1][this->case_parameters.number_of_points_y - 2];
            this->U_faces[this->case_parameters.number_of_points_x - 1][this->case_parameters.number_of_points_y - 1][FACE::NORTH] = this->U[this->case_parameters.number_of_points_x - 1][1];

            for (int i = 1; i < this->case_parameters.number_of_points_x - 1; i++)
            {
                for (int j = 1; j < this->case_parameters.number_of_points_y - 1; j++)
                {
                    for (int variable = 0; variable < 4; variable++)
                    {
                        this->U_faces[i][j][FACE::WEST][variable] = this->U[i - 1][j][variable];
                        this->U_faces[i][j][FACE::EAST][variable] = this->U[i][j][variable];
                        this->U_faces[i][j][FACE::SOUTH][variable] = this->U[i][j - 1][variable];
                        this->U_faces[i][j][FACE::NORTH][variable] = this->U[i][j][variable];
                    }
                }
            }
        }
        else if (this->numerical_scheme == SCHEME::MUSCL)
        {
            this->MUSCL_scheme();
        }
#pragma omp parallel for collapse(2)
        for (int i = 1; i < this->case_parameters.number_of_points_x - 1; i++)
        {
            for (int j = 1; j < this->case_parameters.number_of_points_y - 1; j++)
            {
                auto flux_pair = this->compute_flux(i, j, FACE::WEST);
                this->Rusanov_Riemann_solver(flux_pair.first, flux_pair.second, i, j, FACE::WEST);
                flux_pair = this->compute_flux(i, j, FACE::EAST);
                this->Rusanov_Riemann_solver(flux_pair.first, flux_pair.second, i, j, FACE::EAST);
                flux_pair = this->compute_flux(i, j, FACE::SOUTH);
                this->Rusanov_Riemann_solver(flux_pair.first, flux_pair.second, i, j, FACE::SOUTH);
                flux_pair = this->compute_flux(i, j, FACE::NORTH);
                this->Rusanov_Riemann_solver(flux_pair.first, flux_pair.second, i, j, FACE::NORTH);
            }
        }
        this->integrate_solution_in_time();
        this->case_parameters.time += this->dt;
        this->case_parameters.time_step++;

        std::cout << "Time: " << this->case_parameters.time << ", dt: " << this->dt << "\n";
        this->print_solution();
    }
}

void Euler::print_solution()
{
    std::string output_dir = "../solutions_2d/";
    std::string mkdir_command = "mkdir -p " + output_dir;
    std::system(mkdir_command.c_str());
    std::ostringstream time_step_temp, points_x_temp, points_y_temp;
    time_step_temp << std::setfill('0') << std::setw(6);
    time_step_temp << this->case_parameters.time_step;
    auto time_step = time_step_temp.str();

    points_x_temp << std::setfill('0') << std::setw(6);
    points_x_temp << case_parameters.number_of_points_x;
    auto points_x = points_x_temp.str();

    points_y_temp << std::setfill('0') << std::setw(6);
    points_y_temp << case_parameters.number_of_points_y;
    auto points_y = points_y_temp.str();

    std::string csv_file = output_dir + "solution_" + points_x + "x" + points_y + "_" + time_step + ".csv";
    std::ofstream csv_output;
    csv_output.open(csv_file);
    csv_output << "x,y,rho,u,v,p" << std::endl;
    for (int i = 0; i < case_parameters.number_of_points_x; i++)
    {
        for (int j = 0; j < case_parameters.number_of_points_y; j++)
        {
            double rho = this->U[i][j][0];
            double u = this->U[i][j][1] / rho;
            double v = this->U[i][j][2] / rho;
            double p = (case_parameters.gamma - 1) * (this->U[i][j][3] - 0.5 * rho * (std::pow(u, 2) + std::pow(v, 2)));
            csv_output << this->x[i] << "," << this->y[j] << "," << rho << "," << u << "," << v << "," << p << std::endl;
        }
    }
    csv_output.close();

    std::string vtk_file = output_dir + "solution_" + points_x + "x" + points_y + "_" + time_step + ".vtk";
    std::ofstream vtk_output;
    vtk_output.open(vtk_file);
    vtk_output << "# vtk DataFile Version 2.0\n";
    vtk_output << "2D Euler Solution\n";
    vtk_output << "ASCII\n";
    vtk_output << "DATASET STRUCTURED_GRID\n";
    vtk_output << "DIMENSIONS " << case_parameters.number_of_points_x << " " << case_parameters.number_of_points_y << " 1\n";
    vtk_output << "POINTS " << case_parameters.number_of_points_x * case_parameters.number_of_points_y << " float\n";

    for (int i = 0; i < case_parameters.number_of_points_x; i++)
    {
        for (int j = 0; j < case_parameters.number_of_points_y; j++)
        {
            vtk_output << this->x[i] << " " << this->y[j] << " 0\n";
        }
    }

    vtk_output << "CELL_DATA " << (case_parameters.number_of_points_x - 1) * (case_parameters.number_of_points_y - 1) << "\n";
    vtk_output << "SCALARS rho float 1\nLOOKUP_TABLE default\n";
    for (int i = 0; i < case_parameters.number_of_points_x - 1; i++)
    {
        for (int j = 0; j < case_parameters.number_of_points_y - 1; j++)
        {
            vtk_output << this->U[i][j][0] << "\n";
        }
    }

    vtk_output << "VECTORS velocity float\n";
    for (int i = 0; i < case_parameters.number_of_points_x - 1; i++)
    {
        for (int j = 0; j < case_parameters.number_of_points_y - 1; j++)
        {
            double rho = this->U[i][j][0];
            double u = this->U[i][j][1] / rho;
            double v = this->U[i][j][2] / rho;
            vtk_output << u << " " << v << " 0\n";
        }
    }

    vtk_output << "SCALARS pressure float 1\nLOOKUP_TABLE default\n";
    for (int i = 0; i < case_parameters.number_of_points_x - 1; i++)
    {
        for (int j = 0; j < case_parameters.number_of_points_y - 1; j++)
        {
            double rho = this->U[i][j][0];
            double u = this->U[i][j][1] / rho;
            double v = this->U[i][j][2] / rho;
            double p = (case_parameters.gamma - 1) * (this->U[i][j][3] - 0.5 * rho * (std::pow(u, 2) + std::pow(v, 2)));
            vtk_output << p << "\n";
        }
    }

    vtk_output.close();

    std::cout << "Solution written to: " << csv_file << " and " << vtk_file << std::endl;
}
int main(int argc, char **argv)
{
    std::signal(SIGINT, signal_handler);

#pragma omp parallel
    {
#pragma omp master
        std::cout << "Running with " << omp_get_num_threads() << " threads.\n";
    }
    Euler euler("../config/params_2d.yaml");
    double t1 = omp_get_wtime();
    euler.solve();
    double t2 = omp_get_wtime();
    std::cout << "Simulation took: " << (t2 - t1) << " seconds" << std::endl;

    return 0;
}