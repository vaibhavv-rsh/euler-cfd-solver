#include "euler_1d.hpp"
#include "param_loader.hpp"

Euler::Euler(const std::string &params_file_name)
{
    // load parameters from params_file
    param_loader::ParamLoader &param_loader = param_loader::ParamLoader::getInstance(params_file_name);
    this->case_parameters.number_of_points = param_loader.get_params<int>("case_parameters", "number_of_points");
    this->case_parameters.gamma = param_loader.get_params<double>("case_parameters", "gamma");
    this->case_parameters.domain_length = param_loader.get_params<double>("case_parameters", "domain_length");
    this->case_parameters.end_time = param_loader.get_params<double>("case_parameters", "end_time");
    this->case_parameters.cfl = param_loader.get_params<double>("case_parameters", "cfl");
    this->case_parameters.dx = this->case_parameters.domain_length / (case_parameters.number_of_points - 1);
    this->case_parameters.time = param_loader.get_params<double>("case_parameters", "time");
    this->case_parameters.time_step = param_loader.get_params<int>("case_parameters", "time_step");

    this->numerical_scheme = this->string_to_scheme_enum(param_loader.get_params<std::string>("case_parameters", "numerical_scheme"));
    this->limiter = this->string_to_limiter_enum(param_loader.get_params<std::string>("case_parameters", "limiter"));

    this->x = std::vector<double>(case_parameters.number_of_points);
    this->U = std::vector<std::array<double, 3>>(case_parameters.number_of_points);
    this->U_faces = std::vector<std::array<std::array<double, 3>, 2>>(case_parameters.number_of_points);
    this->F_faces = std::vector<std::array<std::array<double, 3>, 2>>(case_parameters.number_of_points);

    // initiliaze_solution
    this->init_solution();
}
void Euler::init_solution()
{
    // Create/read mesh
    for (int i = 0; i < this->case_parameters.number_of_points; i++)
    {
        x[i] = i * case_parameters.dx;
    }

    double rho = 0.0;
    double u = 0.0;
    double p = 0.0;

    for (int i = 0; i < this->case_parameters.number_of_points; i++)
    {
        if (this->x[i] < 0.5)
        {
            rho = 1.0;
            u = 0.0;
            p = 1.0;
        }
        else
        {
            rho = 0.125;
            u = 0.0;
            p = 0.1;
        }
        U[i][0] = rho;
        U[i][1] = rho * u;
        U[i][2] = p / (case_parameters.gamma - 1) + 0.5 * rho * std::pow(u, 2);
    }
}

void Euler::update_solution()
{
    this->U_old = U;
    double speed_max = 0.0;
    for (int i = 0; i < this->case_parameters.number_of_points; i++)
    {
        double rho = this->U[i][0];
        if (rho < 1e-6)
        {
            std::cout << "rho is zero\n";
        }
        double u = this->U[i][1] / rho;
        double p = (case_parameters.gamma - 1.0) * (U[i][2] - 0.5 * rho * std::pow(u, 2));
        // if(std::isnan(p)){std::cout<<"p is nan\n";}
        double speed_of_sound = std::sqrt(case_parameters.gamma * p / rho);
        if (p < 0)
        {
            std::cout << "p is less than 0\n";
        }
        if (speed_of_sound + std::fabs(u) > speed_max)
            speed_max = speed_of_sound + std::fabs(u);
    }
    this->dt = (case_parameters.cfl * case_parameters.dx) / speed_max;
}

void Euler::MUSCL_scheme()
{
    // use lower‐order scheme near boundaries
    for (int variable = 0; variable < 3; variable++)
    {
        this->U_faces[0][FACE::WEST][variable] = U[0][variable];
        this->U_faces[0][FACE::EAST][variable] = U[0][variable];
        this->U_faces[case_parameters.number_of_points - 1][FACE::WEST][variable] = this->U[case_parameters.number_of_points - 1][variable];
        this->U_faces[case_parameters.number_of_points - 1][FACE::EAST][variable] = this->U[case_parameters.number_of_points - 1][variable];
    }
    for (int i = 1; i < this->case_parameters.number_of_points - 1; i++)
    {
        for (int variable = 0; variable < 3; variable++)
        {
            double du_i_plus_half = this->U[i + 1][variable] - U[i][variable];
            double du_i_minus_half = this->U[i][variable] - U[i - 1][variable];

            double r_left = du_i_minus_half / (du_i_plus_half + 1e-8);
            double r_right = du_i_plus_half / (du_i_minus_half + 1e-8);

            double psi_left = 1.0;
            double psi_right = 1.0;

            if (this->limiter == LIMITER::MINMOD)
            {
                psi_left = std::max(0.0, std::min(1.0, r_left));
                psi_right = std::max(0.0, std::min(1.0, r_right));
            }
            else if (limiter == LIMITER::VANLEER)
            {
                psi_left = (r_left + std::fabs(r_left)) / (1.0 + std::fabs(r_left));
                psi_right = (r_right + std::fabs(r_right)) / (1.0 + std::fabs(r_right));
            }
            this->U_faces[i][FACE::WEST][variable] = U[i][variable] - 0.5 * psi_left * du_i_plus_half;
            this->U_faces[i][FACE::EAST][variable] = U[i][variable] + 0.5 * psi_right * du_i_minus_half;
        }
    }
}
std::pair<std::array<double, 3>, std::array<double, 3>> Euler::compute_flux(int i, int face)
{
    int index_offset = face;
    std::array<double, 3> flux_left, flux_right;
    double rho_left = this->U_faces[i - 1 + index_offset][FACE::EAST][0];
    double u_left = this->U_faces[i - 1 + index_offset][FACE::EAST][1] / rho_left;
    double E_left = this->U_faces[i - 1 + index_offset][FACE::EAST][2];
    double p_left = (this->case_parameters.gamma - 1.0) * (E_left - 0.5 * rho_left * std::pow(u_left, 2));
    double a_left = std::sqrt(case_parameters.gamma * p_left / rho_left);

    double rho_right = this->U_faces[i + index_offset][FACE::WEST][0];
    double u_right = this->U_faces[i + index_offset][FACE::WEST][1] / rho_right;
    double E_right = this->U_faces[i + index_offset][FACE::WEST][2];
    double p_right = (this->case_parameters.gamma - 1.0) * (E_right - 0.5 * rho_right * std::pow(u_right, 2));
    double a_right = std::sqrt(case_parameters.gamma * p_right / rho_right);

    flux_left[0] = rho_left * u_left;
    flux_left[1] = p_left + rho_left * std::pow(u_left, 2);
    flux_left[2] = u_left * (E_left + p_left);

    flux_right[0] = rho_right * u_right;
    flux_right[1] = p_right + rho_right * std::pow(u_right, 2);
    flux_right[2] = u_right * (E_right + p_right);

    this->wave_speed_max = std::max(std::fabs(u_left) + a_left, std::fabs(u_right) + a_right);
    std::cout<<"flux_left: "<<flux_left[0]<<" "<<flux_left[1]<<" "<<flux_left[2]<<std::endl;
    std::cout<<"flux_right: "<<flux_right[0]<<" "<<flux_right[1]<<" "<<flux_right[2]<<std::endl;
    return {flux_left, flux_right};
}
void Euler::Rusanov_Riemann_solver(std::array<double, 3> flux_left, std::array<double, 3> flux_right, int i, int face)
{
    int index_offset = face;
    for (int variable = 0; variable < 3; variable++)
    {
        const double &q_left = this->U_faces[i - 1 + index_offset][FACE::EAST][variable];
        const double &q_right = U_faces[i + index_offset][FACE::WEST][variable];
        const double &f_left = flux_left[variable];
        const double &f_right = flux_right[variable];
        this->F_faces[i][face][variable] = 0.5 * (f_left + f_right) - this->wave_speed_max * (q_right - q_left);
    }
}
void Euler::integrate_solution_in_time()
{
    for (int i = 0; i < this->case_parameters.number_of_points; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            const double &dF = this->F_faces[i][FACE::EAST][j] - this->F_faces[i][FACE::WEST][j];
            this->U[i][j] = this->U_old[i][j] - (this->dt / case_parameters.dx) * dF;
        }
    }
}
void Euler::solve()
{
    std::cout << "numerical_scheme: " << this->numerical_scheme << "\n";
    std::cout << "numerical_limiter: " << this->limiter << "\n";
    while (this->case_parameters.time < case_parameters.end_time)
    {
        this->update_solution();
        if (this->numerical_scheme == SCHEME::CONSTANT)
        { // piecewise_constant reconstruction
            for (int i = 0; i < this->case_parameters.number_of_points; i++)
            {
                for (int variable = 0; variable < 3; variable++)
                {
                    this->U_faces[i][FACE::WEST][variable] = U[i][variable];
                    this->U_faces[i][FACE::EAST][variable] = U[i][variable];
                }
            }
        }
        else if (this->numerical_scheme == SCHEME::MUSCL)
        {
            // use high‐resolution MUSCL scheme on interior nodes / cells
            this->MUSCL_scheme();
        }

        for (int i = 1; i < this->case_parameters.number_of_points - 1; i++)
        {
            if (std::isnan(U[i][0]) || std::isnan(U[i][1]) || std::isnan(U[i][2]))
            {
                std::cout << "NaN detected at cell " << i << "\n";
                return;
            }
            for (int face = FACE::WEST; face <= FACE::EAST; face++)
            {
                // compute fluxes at faces
                auto [flux_left, flux_right] = this->compute_flux(i, face);
                // Rusanov Riemann solver
                this->Rusanov_Riemann_solver(flux_left, flux_right, i, face);
            }
        }
        // update solution
        this->integrate_solution_in_time();
        // update boundary conditions
        double rho_left = 1.0;
        double u_left = 0.0;
        double rho_right = 0.125;
        double u_right = 0.0;

        this->U[0][0] = U[1][0];
        U[0][1] = rho_left * u_left;
        U[0][2] = U[1][2];
        U[this->case_parameters.number_of_points - 1][0] = this->U[this->case_parameters.number_of_points - 2][0];
        U[this->case_parameters.number_of_points - 1][1] = rho_right * u_right;
        U[this->case_parameters.number_of_points - 1][2] = U[this->case_parameters.number_of_points - 2][2];

        std::cout << "Current time: " << std::scientific << std::setw(10) << std::setprecision(3) << case_parameters.time;
        std::cout << ", End time: " << std::scientific << std::setw(10) << std::setprecision(3) << case_parameters.end_time;
        std::cout << ", Current time step: " << std::fixed << std::setw(7) << case_parameters.time_step;
        std::cout << "\r";

        case_parameters.time += this->dt;
        case_parameters.time_step++;
    }
    std::cout << "SIMULATION FINISHED\n";
}
void Euler::print_solution()
{
    this->solve();

    std::ofstream output_file;

    // convert time step and points into 6 digits string with leading zeros
    std::ostringstream time_step_temp, points_temp;
    time_step_temp << std::setfill('0') << std::setw(6);
    time_step_temp << this->case_parameters.time_step;
    auto time_step = time_step_temp.str();
    points_temp << std::setfill('0') << std::setw(6);
    points_temp << case_parameters.number_of_points;
    auto points = points_temp.str();

    output_file.open("/home/vaibhavvikas/euler_cfd/solution_" + points + "_" + time_step + ".csv");
    output_file << "x,rho,u,p" << std::endl;
    for (int i = 0; i < case_parameters.number_of_points; i++)
    {
        double rho = this->U[i][0];
        double u = U[i][1] / rho;
        double p = (case_parameters.gamma - 1) * (U[i][2] - 0.5 * rho * std::pow(u, 2));
        output_file << this->x[i] << "," << rho << "," << u << "," << p << std::endl;
    }
    output_file.close();
}
int main(int argc, char **argv)
{
    Euler euler("/home/vaibhavvikas/euler_cfd/config/params.yaml");
    euler.print_solution();

    return 0;
}