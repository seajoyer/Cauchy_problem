#include <iostream>
#include <iomanip>
#include <cmath>
#include "ode_solver.hpp"

double exact_solution(double x) {
    return std::sqrt(1 + x*x) + std::exp(-2*x);
}

int main() {
    ode::SecondOrderODE ode = [](double x, double y, double dy) {
        return (3 - 2*x + 4*x*x)*std::exp(-2*x)/(1 + x*x)
               - x/(1 + x*x)*dy
               - 1/(1 + x*x)*y;
    };

    double x0 = 0.0;
    double y0 = 2.0;
    double dy0 = -2.0;

    ode::ODESolver solver(ode, x0, y0, dy0);

    double h_fine = 0.05;
    double h_coarse = 0.1;
    double x_end = 1.0;

    auto solution_rk4_fine = solver.solveRK4(h_fine, x_end);
    auto solution_rk4_coarse = solver.solveRK4(h_coarse, x_end);
    auto solution_euler = solver.solveEuler(h_fine, x_end);
    auto solution_adams = solver.solveAdams3(h_fine, x_end);

    auto runge_errors = ode::ODESolver::calculateRungeError(
        solution_rk4_fine, solution_rk4_coarse, 4
    );

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "x\t\tТочное\t\tRK4\t\tEuler\t\tAdams3\t\tОценка Рунге\n";
    std::cout << std::string(100, '-') << "\n";

    for (size_t i = 0; i < solution_rk4_fine.size(); i += 2) {
        double x = solution_rk4_fine[i].x;
        double exact = exact_solution(x);

        std::cout << x << "\t"
                 << exact << "\t"
                 << solution_rk4_fine[i].y << "\t"
                 << solution_euler[i].y << "\t"
                 << solution_adams[i].y << "\t";

        if (i/2 < runge_errors.size()) {
            std::cout << runge_errors[i/2];
        }
        std::cout << "\n";
    }

    return 0;
}
