#pragma once

#include <functional>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace ode {

struct State {
    double x;
    double y;
    double dy;

    State(double x = 0.0, double y = 0.0, double dy = 0.0)
        : x(x), y(y), dy(dy) {}
};

using SecondOrderODE = std::function<double(double x, double y, double dy)>;

class ODESolver {
public:
    ODESolver(SecondOrderODE ode, double x0, double y0, double dy0)
        : m_ode(ode), m_initial_state(x0, y0, dy0) {}

    std::vector<State> solveRK4(double h, double x_end) {
        std::vector<State> solution;
        solution.reserve(static_cast<size_t>((x_end - m_initial_state.x) / h) + 1);

        State current = m_initial_state;
        solution.push_back(current);

        while (current.x < x_end - h/2) {
            current = rk4Step(current, h);
            solution.push_back(current);
        }

        return solution;
    }

    std::vector<State> solveEuler(double h, double x_end) {
        std::vector<State> solution;
        solution.reserve(static_cast<size_t>((x_end - m_initial_state.x) / h) + 1);

        State current = m_initial_state;
        solution.push_back(current);

        while (current.x < x_end - h/2) {
            current = eulerStep(current, h);
            solution.push_back(current);
        }

        return solution;
    }

    std::vector<State> solveAdams3(double h, double x_end) {
        std::vector<State> solution;
        solution.reserve(static_cast<size_t>((x_end - m_initial_state.x) / h) + 1);

        for (int i = 0; i < 3; ++i) {
            if (i == 0) {
                solution.push_back(m_initial_state);
            } else {
                solution.push_back(rk4Step(solution.back(), h));
            }
        }

        while (solution.back().x < x_end - h/2) {
            solution.push_back(adams3Step(solution, h));
        }

        return solution;
    }

    static std::vector<double> calculateRungeError(
        const std::vector<State>& fine,
        const std::vector<State>& coarse,
        double p
    ) {
        std::vector<double> errors;
        size_t coarse_idx = 0;

        for (size_t i = 0; i < fine.size(); i += 2) {
            if (coarse_idx >= coarse.size()) break;

            double error = std::abs(fine[i].y - coarse[coarse_idx].y) /
                          (1 - std::pow(2.0, -p));
            errors.push_back(error);
            coarse_idx++;
        }

        return errors;
    }

private:
    SecondOrderODE m_ode;
    State m_initial_state;

    std::pair<double, double> system(const State& state) {
        return {
            state.dy,
            m_ode(state.x, state.y, state.dy)
        };
    }

    State rk4Step(const State& current, double h) {
        auto [k1_dy, k1_d2y] = system(current);

        State k2_state{
            current.x + h/2,
            current.y + h/2 * k1_dy,
            current.dy + h/2 * k1_d2y
        };
        auto [k2_dy, k2_d2y] = system(k2_state);

        State k3_state{
            current.x + h/2,
            current.y + h/2 * k2_dy,
            current.dy + h/2 * k2_d2y
        };
        auto [k3_dy, k3_d2y] = system(k3_state);

        State k4_state{
            current.x + h,
            current.y + h * k3_dy,
            current.dy + h * k3_d2y
        };
        auto [k4_dy, k4_d2y] = system(k4_state);

        return State{
            current.x + h,
            current.y + h/6 * (k1_dy + 2*k2_dy + 2*k3_dy + k4_dy),
            current.dy + h/6 * (k1_d2y + 2*k2_d2y + 2*k3_d2y + k4_d2y)
        };
    }

    // Single step of Euler method
    State eulerStep(const State& current, double h) {
        auto [dy, d2y] = system(current);
        return State{
            current.x + h,
            current.y + h * dy,
            current.dy + h * d2y
        };
    }

    // Single step of Adams third-order method
    State adams3Step(const std::vector<State>& prev_states, double h) {
        const auto& current = prev_states.back();
        const auto& prev1 = prev_states[prev_states.size() - 2];
        const auto& prev2 = prev_states[prev_states.size() - 3];

        auto [f0_dy, f0_d2y] = system(current);
        auto [f1_dy, f1_d2y] = system(prev1);
        auto [f2_dy, f2_d2y] = system(prev2);

        return State{
            current.x + h,
            current.y + h * (23*f0_dy - 16*f1_dy + 5*f2_dy) / 12,
            current.dy + h * (23*f0_d2y - 16*f1_d2y + 5*f2_d2y) / 12
        };
    }
};

}
