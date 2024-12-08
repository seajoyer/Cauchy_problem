import numpy as np
import matplotlib.pyplot as plt
from ode_solver import ODESolver

def main():
    def func(x: float, y: float, dy: float) -> float:
        return -(x/(1+x**2))*dy + (1/(1+x**2))*y + (3-2*x+4*x**2)/(1+x**2) * np.exp(-2*x)

    x0 = 0
    y0 = (2, -2)  # (y(0), y'(0))

    def exact_solution(x):
        return np.sqrt(1 + x**2) + np.exp(-2*x)

    h1 = 0.1
    h2 = 0.05
    x_end = 1.0

    solver_h1 = ODESolver(func, x0, y0, h1, x_end)
    solver_h2 = ODESolver(func, x0, y0, h2, x_end)

    # Solve using different methods
    x_euler_h1, y_euler_h1 = solver_h1.euler_method()
    x_rk4_h1, y_rk4_h1 = solver_h1.runge_kutta_4()
    x_adams_h1, y_adams_h1 = solver_h1.adams_3()

    x_rk4_h2, y_rk4_h2 = solver_h2.runge_kutta_4()


    x_exact = np.linspace(0, 1, 1000)
    y_exact = exact_solution(x_exact)

    runge_error = solver_h2.estimate_runge_error(y_rk4_h2, y_rk4_h1)
    print(f"Estimated maximum error (Runge estimation): {runge_error:.2e}")

    plt.figure(figsize=(12, 8))
    plt.plot(x_exact, y_exact, 'k-', label='Exact solution')
    plt.plot(x_euler_h1, y_euler_h1, 'b--', label=f'Euler (h={h1})')
    plt.plot(x_rk4_h1, y_rk4_h1, 'r--', label=f'RK4 (h={h1})')
    plt.plot(x_adams_h1, y_adams_h1, 'g--', label=f'Adams (h={h1})')
    plt.plot(x_rk4_h2, y_rk4_h2, 'r:', label=f'RK4 (h={h2})')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Comparison of Numerical Methods')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
