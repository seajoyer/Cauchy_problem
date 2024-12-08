import numpy as np
from typing import Callable, Tuple, List

class ODESolver:
    """Class for solving ordinary differential equations using various numerical methods."""

    def __init__(self,
                 func: Callable[[float, float, float], float],
                 x0: float,
                 y0: Tuple[float, float],
                 h: float,
                 x_end: float):
        """
        Initialize the ODE solver.

        Args:
            func: Function that returns the second derivative y''
            x0: Initial x value
            y0: Tuple of initial values (y(0), y'(0))
            h: Step size
            x_end: End point of the interval
        """
        self.func = func
        self.x0 = x0
        self.y0  = y0[0]
        self.dy0 = y0[1]
        self.h = h
        self.x_end = x_end

    def _convert_to_system(self, x: float, y: List[float]) -> List[float]:
        """
        Convert second-order ODE to system of first-order ODEs.

        Args:
            x: Current x value
            y: List of [y, y'] values

        Returns:
            List of derivatives [y', y'']
        """
        return [y[1], self.func(x, y[0], y[1])]

    def euler_method(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve ODE using Euler's method.

        Returns:
            Tuple of arrays (x values, y values)
        """
        n = int((self.x_end - self.x0) / self.h)
        x = np.linspace(self.x0, self.x_end, n+1)
        y = np.zeros((n+1, 2))
        y[0] = [self.y0, self.dy0]

        for i in range(n):
            derivatives = self._convert_to_system(x[i], y[i])
            y[i+1] = y[i] + np.array(derivatives) * self.h

        return x, y[:, 0]

    def runge_kutta_4(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve ODE using 4th order Runge-Kutta method.

        Returns:
            Tuple of arrays (x values, y values)
        """
        n = int((self.x_end - self.x0) / self.h)
        x = np.linspace(self.x0, self.x_end, n+1)
        y = np.zeros((n+1, 2))
        y[0] = [self.y0, self.dy0]

        for i in range(n):
            k1 = np.array(self._convert_to_system(x[i], y[i]))
            k2 = np.array(self._convert_to_system(x[i] + self.h/2, y[i] + k1*self.h/2))
            k3 = np.array(self._convert_to_system(x[i] + self.h/2, y[i] + k2*self.h/2))
            k4 = np.array(self._convert_to_system(x[i] + self.h, y[i] + k3*self.h))

            y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) * self.h/6

        return x, y[:, 0]

    def adams_3(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve ODE using Adams third-order method.
        First few points are calculated using RK4, then Adams method is applied.

        Returns:
            Tuple of arrays (x values, y values)
        """
        n = int((self.x_end - self.x0) / self.h)
        x = np.linspace(self.x0, self.x_end, n+1)
        y = np.zeros((n+1, 2))
        y[0] = [self.y0, self.dy0]

        # Calculate first 3 points using RK4
        for i in range(3):
            k1 = np.array(self._convert_to_system(x[i], y[i]))
            k2 = np.array(self._convert_to_system(x[i] + self.h/2, y[i] + k1*self.h/2))
            k3 = np.array(self._convert_to_system(x[i] + self.h/2, y[i] + k2*self.h/2))
            k4 = np.array(self._convert_to_system(x[i] + self.h, y[i] + k3*self.h))

            y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) * self.h/6

        # Use Adams method for remaining points
        for i in range(3, n):
            f = np.zeros((4, 2))
            for j in range(4):
                f[j] = self._convert_to_system(x[i-j], y[i-j])

            y[i+1] = y[i] + self.h * (
                23*f[0] - 16*f[1] + 5*f[2]
            ) / 12

        return x, y[:, 0]

    def estimate_runge_error(self, y_h: np.ndarray, y_2h: np.ndarray) -> float:
        """
        Estimate error using Runge's method for RK4.

        Args:
            y_h: Solution with step h
            y_2h: Solution with step 2h

        Returns:
            Maximum estimated error
        """
        # Take every second point from y_h to match y_2h
        y_h = y_h[::2]
        return np.max(np.abs(y_h - y_2h) / 15)
