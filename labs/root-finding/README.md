# Lab 04 — Root Finding

MATLAB implementations of four nonlinear solvers:

- **False Position (Regula Falsi):** Bracketed method using `[xL, xU]`; interpolates a secant line to update the bracketing root.
- **Secant:** Derivative-free method using two initial guesses; updates via successive secant lines.
- **Müller:** Uses a quadratic fit through three points to predict the root; works with three starting values.
- **Modified Secant (Multi-D):** Solves a system by approximating the Jacobian with finite differences and updating all variables simultaneously.
