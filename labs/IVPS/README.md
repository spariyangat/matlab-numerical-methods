# ABM4 IVP Solver (RK4-seeded)

Numerically integrates a system of first-order ODEs using **Adams–Bashforth–Moulton (4th order)** with **RK4** used to seed the necessary initial history.

**Inputs**
- Cell array of `dy/dt` function handles (one per state, ordered)
- Initial **augmented state** row vector: `[t0, y1_0, y2_0, …]`
- Step size `h`
- Total steps `N` (RK4 + ABM combined)
- Convergence tolerance `ε_s` (decimal)

**Behavior**
- RK4 computes the minimal initial points required by ABM4
- ABM4 predictor–corrector advances remaining steps, storing only final (converged) states

**Outputs**
- Writes full time history with columns `[t, y1, y2, …]` using writematrix
- Plots each state vs. time in its own figure titled `Var1`, `Var2`, …

