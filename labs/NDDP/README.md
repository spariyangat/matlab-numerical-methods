# Lab 05 — NDDP (Newton’s Divided Difference Interpolation)

Interpolates `y(x_t)` from unevenly spaced data using Newton’s divided differences to a user-specified tolerance.

**Inputs**
- Data file name (string): two columns `[x, y]` sorted by ascending `x`
- Target `x_t` (scalar)
- Convergence `ε_s` (decimal)

**Behavior**
- Builds the divided-differences table incrementally using pointers
- Updates the approximation and estimated error after each added point
- Prints each approximation and the final divided-differences table
- Continues until convergence or until the data is exhausted

**Return value** (row vector)
- Success: `[1, y_hat, achieved_eps]`
- Data exhausted before convergence: `[0, y_hat, achieved_eps]`
