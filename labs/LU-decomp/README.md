# LU Decomposition (with Partial Pivoting)

Factorizes a square matrix `A` into `L` (unit lower) and `U` (upper) using formula-based LU with partial pivoting. Tracks row swaps via a generic permutation vector `b_gen`.

**Inputs**
- Square, non-singular matrix `A`
- Output filename (string)

**Behavior**
- Computes `L` and `U` via dot-product formulas
- Applies partial pivoting (row swaps only when needed) and updates `b_gen`
- Concatenates results into an `N × (2N+1)` array `[L | U | b_gen]`

**Output**
- Writes the concatenated array to the specified ASCII file with `writematrix`
- No return value; no command-window output
