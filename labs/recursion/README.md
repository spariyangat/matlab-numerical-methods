# Recursive Gaussian Elimination

Solves `A x = b` using a **recursive** implementation of Gaussian elimination with **partial pivoting**.

**Inputs**
- Square, non-singular matrix `A`
- Column vector `b` (rows match `A`)

**Behavior**
- Displays the received `A` and `b`
- Each recursion: select pivot (partial pivoting) → eliminate first column to form a reduced `(n−1)×(n−1)` system → recurse on the reduced system → back-substitute to assemble `x`

**Returns**
- Solution vector `x`
