# Chapter 2: Elimination, Matrix Algebra, and LU Decomposition

## 🎯 Core Thesis
To solve a linear system of equations $Ax = b$, row operations are applied to transform the system into an upper triangular form $Ux = c$. Gilbert Strang explains that this physical row elimination is algebraically represented as multiplying $A$ by a sequence of lower-triangular **elimination matrices** ($E_{ij}$). 

When these elimination steps are inverted and compiled, they naturally yield the **$A = LU$ factorization**, where $L$ is a lower-triangular matrix containing the row multipliers and $U$ is the final upper-triangular matrix containing the pivots. 

By factoring the pivots out of $U$, we can further refine this to the highly symmetric **$A = LDU$ decomposition**. Understanding these factorizations is essential because they decouple the matrix operations on $A$ from the vector $b$, allowing rapid solving of multiple systems using forward and back substitution.

---

## 🔑 Key Mathematical Formulations

### 1. Matrix Multiplication Views
Gilbert Strang highlights four distinct ways to view the matrix multiplication $C = AB$ (where $A$ is $m \times n$ and $B$ is $n \times p$):

*   **Dot Product View (Standard)**: Each element $C_{ij}$ is the dot product of Row $i$ of $A$ and Column $j$ of $B$:
    $$C_{ij} = \sum_{k=1}^{n} A_{ik} B_{kj}$$
*   **Column View**: Columns of $C$ are linear combinations of the columns of $A$:
    $$\mathbf{c}_j = A \mathbf{b}_j = b_{1j}\mathbf{a}_1 + b_{2j}\mathbf{a}_2 + \dots + b_{nj}\mathbf{a}_n$$
*   **Row View**: Rows of $C$ are linear combinations of the rows of $B$:
    $$\mathbf{r}_{i, C} = \mathbf{a}_{i, A} B = a_{i1}\mathbf{r}_{1, B} + a_{i2}\mathbf{r}_{2, B} + \dots + a_{in}\mathbf{r}_{n, B}$$
*   **Outer Product View**: $C$ is the sum of $n$ rank-1 matrices formed by multiplying columns of $A$ by rows of $B$:
    $$C = \sum_{k=1}^{n} (\text{Col } k \text{ of } A) \times (\text{Row } k \text{ of } B)$$

---

### 2. Elimination & Permutations
*   **Elimination Matrices ($E_{ij}$)**: To subtract a multiplier $l_{ij}$ of Row $j$ from Row $i$, we multiply the target matrix by an elimination matrix $E_{ij}$. $E_{ij}$ is equal to the Identity matrix except for a non-zero element $-l_{ij}$ in the $(i, j)$ position.
*   **Permutation Matrices ($P$)**: If a pivot slot contains a zero, we must swap rows. A Permutation Matrix $P$ is formed by reordering the rows of the Identity matrix.
    *   **Orthonormal Properties**: Permutation matrices are orthogonal. Thus:
        $$P^T P = I \implies P^T = P^{-1}$$
    *   **Determinant**: The determinant of any permutation matrix is either $+1$ (for an even number of row swaps) or $-1$ (for an odd number of row swaps):
        $$\det(P) = \pm 1$$

---

### 3. The $A = LU$ and $A = LDU$ Factorizations
If no row exchanges are required during elimination:
$$(E_{32} E_{31} E_{21}) A = U \implies A = L U$$

*   **The Multiplier Matrix $L$**: The multipliers $l_{ij}$ go directly into their respective slots below the diagonal in $L$:
    $$L = \begin{bmatrix} 1 & 0 & 0 \\ l_{21} & 1 & 0 \\ l_{31} & l_{32} & 1 \end{bmatrix}$$
*   **Extracting Pivots ($A = LDU$)**: In $A = LU$, the upper triangular matrix $U$ has pivots $d_1, d_2, \dots, d_n$ on its diagonal. We can factor $U$ into a diagonal pivot matrix $D$ and a new upper triangular matrix with $1$s on the diagonal (let's call it $U_{new}$):
    $$U = \begin{bmatrix} d_1 & u_{12} & u_{13} \\ 0 & d_2 & u_{23} \\ 0 & 0 & d_3 \end{bmatrix} = \begin{bmatrix} d_1 & 0 & 0 \\ 0 & d_2 & 0 \\ 0 & 0 & d_3 \end{bmatrix} \begin{bmatrix} 1 & u_{12}/d_1 & u_{13}/d_1 \\ 0 & 1 & u_{23}/d_2 \\ 0 & 0 & 1 \end{bmatrix} = D U_{new}$$
    Substituting this back yields the symmetric **$LDU$ factorization**:
    $$A = LDU$$

---

## ✍️ Step-by-Step Worked Example: Factoring $A = LDU$

### Problem:
Given the matrix $A$:

$$A = \begin{bmatrix} 2 & 4 & 2 \\ 4 & 11 & 1 \\ 2 & 1 & 15 \end{bmatrix}$$

1.  Compute the standard $A = LU$ factorization.
2.  Extract the diagonal pivots to construct the symmetric $A = LDU$ decomposition.

---

### Step 1: Row Elimination to find $U$
*   **Eliminate $A_{21}$**: The first pivot is $d_1 = 2$. Multiplier $l_{21} = \frac{4}{2} = 2$.
    Subtract $2$ times Row 1 from Row 2 ($R_2 \leftarrow R_2 - 2R_1$):
    $$\begin{bmatrix} 2 & 4 & 2 \\ 0 & 3 & -3 \\ 2 & 1 & 15 \end{bmatrix}$$
*   **Eliminate $A_{31}$**: Multiplier $l_{31} = \frac{2}{2} = 1$.
    Subtract Row 1 from Row 3 ($R_3 \leftarrow R_3 - R_1$):
    $$\begin{bmatrix} 2 & 4 & 2 \\ 0 & 3 & -3 \\ 0 & -3 & 13 \end{bmatrix}$$
*   **Eliminate $A_{32}$**: The second pivot is $d_2 = 3$. Multiplier $l_{32} = \frac{-3}{3} = -1$.
    Subtract $-1$ times Row 2 from Row 3 ($R_3 \leftarrow R_3 + R_2$):
    $$\begin{bmatrix} 2 & 4 & 2 \\ 0 & 3 & -3 \\ 0 & 0 & 10 \end{bmatrix} = U$$

Our upper triangular matrix $U$ has pivots **$2$, $3$, and $10$** on the diagonal.

---

### Step 2: Construct the Multiplier Matrix $L$
Using the multipliers $l_{21} = 2$, $l_{31} = 1$, and $l_{32} = -1$:

$$L = \begin{bmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 1 & -1 & 1 \end{bmatrix}$$

---

### Step 3: Extract the Diagonal Pivots to Construct $D$ and $U_{new}$
We build the diagonal pivot matrix $D$:

$$D = \begin{bmatrix} 2 & 0 & 0 \\ 0 & 3 & 0 \\ 0 & 0 & 10 \end{bmatrix}$$

To find $U_{new}$, we divide each row of $U$ by its respective pivot:
*   Row 1 of $U$ divided by $2$: $\begin{bmatrix} 1 & 2 & 1 \end{bmatrix}$
*   Row 2 of $U$ divided by $3$: $\begin{bmatrix} 0 & 1 & -1 \end{bmatrix}$
*   Row 3 of $U$ divided by $10$: $\begin{bmatrix} 0 & 0 & 1 \end{bmatrix}$

This yields the normalized upper triangular matrix:

$$U_{new} = \begin{bmatrix} 1 & 2 & 1 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{bmatrix}$$

*   *Note the Symmetry:* Because the original matrix $A$ was symmetric ($A^T = A$), our $LDU$ factorization exhibits beautiful mathematical symmetry: **$U_{new}$ is the exact transpose of $L$** ($U_{new} = L^T$).
    $$A = L D L^T = \begin{bmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 1 & -1 & 1 \end{bmatrix} \begin{bmatrix} 2 & 0 & 0 \\ 0 & 3 & 0 \\ 0 & 0 & 10 \end{bmatrix} \begin{bmatrix} 1 & 2 & 1 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{bmatrix}$$

---

## 📈 Quant & Data Science Bridging

### The Computational Efficiency of Forward & Back Substitution
In both quantitative research and data science, we frequently solve the linear system $Ax = b$ for multiple different vectors $b$ (e.g. evaluating different portfolio targets or updating model weights).

If we solve $Ax = b$ by calculating the explicit inverse $A^{-1}$:
1.  Computing $A^{-1}$ takes approximately **$2n^3$ flops** ($O(n^3)$ complexity).
2.  Multiplying $A^{-1}b$ takes $2n^2$ flops ($O(n^2)$ complexity).
3.  *Warning:* Computing explicit matrix inverses is numerically unstable and highly prone to floating-point truncation errors.

By utilizing **LU Decomposition**, we achieve a massive computational advantage:
1.  We factor $A = LU$ once, which requires $\frac{2}{3}n^3$ flops ($O(n^3)$ complexity).
2.  We solve the system $LUx = b$ by breaking it into two triangular steps:
    *   **Forward Substitution**: Solve $Ly = b$ for $y$ (requires $n^2$ flops, $O(n^2)$ complexity).
    *   **Back Substitution**: Solve $Ux = y$ for $x$ (requires $n^2$ flops, $O(n^2)$ complexity).
3.  When $b$ changes (which happens continuously in real-time trading systems or batch gradient updates), **we do not re-factor the matrix**. We simply re-run the forward and back substitutions in $O(n^2)$ time, bypassing the expensive $O(n^3)$ calculations and securing numerical stability.
