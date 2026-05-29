# Chapter 2: Elimination, Matrix Algebra, and LU Decomposition

## 🎯 Core Thesis
To solve a linear system of equations $Ax = b$, row operations are applied to transform the system into an upper triangular form $Ux = c$. Gilbert Strang explains that this physical row elimination is algebraically represented as multiplying $A$ by a sequence of lower-triangular **elimination matrices** ($E_{ij}$). 

When these elimination steps are inverted and compiled, they naturally yield the **$A = LU$ factorization**, where $L$ is a lower-triangular matrix containing the row multipliers and $U$ is the final upper-triangular matrix containing the pivots. 

Understanding $A = LU$ is essential because it decouples the matrix operations on $A$ from the vector $b$, allowing rapid solving of multiple systems using forward and back substitution.

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

### 2. Elimination Matrices ($E_{ij}$)
To subtract a multiplier $l_{ij}$ of Row $j$ from Row $i$, we multiply the target matrix by an elimination matrix $E_{ij}$. $E_{ij}$ is equal to the Identity matrix except for a non-zero element $-l_{ij}$ in the $(i, j)$ position.
*   *Example:* To subtract $2$ times Row 1 from Row 2 in a $3 \times 3$ system, we use:
    $$E_{21} = \begin{bmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}$$

### 3. The $A = LU$ Factorization
If no row exchanges are required during elimination, the sequence of row operations yields:
$$(E_{32} E_{31} E_{21}) A = U$$

Inverting this sequence to isolate $A$ yields:
$$A = (E_{21}^{-1} E_{31}^{-1} E_{32}^{-1}) U = L U$$

*   **No Multiplier Interference**: While the product of elimination matrices $(E_{32}E_{31}E_{21})$ contains complex mixed terms, the inverse product $L$ has a beautiful structure: the multipliers $l_{ij}$ go directly into their respective slots below the diagonal without any modification:
    $$L = \begin{bmatrix} 1 & 0 & 0 \\ l_{21} & 1 & 0 \\ l_{31} & l_{32} & 1 \end{bmatrix}$$
*   **Permutations**: If row exchanges are required to bypass a zero pivot, a permutation matrix $P$ is pre-multiplied, yielding:
    $$PA = LU$$

---

## ✍️ Step-by-Step Worked Example: Factoring $A = LU$

### Problem:
Perform a complete $A = LU$ factorization for the following $3 \times 3$ matrix:

$$A = \begin{bmatrix} 2 & 1 & 1 \\ 4 & 5 & 0 \\ 2 & 7 & 8 \end{bmatrix}$$

### Step 1: Eliminate $A_{21}$ (Row 2, Column 1)
The first pivot is $d_1 = 2$. The multiplier to eliminate the value $4$ in the $(2, 1)$ slot is:
$$l_{21} = \frac{4}{2} = 2$$

We subtract $2$ times Row 1 from Row 2 ($R_2 \leftarrow R_2 - 2R_1$):
$$E_{21} A = \begin{bmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} 2 & 1 & 1 \\ 4 & 5 & 0 \\ 2 & 7 & 8 \end{bmatrix} = \begin{bmatrix} 2 & 1 & 1 \\ 0 & 3 & -2 \\ 2 & 7 & 8 \end{bmatrix}$$

### Step 2: Eliminate $A_{31}$ (Row 3, Column 1)
The multiplier to eliminate the value $2$ in the $(3, 1)$ slot is:
$$l_{31} = \frac{2}{2} = 1$$

We subtract Row 1 from Row 3 ($R_3 \leftarrow R_3 - R_1$):
$$E_{31} (E_{21} A) = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -1 & 0 & 1 \end{bmatrix} \begin{bmatrix} 2 & 1 & 1 \\ 0 & 3 & -2 \\ 2 & 7 & 8 \end{bmatrix} = \begin{bmatrix} 2 & 1 & 1 \\ 0 & 3 & -2 \\ 0 & 6 & 7 \end{bmatrix}$$

### Step 3: Eliminate $A_{32}$ (Row 3, Column 2)
The second pivot is $d_2 = 3$. The multiplier to eliminate the value $6$ in the $(3, 2)$ slot is:
$$l_{32} = \frac{6}{3} = 2$$

We subtract $2$ times Row 2 from Row 3 ($R_3 \leftarrow R_3 - 2R_2$):
$$E_{32} (E_{31} E_{21} A) = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -2 & 1 \end{bmatrix} \begin{bmatrix} 2 & 1 & 1 \\ 0 & 3 & -2 \\ 0 & 6 & 7 \end{bmatrix} = \begin{bmatrix} 2 & 1 & 1 \\ 0 & 3 & -2 \\ 0 & 0 & 11 \end{bmatrix} = U$$

Our upper triangular matrix $U$ has three non-zero pivots: $2$, $3$, and $11$.

### Step 4: Construct the Multiplier Matrix $L$
We construct $L$ using our calculated multipliers $l_{21} = 2$, $l_{31} = 1$, and $l_{32} = 2$:

$$L = \begin{bmatrix} 1 & 0 & 0 \\ l_{21} & 1 & 0 \\ l_{31} & l_{32} & 1 \end{bmatrix} = \begin{bmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 1 & 2 & 1 \end{bmatrix}$$

### Step 5: Verification
We verify that $LU = A$:

$$LU = \begin{bmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 1 & 2 & 1 \end{bmatrix} \begin{bmatrix} 2 & 1 & 1 \\ 0 & 3 & -2 \\ 0 & 0 & 11 \end{bmatrix} = \begin{bmatrix} (2) & (1) & (1) \\ (4) & (2+3) & (2-2) \\ (2) & (1+6) & (1-4+11) \end{bmatrix} = \begin{bmatrix} 2 & 1 & 1 \\ 4 & 5 & 0 \\ 2 & 7 & 8 \end{bmatrix} = A$$

The factorization is correct.

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
