# Chapter 4: Orthogonality, Projections, and Ordinary Least Squares (OLS)

## 🎯 Core Thesis
When a linear system of equations $Ax = b$ has more equations than variables ($m > n$), the system is **overdetermined**, and typically no exact solution exists because $b$ lies outside the column space $C(A)$. 

Gilbert Strang shows that the optimal proxy solution $\hat{x}$ is found by projecting the vector $b$ orthogonally onto the column space $C(A)$. The projection error $\mathbf{e} = b - A\hat{x}$ must be perpendicular to $C(A)$, meaning it lies in the left nullspace $N(A^T)$. 

This geometric requirement leads directly to the **Ordinary Least Squares (OLS) Normal Equation**:

$$A^T A \hat{x} = A^T b$$

By utilizing the **Gram-Schmidt process**, we can transform the columns of $A$ into an orthonormal basis, yielding the **$A = QR$ decomposition**, which serves as the computationally stable foundation for high-performance linear regressions.

---

## 🔑 Key Mathematical Formulations

### 1. Orthogonality of the Subspaces
The four fundamental subspaces are perpendicular in pairs:
*   **In $\mathbb{R}^n$**: The row space is the orthogonal complement of the nullspace:
    $$C(A^T) = N(A)^\perp \quad \text{and} \quad C(A^T) \perp N(A)$$
*   **In $\mathbb{R}^m$**: The column space is the orthogonal complement of the left nullspace:
    $$C(A) = N(A^T)^\perp \quad \text{and} \quad C(A) \perp N(A^T)$$

### 2. Projection Matrices
To project a vector $b$ onto a subspace spanned by the columns of a full-rank matrix $A$:

$$P = A(A^T A)^{-1} A^T$$

*   **Properties of Projection Matrices**:
    1.  **Symmetric**: $P^T = P$
    2.  **Idempotent**: $P^2 = P$ (projecting a second time does not change the position)

### 3. The Algebraic Derivation of OLS
When $Ax = b$ has no solution, we seek $\hat{x}$ that minimizes the squared Euclidean norm of the error vector $\|\mathbf{e}\|^2 = \|b - Ax\|^2$.

```
              b
             /|
            / |
           /  |  e = b - A*x_hat
          /   |
         /    | 
        /     ▼
       └───────► C(A)
         A*x_hat
```

1.  The error vector $\mathbf{e}$ must be orthogonal to every column of $A$:
    $$\mathbf{a}_i^T (b - A\hat{x}) = 0 \quad \text{for all } i$$
2.  In matrix form, this means the error vector must lie in the left nullspace $N(A^T)$:
    $$A^T (b - A\hat{x}) = \mathbf{0}$$
3.  Expanding this yields the **Normal Equation**:
    $$A^T A \hat{x} = A^T b \implies \hat{x} = (A^T A)^{-1} A^T b$$

### 4. Orthonormal Bases and the $A = QR$ Factorization
A set of vectors $\mathbf{q}_1, \mathbf{q}_2, \dots, \mathbf{q}_n$ is **orthonormal** if they are mutually perpendicular and have unit length:
$$\mathbf{q}_i^T \mathbf{q}_j = \begin{cases} 0 & \text{if } i \neq  j \\ 1 & \text{if } i = j \end{cases}$$

If $Q$ is a matrix with orthonormal columns, then $Q^T Q = I$. If we factor $A = QR$ (where $R$ is upper triangular):
$$A^T A = (QR)^T (QR) = R^T Q^T Q R = R^T R$$

The Normal Equation simplifies to:
$$R^T R \hat{x} = (QR)^T b = R^T Q^T b \implies R \hat{x} = Q^T b$$

---

## ✍️ Step-by-Step Worked Example: Ordinary Least Squares (OLS)

### Problem:
Find the best-fitting line $y = C + Dt$ for the following three data points:
$$(t_1, y_1) = (0, 6), \quad (t_2, y_2) = (1, 0), \quad (t_3, y_3) = (2, 0)$$

---

### Step 1: Set up the Overdetermined System $Ax = b$
Substituting the data points into the line equation $C + Dt = y$:
*   Point 1: $C + D(0) = 6$
*   Point 2: $C + D(1) = 0$
*   Point 3: $C + D(2) = 0$

This yields the matrix equation $Ax = b$ (where $x = \begin{bmatrix} C & D \end{bmatrix}^T$):

$$\begin{bmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{bmatrix} \begin{bmatrix} C \\ D \end{bmatrix} = \begin{bmatrix} 6 \\ 0 \\ 0 \end{bmatrix}$$

This system has no exact solution (the columns of $A$ cannot combine to yield $b$). We must solve via OLS.

---

### Step 2: Compute the Normal Equation Components
1.  Compute $A^T A$:
    $$A^T A = \begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \end{bmatrix} \begin{bmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{bmatrix} = \begin{bmatrix} 3 & 3 \\ 3 & 5 \end{bmatrix}$$
2.  Compute $A^T b$:
    $$A^T b = \begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \end{bmatrix} \begin{bmatrix} 6 \\ 0 \\ 0 \end{bmatrix} = \begin{bmatrix} 6 \\ 0 \end{bmatrix}$$

---

### Step 3: Solve $A^T A \hat{x} = A^T b$
We set up the system:

$$\begin{bmatrix} 3 & 3 \\ 3 & 5 \end{bmatrix} \begin{bmatrix} \hat{C} \\ \hat{D} \end{bmatrix} = \begin{bmatrix} 6 \\ 0 \end{bmatrix}$$

Using elimination or matrix inversion:
1.  Subtract Row 1 from Row 2 ($R_2 \leftarrow R_2 - R_1$):
    $$\begin{bmatrix} 3 & 3 \\ 0 & 2 \end{bmatrix} \begin{bmatrix} \hat{C} \\ \hat{D} \end{bmatrix} = \begin{bmatrix} 6 \\ -6 \end{bmatrix}$$
2.  Solve for $\hat{D}$:
    $$2\hat{D} = -6 \implies \hat{D} = -3$$
3.  Solve for $\hat{C}$ via back substitution:
    $$3\hat{C} + 3(-3) = 6 \implies 3\hat{C} - 9 = 6 \implies 3\hat{C} = 15 \implies \hat{C} = 5$$

---

### Step 4: Verification of the Error Vector
The best-fit line is $y = 5 - 3t$. Let's compute the projection $p = A\hat{x}$ and the error vector $\mathbf{e}$:

$$p = A\hat{x} = \begin{bmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{bmatrix} \begin{bmatrix} 5 \\ -3 \end{bmatrix} = \begin{bmatrix} 5 \\ 2 \\ -1 \end{bmatrix}$$

$$\mathbf{e} = b - p = \begin{bmatrix} 6 \\ 0 \\ 0 \end{bmatrix} - \begin{bmatrix} 5 \\ 2 \\ -1 \end{bmatrix} = \begin{bmatrix} 1 \\ -2 \\ 1 \end{bmatrix}$$

We verify that the error vector $\mathbf{e}$ is orthogonal to the column space $C(A)$ by checking $A^T \mathbf{e} = \mathbf{0}$:

$$A^T \mathbf{e} = \begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \end{bmatrix} \begin{bmatrix} 1 \\ -2 \\ 1 \end{bmatrix} = \begin{bmatrix} 1 - 2 + 1 \\ 0 - 2 + 2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$$

The error vector is perfectly orthogonal to the column space.

---

## 📈 Quant & Data Science Bridging

### 1. Bypassing the Multicollinearity Instability Vector
In quantitative research and data science, we estimate regression coefficients $\hat{x}$ using massive feature matrices $A$. 
*   **The Trap**: Standard statistical software calculates $\hat{x} = (A^T A)^{-1} A^T b$. If your features are highly correlated (known as **multicollinearity**), the matrix $A^T A$ becomes near-singular. Its **condition number** spikes, leading to massive numerical instability where small rounding errors in your double-precision floats completely distort your coefficients.

### 2. The Computational Solution: $A = QR$ Decomposition
To secure numerical stability in production algorithmic trading engines, we never compute $(A^T A)^{-1}$ directly. Instead, we compute the **$A = QR$ decomposition** using Householder reflections:

1.  Decompose $A = QR$. Since $Q$ is orthogonal, it preserves Euclidean lengths and does not amplify numerical errors.
2.  Substitute into the OLS formulation:
    $$R \hat{x} = Q^T b$$
3.  Because $R$ is upper triangular, we solve for $\hat{x}$ instantly using **back substitution** without any matrix inversion. This bypasses the inversion bottleneck, guarantees numerical stability, and protects the system from multicollinearity float truncation.
