# Chapter 5: Eigenvalues, Diagonalization, and Symmetric Matrices

## 🎯 Core Thesis
Most matrices act as operations that rotate and stretch vectors. However, for any square matrix $A$, there exist special directions where the matrix's transformation only scales the vector without changing its direction. These vectors are **eigenvectors** ($\mathbf{x}$), and their scaling factors are **eigenvalues** ($\lambda$).

Gilbert Strang explains that by compiling these eigenvectors, we can **diagonalize** the matrix:

$$A = S \Lambda S^{-1}$$

For real symmetric matrices ($A^T = A$), this factorization reaches its pinnacle in the **Spectral Theorem**: the eigenvectors are mutually orthogonal, yielding $A = Q \Lambda Q^T$. Symmetric matrices that are **Positive Definite** are of paramount importance because they define stable energy states and physically valid covariance matrices.

---

## 🔑 Key Mathematical Formulations

### 1. The Eigenvalue Equation & Characteristic Polynomial
To find the eigenvalues of an $n \times n$ matrix $A$, we solve the homogeneous equation:

$$(A - \lambda I)\mathbf{x} = \mathbf{0}$$

For a non-trivial solution ($\mathbf{x} \neq \mathbf{0}$), the matrix $(A - \lambda I)$ must be singular, meaning its determinant must be zero. This yields the **Characteristic Equation**:

$$\det(A - \lambda I) = 0$$

*   **Trace Property**: The sum of the eigenvalues equals the trace (sum of diagonal entries) of $A$:
    $$\sum_{i=1}^{n} \lambda_i = \text{tr}(A)$$
*   **Determinant Property**: The product of the eigenvalues equals the determinant of $A$:
    $$\prod_{i=1}^{n} \lambda_i = \det(A)$$

### 2. Matrix Diagonalization
If an $n \times n$ matrix $A$ has $n$ linearly independent eigenvectors, they can be placed into the columns of an eigenvector matrix $S$. Multiplying $AS$ yields:

$$AS = A \begin{bmatrix} \mathbf{x}_1 & \dots & \mathbf{x}_n \end{bmatrix} = \begin{bmatrix} \lambda_1\mathbf{x}_1 & \dots & \lambda_n\mathbf{x}_n \end{bmatrix} = S \Lambda \implies A = S \Lambda S^{-1}$$

*   Here, $\Lambda$ is a diagonal matrix containing the eigenvalues $\lambda_i$.

### 3. The Spectral Theorem for Symmetric Matrices
If $S$ is a real symmetric matrix ($S^T = S$):
1.  All eigenvalues of $S$ are **real numbers**.
2.  The eigenvectors of $S$ are **mutually orthogonal** (or can be chosen to be orthonormal), meaning the eigenvector matrix is orthogonal ($Q^T Q = I$).
3.  This yields the **Spectral Decomposition** (or Principal Axis Theorem):
    $$S = Q \Lambda Q^T = \sum_{i=1}^{n} \lambda_i \mathbf{q}_i \mathbf{q}_i^T$$

### 4. Symmetric Positive Definite (SPD) Matrices
A symmetric matrix $A$ is **Positive Definite** if it satisfies any of these four equivalent mathematical conditions:
1.  **Eigenvalue Test**: All eigenvalues are strictly positive ($\lambda_i > 0$).
2.  **Energy Test**: The quadratic form has positive energy for all non-zero vectors $x \neq \mathbf{0}$:
    $$x^T A x > 0$$
3.  **Pivot Test**: All pivots in Gaussian elimination are strictly positive ($d_i > 0$).
4.  **Determinant Test**: All upper-left subdeterminants (leading principal minors) are positive.

---

## ✍️ Step-by-Step Worked Example: Diagonalizing and Testing Positivity

### Problem:
Given the symmetric matrix $A$:

$$A = \begin{bmatrix} 4 & 2 \\ 2 & 7 \end{bmatrix}$$

1.  Find the eigenvalues and orthonormal eigenvectors of $A$.
2.  Construct the Spectral Decomposition $A = Q \Lambda Q^T$.
3.  Verify if $A$ is Symmetric Positive Definite.

---

### Step 1: Solve the Characteristic Equation
$$\det(A - \lambda I) = \det\begin{bmatrix} 4-\lambda & 2 \\ 2 & 7-\lambda \end{bmatrix} = 0$$

$$\det(A - \lambda I) = (4-\lambda)(7-\lambda) - 4 = \lambda^2 - 11\lambda + 28 - 4 = \lambda^2 - 11\lambda + 24 = 0$$

Factoring the quadratic equation:
$$(\lambda - 8)(\lambda - 3) = 0 \implies \lambda_1 = 8, \quad \lambda_2 = 3$$

*   *Verification:*
    *   $\text{tr}(A) = 4 + 7 = 11$; $\lambda_1 + \lambda_2 = 8 + 3 = 11$ (Correct).
    *   $\det(A) = (4)(7) - (2)(2) = 24$; $\lambda_1 \lambda_2 = (8)(3) = 24$ (Correct).

---

### Step 2: Find the Eigenvectors
*   **For $\lambda_1 = 8$**: Solve $(A - 8I)\mathbf{x} = \mathbf{0}$:
    $$\begin{bmatrix} 4-8 & 2 \\ 2 & 7-8 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} -4 & 2 \\ 2 & -1 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$$
    This yields the equation $2x_1 - x_2 = 0 \implies x_2 = 2x_1$. A basic eigenvector is:
    $$\mathbf{x}_1 = \begin{bmatrix} 1 \\ 2 \end{bmatrix}$$
    Normalizing to unit length ($L_2$ norm):
    $$\mathbf{q}_1 = \frac{1}{\sqrt{1^2 + 2^2}} \begin{bmatrix} 1 \\ 2 \end{bmatrix} = \frac{1}{\sqrt{5}} \begin{bmatrix} 1 \\ 2 \end{bmatrix}$$

*   **For $\lambda_2 = 3$**: Solve $(A - 3I)\mathbf{x} = \mathbf{0}$:
    $$\begin{bmatrix} 4-3 & 2 \\ 2 & 7-3 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 1 & 2 \\ 2 & 4 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$$
    This yields the equation $x_1 + 2x_2 = 0 \implies x_1 = -2x_2$. A basic eigenvector is:
    $$\mathbf{x}_2 = \begin{bmatrix} -2 \\ 1 \end{bmatrix}$$
    Normalizing to unit length:
    $$\mathbf{q}_2 = \frac{1}{\sqrt{(-2)^2 + 1^2}} \begin{bmatrix} -2 \\ 1 \end{bmatrix} = \frac{1}{\sqrt{5}} \begin{bmatrix} -2 \\ 1 \end{bmatrix}$$

*   *Verification of Orthogonality*: Check $\mathbf{q}_1^T \mathbf{q}_2$:
    $$\mathbf{q}_1^T \mathbf{q}_2 = \left(\frac{1}{\sqrt{5}} \begin{bmatrix} 1 \\ 2 \end{bmatrix}\right)^T \left(\frac{1}{\sqrt{5}} \begin{bmatrix} -2 \\ 1 \end{bmatrix}\right) = \frac{1}{5} ( (1)(-2) + (2)(1) ) = 0 \quad (\text{Orthonormal!})$$

---

### Step 3: Construct the Spectral Decomposition $Q\Lambda Q^T$
$$Q = \frac{1}{\sqrt{5}} \begin{bmatrix} 1 & -2 \\ 2 & 1 \end{bmatrix}, \quad \Lambda = \begin{bmatrix} 8 & 0 \\ 0 & 3 \end{bmatrix}$$

$$Q\Lambda Q^T = \left(\frac{1}{\sqrt{5}} \begin{bmatrix} 1 & -2 \\ 2 & 1 \end{bmatrix}\right) \begin{bmatrix} 8 & 0 \\ 0 & 3 \end{bmatrix} \left(\frac{1}{\sqrt{5}} \begin{bmatrix} 1 & 2 \\ -2 & 1 \end{bmatrix}\right)$$

$$Q\Lambda Q^T = \frac{1}{5} \begin{bmatrix} 8 & -6 \\ 16 & 3 \end{bmatrix} \begin{bmatrix} 1 & 2 \\ -2 & 1 \end{bmatrix} = \frac{1}{5} \begin{bmatrix} 8+12 & 16-6 \\ 16-6 & 32+3 \end{bmatrix} = \frac{1}{5} \begin{bmatrix} 20 & 10 \\ 10 & 35 \end{bmatrix} = \begin{bmatrix} 4 & 2 \\ 2 & 7 \end{bmatrix} = A$$

The Spectral Decomposition is correct.

---

### Step 4: Verify if $A$ is Symmetric Positive Definite
We execute the four fundamental tests:
1.  **Eigenvalue Test**: The eigenvalues are $\lambda_1 = 8 > 0$ and $\lambda_2 = 3 > 0$. (Passes).
2.  **Determinant Test**: 
    *   Upper-left $1 \times 1$ determinant = $4 > 0$.
    *   Upper-left $2 \times 2$ determinant = $(4)(7) - (2)(2) = 24 > 0$. (Passes).
3.  **Pivot Test**: Under row elimination, the pivots are:
    *   Pivot 1 = $4 > 0$.
    *   Pivot 2 = $7 - \frac{2}{4}(2) = 6 > 0$. (Passes).

Matrix $A$ is **Symmetric Positive Definite (SPD)**.

---

## 📈 Quant & Data Science Bridging

### Why Portfolio Variance Demands Positive Definite Covariance Matrices
In quantitative portfolio optimization (Modern Portfolio Theory), let $\mathbf{w} \in \mathbb{R}^n$ represent the portfolio allocation weights vector, and let $\Sigma \in \mathbb{R}^{n \times n}$ represent the historical covariance matrix of asset returns.

The overall portfolio return variance $\sigma^2_p$ is calculated as a **quadratic form**:

$$\sigma^2_p = \mathbf{w}^T \Sigma \mathbf{w}$$

*   **The Physical Constraint**: By definition, a physical variance $\sigma^2_p$ **cannot be negative**. To guarantee that any arbitrary long/short portfolio allocation $\mathbf{w}$ has a strictly positive variance ($\sigma^2_p > 0$ for all $\mathbf{w} \neq \mathbf{0}$), the covariance matrix $\Sigma$ must be **Symmetric Positive Definite**.
*   **The Trap**: If your historical return series has fewer historical observations than the number of assets ($T < n$), or contains duplicate assets, the covariance matrix $\Sigma$ becomes **positive semi-definite** (having a zero eigenvalue). This means there exists a non-zero portfolio allocation $\mathbf{w}$ in the nullspace of $\Sigma$ where $\mathbf{w}^T \Sigma \mathbf{w} = 0$, implying you can trade a portfolio with **infinite leverage and absolute zero variance**—a physical impossibility that breaks optimizer math.
*   **The Safeguard**: Quants use **eigenvalue shrinkage** or regularization (like Ridge regression or Ledoit-Wolf estimation) to artificially shift the eigenvalues of the historical covariance matrix away from zero, guaranteeing it is strictly positive definite.
