# Chapter 5: Eigenvalues, Diagonalization, and Symmetric Matrices

## 🎯 Core Thesis
Most matrices act as operations that rotate and stretch vectors. However, for any square matrix $A$, there exist special directions where the matrix's transformation only scales the vector without changing its direction. These vectors are **eigenvectors** ($\mathbf{x}$), and their scaling factors are **eigenvalues** ($\lambda$).

Gilbert Strang explains that by compiling these eigenvectors, we can **diagonalize** the matrix:

$$A = S \Lambda S^{-1}$$

For real symmetric matrices ($A^T = A$), this factorization reaches its pinnacle in the **Spectral Theorem**: the eigenvectors are mutually orthogonal, yielding $A = Q \Lambda Q^T$. Symmetric matrices that are **Positive Definite** are of paramount importance because they define stable energy states, ellipsoid geometries, and physically valid covariance matrices.

---

## 🔑 Key Mathematical Formulations

### 1. Eigenvalue Algebra & Similar Matrices
To find the eigenvalues of an $n \times n$ matrix $A$, we solve the **Characteristic Equation**:

$$\det(A - \lambda I) = 0 \implies \lambda^2 - \text{tr}(A)\lambda + \det(A) = 0 \quad (\text{for a } 2 \times 2 \text{ matrix})$$

*   **Similarity Transformations**: Two matrices $A$ and $B$ are **similar** if there exists an invertible matrix $M$ such that:
    $$B = M^{-1} A M$$
    *   *Theorem:* Similar matrices share the **exact same eigenvalues**.
    *   *Proof:* 
        $$\det(B - \lambda I) = \det(M^{-1} A M - \lambda M^{-1} I M) = \det(M^{-1} (A - \lambda I) M) = \det(M^{-1}) \det(A - \lambda I) \det(M) = \det(A - \lambda I)$$
        Since their characteristic polynomials are identical, their eigenvalues are identical.
*   **Cayley-Hamilton Theorem**: Every square matrix satisfies its own characteristic polynomial. If $p(\lambda) = \det(A - \lambda I) = 0$, then:
    $$p(A) = \mathbf{0}$$

---

### 2. Matrix Diagonalization & Spectral Theorem
If an $n \times n$ matrix $A$ has $n$ linearly independent eigenvectors:

$$A = S \Lambda S^{-1}$$

For a real symmetric matrix $S$ ($S^T = S$), the eigenvectors are mutually orthogonal, yielding the **Spectral Decomposition**:

$$S = Q \Lambda Q^T = \sum_{i=1}^{n} \lambda_i \mathbf{q}_i \mathbf{q}_i^T$$

---

### 3. Quadratic Forms & Ellipsoid Geometry
For a symmetric matrix $A$, the quadratic form is $f(x) = x^T A x$. The level set equation $x^T A x = 1$ defines an **ellipsoid** in $\mathbb{R}^n$ if and only if $A$ is positive definite:

$$\mathbf{x}^T A \mathbf{x} = \lambda_1 (\mathbf{x}^T \mathbf{q}_1)^2 + \lambda_2 (\mathbf{x}^T \mathbf{q}_2)^2 + \dots + \lambda_n (\mathbf{x}^T \mathbf{q}_n)^2 = 1$$

*   **Principal Axes**: The principal axes of the ellipsoid are aligned exactly with the **orthonormal eigenvectors $\mathbf{q}_i$** of $A$.
*   **Semi-Axes Lengths**: The half-lengths of the ellipsoid's axes are inversely proportional to the square root of the eigenvalues:
    $$\text{Length of semi-axis } i = \frac{1}{\sqrt{\lambda_i}}$$
    *   *Intuition:* Large eigenvalues correspond to high "stiffness" or energy, squeezing the ellipsoid tight along that axis (making it short). Small eigenvalues stretch the ellipsoid out.

```
                    q2 (Axis Length = 1/sqrt(lambda_2))
                     ▲
                     │    * * *
                     │  *       *
                     │*           *
                     └───────────────► q1 (Axis Length = 1/sqrt(lambda_1))
                    *               *
                     *             *
                       * * * * * *
```

---

### 4. The Cholesky Decomposition Algorithm ($A = L L^T$)
If a symmetric matrix $A$ is positive definite, it can be factored uniquely into a lower triangular matrix $L$ and its transpose $L^T$ (where $L$ has positive diagonal entries):

$$\begin{bmatrix} A_{11} & A_{12} & A_{13} \\ A_{21} & A_{22} & A_{23} \\ A_{31} & A_{32} & A_{33} \end{bmatrix} = \begin{bmatrix} L_{11} & 0 & 0 \\ L_{21} & L_{22} & 0 \\ L_{31} & L_{32} & L_{33} \end{bmatrix} \begin{bmatrix} L_{11} & L_{21} & L_{31} \\ 0 & L_{22} & L_{32} \\ 0 & 0 & L_{33} \end{bmatrix}$$

The algebraic algorithm solves for the elements of $L$ column-by-column:
1.  **For Column 1**:
    $$L_{11} = \sqrt{A_{11}}$$
    $$L_{i1} = \frac{A_{i1}}{L_{11}} \quad (\text{for } i > 1)$$
2.  **For Column $j$ ($j > 1$)**:
    $$L_{jj} = \sqrt{A_{jj} - \sum_{k=1}^{j-1} L_{jk}^2}$$
    $$L_{ij} = \frac{1}{L_{jj}} \left( A_{ij} - \sum_{k=1}^{j-1} L_{ik} L_{jk} \right) \quad (\text{for } i > j)$$

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
$$\det(A - \lambda I) = \det\begin{bmatrix} 4-\lambda & 2 \\ 2 & 7-\lambda \end{bmatrix} = \lambda^2 - 11\lambda + 24 = 0 \implies (\lambda - 8)(\lambda - 3) = 0$$

$$\lambda_1 = 8, \quad \lambda_2 = 3 \quad (\text{Both are positive!})$$

---

### Step 2: Find the Orthonormal Eigenvectors
*   **For $\lambda_1 = 8$**: Solve $(A - 8I)\mathbf{x} = \mathbf{0}$:
    $$\begin{bmatrix} -4 & 2 \\ 2 & -1 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix} \implies x_2 = 2x_1 \implies \mathbf{q}_1 = \frac{1}{\sqrt{5}} \begin{bmatrix} 1 \\ 2 \end{bmatrix}$$
*   **For $\lambda_2 = 3$**: Solve $(A - 3I)\mathbf{x} = \mathbf{0}$:
    $$\begin{bmatrix} 1 & 2 \\ 2 & 4 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix} \implies x_1 = -2x_2 \implies \mathbf{q}_2 = \frac{1}{\sqrt{5}} \begin{bmatrix} -2 \\ 1 \end{bmatrix}$$

---

### Step 3: Construct the Spectral Decomposition $Q\Lambda Q^T$
$$Q = \frac{1}{\sqrt{5}} \begin{bmatrix} 1 & -2 \\ 2 & 1 \end{bmatrix}, \quad \Lambda = \begin{bmatrix} 8 & 0 \\ 0 & 3 \end{bmatrix}$$

$$Q\Lambda Q^T = \frac{1}{5} \begin{bmatrix} 8 & -6 \\ 16 & 3 \end{bmatrix} \begin{bmatrix} 1 & 2 \\ -2 & 1 \end{bmatrix} = \frac{1}{5} \begin{bmatrix} 20 & 10 \\ 10 & 35 \end{bmatrix} = \begin{bmatrix} 4 & 2 \\ 2 & 7 \end{bmatrix} = A \quad (\text{Verified})$$

---

## 📈 Quant & Data Science Bridging

### 1. Portfolio Variance Demands Positive Definite Covariance Matrices
In quantitative portfolio optimization, let $\mathbf{w} \in \mathbb{R}^n$ represent portfolio allocation weights, and let $\Sigma \in \mathbb{R}^{n \times n}$ represent the historical covariance matrix of asset returns. The overall portfolio variance is:

$$\sigma^2_p = \mathbf{w}^T \Sigma \mathbf{w}$$

*   **The Physical Constraint**: Covariance matrices must be Symmetric Positive Definite to ensure portfolio variance is strictly positive ($\sigma^2_p > 0$) for all non-zero allocations.

### 2. Monte Carlo Correlated Path Generations (Cholesky Engine)
In quantitative finance derivatives trading (e.g. pricing basket options or simulating multi-asset portfolios), we must generate thousands of simulated future asset price paths that reflect their historical correlations.

1.  We compile the historical asset covariance matrix $\Sigma$ (which is Symmetric Positive Definite).
2.  We perform a **Cholesky Decomposition** to find the lower triangular factor $L$:
    $$\Sigma = L L^T$$
3.  Let $\mathbf{z} \sim \mathcal{N}(\mathbf{0}, I)$ represent a vector of independent, uncorrelated standard Gaussian random variables (easily generated by standard pseudo-random algorithms).
4.  We multiply the independent random vector $\mathbf{z}$ by the Cholesky factor $L$:
    $$\mathbf{x} = L \mathbf{z}$$
5.  *The Proof of Correlation:* The covariance of the transformed vector $\mathbf{x}$ matches our target covariance matrix $\Sigma$ exactly:
    $$\text{Cov}(\mathbf{x}) = \mathbb{E}[\mathbf{x}\mathbf{x}^T] = \mathbb{E}[(L\mathbf{z})(L\mathbf{z})^T] = \mathbb{E}[L\mathbf{z}\mathbf{z}^TL^T] = L\mathbb{E}[\mathbf{z}\mathbf{z}^T]L^T = L(I)L^T = L L^T = \Sigma$$
By using the Cholesky factor as a linear transformation, we generate correlated Monte Carlo paths instantly.
