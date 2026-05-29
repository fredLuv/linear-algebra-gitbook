# Chapter 6: Singular Value Decomposition (SVD) and PCA Derivations

## 🎯 Core Thesis
While matrix diagonalization ($A = S \Lambda S^{-1}$) requires a square matrix with linearly independent eigenvectors, the **Singular Value Decomposition (SVD)** is the ultimate matrix factorization because it applies to **any matrix of any shape** ($m \times n$, rank $r$). 

Gilbert Strang describes SVD as finding an orthonormal basis in the input space $\mathbb{R}^n$ that maps to an orthonormal basis in the output space $\mathbb{R}^m$. The SVD factors a matrix $A$ into:

$$A = U \Sigma V^T$$

In data science and quant research, SVD serves as the underlying mathematical engine of **Principal Component Analysis (PCA)**, which projects high-dimensional data onto orthogonal axes of maximum variance to isolate risk signals from noise.

---

## 🔑 Key Mathematical Formulations

Let $A$ be an $m \times n$ matrix with rank $r$.

$$\begin{bmatrix} A \end{bmatrix} = \begin{bmatrix} & & \\ & U & \\ & & \end{bmatrix} \begin{bmatrix} \sigma_1 & & \\ & \sigma_2 & \\ & & \ddots \end{bmatrix} \begin{bmatrix} & & \\ & V^T & \\ & & \end{bmatrix}$$

### 1. SVD Component Definitions
*   **$U$ (Left Singular Vectors)**: An $m \times m$ orthogonal matrix ($U^T U = I_m$). Its columns are the orthonormal eigenvectors of the symmetric matrix **$A A^T$**.
*   **$V$ (Right Singular Vectors)**: An $n \times n$ orthogonal matrix ($V^T V = I_n$). Its columns are the orthonormal eigenvectors of the symmetric matrix **$A^T A$**.
*   **$\Sigma$ (Singular Values)**: An $m \times n$ diagonal matrix containing the singular values $\sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_r > 0$. The singular values are the square roots of the non-zero eigenvalues of both $A^T A$ and $A A^T$:
    $$\sigma_i = \sqrt{\lambda_i}$$

### 2. Geometric Mapping View
SVD shows that any linear transformation $A$ can be decomposed into three sequential geometric steps:
1.  **Rotation/Reflection ($V^T$)**: Rotate the input vector to align with the principal directions.
2.  **Stretching ($\Sigma$)**: Stretch the vector along the coordinate axes by factors of $\sigma_i$.
3.  **Second Rotation ($U$)**: Rotate the stretched coordinates into their final output space.

### 3. Low-Rank Matrix Approximation (Eckart-Young-Mirsky Theorem)
Any matrix $A$ can be expressed as a sum of rank-1 matrices:
$$A = \sum_{i=1}^{r} \sigma_i \mathbf{u}_i \mathbf{v}_i^T$$

*   **The Theorem**: If we want the best rank-$k$ approximation of $A$ (denoted $A_k$, where $k < r$) that minimizes the Frobenius norm $\|A - A_k\|_F$, we simply truncate the SVD sum to the first $k$ terms:
    $$A_k = \sum_{i=1}^{k} \sigma_i \mathbf{u}_i \mathbf{v}_i^T$$
    This is the mathematical foundation for data compression and noise filtration.

---

## ✍️ Step-by-Step Worked Example: Singular Value Decomposition

### Problem:
Compute the complete Singular Value Decomposition $A = U \Sigma V^T$ for the following asymmetric matrix:

$$A = \begin{bmatrix} 3 & 0 \\ 4 & 0 \end{bmatrix}$$

---

### Step 1: Compute $A^T A$ and Find its Eigenvalues and Eigenvectors ($V$)
1.  Calculate $A^T A$:
    $$A^T A = \begin{bmatrix} 3 & 4 \\ 0 & 0 \end{bmatrix} \begin{bmatrix} 3 & 0 \\ 4 & 0 \end{bmatrix} = \begin{bmatrix} 9+16 & 0 \\ 0 & 0 \end{bmatrix} = \begin{bmatrix} 25 & 0 \\ 0 & 0 \end{bmatrix}$$
2.  Find eigenvalues of $A^T A$:
    Since it is diagonal, the eigenvalues are:
    $$\lambda_1 = 25, \quad \lambda_2 = 0$$
3.  Compute the singular values ($\Sigma$):
    $$\sigma_1 = \sqrt{25} = 5, \quad \sigma_2 = 0$$
    The singular value matrix is:
    $$\Sigma = \begin{bmatrix} 5 & 0 \\ 0 & 0 \end{bmatrix}$$
4.  Find eigenvectors of $A^T A$ (columns of $V$):
    *   For $\lambda_1 = 25$:
        $$\begin{bmatrix} 25-25 & 0 \\ 0 & 0-25 \end{bmatrix} \begin{bmatrix} v_1 \\ v_2 \end{bmatrix} = \begin{bmatrix} 0 & 0 \\ 0 & -25 \end{bmatrix} \begin{bmatrix} v_1 \\ v_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix} \implies \mathbf{v}_1 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}$$
    *   For $\lambda_2 = 0$:
        $$\begin{bmatrix} 25-0 & 0 \\ 0 & 0-0 \end{bmatrix} \begin{bmatrix} v_1 \\ v_2 \end{bmatrix} = \begin{bmatrix} 25 & 0 \\ 0 & 0 \end{bmatrix} \begin{bmatrix} v_1 \\ v_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix} \implies \mathbf{v}_2 = \begin{bmatrix} 0 \\ 1 \end{bmatrix}$$
    The right singular vector matrix is:
    $$V = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} \implies V^T = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}$$

---

### Step 2: Compute $A A^T$ and Find its Eigenvectors ($U$)
1.  Calculate $A A^T$:
    $$A A^T = \begin{bmatrix} 3 & 0 \\ 4 & 0 \end{bmatrix} \begin{bmatrix} 3 & 4 \\ 0 & 0 \end{bmatrix} = \begin{bmatrix} 9 & 12 \\ 12 & 16 \end{bmatrix}$$
2.  Find eigenvalues of $A A^T$:
    $$\det(A A^T - \lambda I) = (9-\lambda)(16-\lambda) - 144 = \lambda^2 - 25\lambda + 144 - 144 = \lambda(\lambda - 25) = 0$$
    The non-zero eigenvalue is $\lambda_1 = 25$ (exactly matching $A^T A$).
3.  Find the normalized eigenvectors:
    *   For $\lambda_1 = 25$:
        $$\begin{bmatrix} 9-25 & 12 \\ 12 & 16-25 \end{bmatrix} \begin{bmatrix} u_1 \\ u_2 \end{bmatrix} = \begin{bmatrix} -16 & 12 \\ 12 & -9 \end{bmatrix} \begin{bmatrix} u_1 \\ u_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix} \implies 4u_1 - 3u_2 = 0 \implies \mathbf{u}_1 = \frac{1}{5}\begin{bmatrix} 3 \\ 4 \end{bmatrix}$$
    *   For $\lambda_2 = 0$:
        The second eigenvector must be orthogonal to $\mathbf{u}_1$ in $\mathbb{R}^2$:
        $$\mathbf{u}_2 = \frac{1}{5}\begin{bmatrix} -4 \\ 3 \end{bmatrix}$$
    The left singular vector matrix is:
    $$U = \frac{1}{5}\begin{bmatrix} 3 & -4 \\ 4 & 3 \end{bmatrix}$$

---

### Step 3: Verify the Decomposition $A = U \Sigma V^T$
$$U \Sigma V^T = \left(\frac{1}{5}\begin{bmatrix} 3 & -4 \\ 4 & 3 \end{bmatrix}\right) \begin{bmatrix} 5 & 0 \\ 0 & 0 \end{bmatrix} \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}$$

$$U \Sigma V^T = \frac{1}{5} \begin{bmatrix} 15 & 0 \\ 20 & 0 \end{bmatrix} \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} = \begin{bmatrix} 3 & 0 \\ 4 & 0 \end{bmatrix} = A$$

The Singular Value Decomposition is correct.

---

## 📈 Quant & Data Science Bridging

### The Mathematical Derivation of PCA
Principal Component Analysis (PCA) is widely used in data science and quantitative finance to compress highly dimensional features. Here is the mathematical derivation of PCA using the SVD of the covariance matrix:

1.  Let $X \in \mathbb{R}^{m \times n}$ represent a centered data matrix of $m$ observations (e.g. trading days) and $n$ features (e.g. assets), where each column has a mean of zero.
2.  The $n \times n$ **empirical covariance matrix** $\Sigma_C$ is:
    $$\Sigma_C = \frac{1}{m-1} X^T X$$
3.  We apply SVD to the data matrix $X$:
    $$X = U \Sigma V^T$$
4.  Substitute this SVD into the covariance matrix:
    $$\Sigma_C = \frac{1}{m-1} (U \Sigma V^T)^T (U \Sigma V^T) = \frac{1}{m-1} (V \Sigma^T U^T U \Sigma V^T)$$
5.  Since $U$ is orthogonal, $U^T U = I_m$:
    $$\Sigma_C = \frac{1}{m-1} V \Sigma^2 V^T = V \left( \frac{\Sigma^2}{m-1} \right) V^T$$
6.  This is exactly the **Spectral Decomposition** ($Q \Lambda Q^T$) of the covariance matrix $\Sigma_C$:
    *   The columns of $V$ (the right singular vectors of $X$) are the **Principal Components** (or loading vectors).
    *   The eigenvalues $\lambda_i$ of the covariance matrix are proportional to the squared singular values $\sigma_i^2$:
        $$\lambda_i = \frac{\sigma_i^2}{m-1}$$
        These represent the **variance explained** along each principal component axis.

### Quant Application: Term Structure of Interest Rates (Yield Curve PCA)
In quantitative fixed income trading, interest rates across different maturities (1-month, 2-year, 5-year, 10-year, 30-year treasury bonds) are highly correlated. By running PCA on the interest rate covariance matrix, quants compress these movements into three distinct principal components that explain $98\%+$ of the yield curve's variance:

```
Yield Curve Movements:
1. Shift (PC1 - ~90% variance)   : Linear parallel movement across all maturities.
2. Slope (PC2 - ~7% variance)    : Twist where short-term rates rise while long-term rates fall.
3. Curvature (PC3 - ~1% variance): Butterfly bending of the middle of the curve.
```

*   **Trading Execution**: Rather than trading 15 different maturities, fixed-income arbitrage desks trade three synthetic portfolios mapped directly to **Shift, Slope, and Curvature**, hedging out unwanted interest rate exposure while isolating mispriced yield curve anomalies.
