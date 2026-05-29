# Chapter 6: Singular Value Decomposition (SVD) and PCA Derivations

## 🎯 Core Thesis
While matrix diagonalization ($A = S \Lambda S^{-1}$) requires a square matrix with linearly independent eigenvectors, the **Singular Value Decomposition (SVD)** is the ultimate matrix factorization because it applies to **any matrix of any shape** ($m \times n$, rank $r$). 

Gilbert Strang describes SVD as finding an orthonormal basis in the input space $\mathbb{R}^n$ that maps to an orthonormal basis in the output space $\mathbb{R}^m$. The SVD factors a matrix $A$ into:

$$A = U \Sigma V^T$$

In data science and quant research, SVD serves as the underlying mathematical engine of **Principal Component Analysis (PCA)**, which projects high-dimensional data onto orthogonal axes of maximum variance to isolate risk signals from noise.

---

## 🔑 Key Mathematical Formulations

Let $A$ be an $m \times n$ matrix with rank $r$.

### 1. SVD Component Definitions
*   **$U$ (Left Singular Vectors)**: An $m \times m$ orthogonal matrix ($U^T U = I_m$). Its columns are the orthonormal eigenvectors of the symmetric matrix **$A A^T$**.
*   **$V$ (Right Singular Vectors)**: An $n \times n$ orthogonal matrix ($V^T V = I_n$). Its columns are the orthonormal eigenvectors of the symmetric matrix **$A^T A$**.
*   **$\Sigma$ (Singular Values)**: An $m \times n$ diagonal matrix containing the singular values $\sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_r > 0$. The singular values are the square roots of the non-zero eigenvalues of both $A^T A$ and $A A^T$:

$$\sigma_i = \sqrt{\lambda_i}$$

---

### 2. Advanced Decompositions & Operators

#### The Moore-Penrose Pseudoinverse ($A^\dagger$)
If a matrix $A$ is rectangular or singular (not invertible), we can construct its unique pseudoinverse $A^\dagger$ (size $n \times m$) using its SVD components:

$$A^\dagger = V \Sigma^\dagger U^T$$

*   **Constructing $\Sigma^\dagger$**: Take the transpose of $\Sigma$ and replace every non-zero singular value $\sigma_i$ with its reciprocal $\frac{1}{\sigma_i}$:

$$\Sigma^\dagger = \begin{bmatrix} 1/\sigma_1 & 0 & 0 \\ 0 & 1/\sigma_2 & 0 \\ 0 & 0 & 0 \end{bmatrix}$$

*   **Optimal Solver**: For any linear system $Ax = b$, the vector $\hat{x} = A^\dagger b$ is the absolute optimal solution that:
    1.  Minimizes the residual error $\|Ax - b\|_2$ (the Least Squares solution).
    2.  If the system has infinite solutions, $\hat{x}$ selects the unique solution with the **minimum Euclidean norm** $\|x\|_2$.

#### Polar Decomposition ($A = QS$)
Analogous to representing a complex number in polar coordinates ($z = re^{i\theta}$), any real square matrix $A$ can be factored into:

$$A = Q S$$

*   **$Q$ (Orthogonal)**: An orthogonal matrix ($Q^T Q = I$) representing pure rotation and reflection.
*   **$S$ (Symmetric Positive Semi-Definite)**: A symmetric matrix ($S^T = S$) representing pure stretching along orthogonal directions.
*   *Derivation from SVD:*

$$A = U \Sigma V^T = (U V^T)(V \Sigma V^T) = Q S$$

    Where $Q = U V^T$ and $S = V \Sigma V^T$.

---

### 3. Low-Rank Matrix Approximation (Eckart-Young-Mirsky Theorem)
Any matrix $A$ can be expressed as a sum of rank-1 matrices:
$$A = \sum_{i=1}^{r} \sigma_i \mathbf{u}_i \mathbf{v}_i^T$$

*   **The Theorem**: If we want the best rank-$k$ approximation of $A$ (denoted $A_k$, where $k < r$) that minimizes the Frobenius norm $\|A - A_k\|_F$, we simply truncate the SVD sum to the first $k$ terms:

$$A_k = \sum_{i=1}^{k} \sigma_i \mathbf{u}_i \mathbf{v}_i^T$$

---

## ✍️ Step-by-Step Worked Example: Singular Value Decomposition

### Problem:
Compute the complete Singular Value Decomposition $A = U \Sigma V^T$ for the following asymmetric matrix:

$$A = \begin{bmatrix} 3 & 0 \\ 4 & 0 \end{bmatrix}$$

---

### Step 1: Compute $A^T A$ and Find its Eigenvalues and Eigenvectors ($V$)
1.  Calculate $A^T A$:

$$A^T A = \begin{bmatrix} 3 & 4 \\ 0 & 0 \end{bmatrix} \begin{bmatrix} 3 & 0 \\ 4 & 0 \end{bmatrix} = \begin{bmatrix} 25 & 0 \\ 0 & 0 \end{bmatrix}$$

2.  Find eigenvalues of $A^T A$:

$$\lambda_1 = 25, \quad \lambda_2 = 0 \implies \sigma_1 = 5, \quad \sigma_2 = 0 \implies \Sigma = \begin{bmatrix} 5 & 0 \\ 0 & 0 \end{bmatrix}$$

3.  Right singular vectors (columns of $V$):

$$\mathbf{v}_1 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}, \quad \mathbf{v}_2 = \begin{bmatrix} 0 \\ 1 \end{bmatrix} \implies V = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}$$

---

### Step 2: Compute $A A^T$ and Find its Eigenvectors ($U$)
1.  Calculate $A A^T$:

$$A A^T = \begin{bmatrix} 3 & 0 \\ 4 & 0 \end{bmatrix} \begin{bmatrix} 3 & 4 \\ 0 & 0 \end{bmatrix} = \begin{bmatrix} 9 & 12 \\ 12 & 16 \end{bmatrix}$$

2.  Find eigenvalues of $A A^T$:

$$\det(A A^T - \lambda I) = \lambda(\lambda - 25) = 0 \implies \lambda_1 = 25$$

3.  Find the normalized eigenvectors:
    *   For $\lambda_1 = 25$:

$$\begin{bmatrix} -16 & 12 \\ 12 & -9 \end{bmatrix} \begin{bmatrix} u_1 \\ u_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix} \implies 4u_1 - 3u_2 = 0 \implies \mathbf{u}_1 = \frac{1}{5}\begin{bmatrix} 3 \\ 4 \end{bmatrix}$$

    *   For $\lambda_2 = 0$:

$$\mathbf{u}_2 = \frac{1}{5}\begin{bmatrix} -4 \\ 3 \end{bmatrix}$$

    The left singular vector matrix is:

$$U = \frac{1}{5}\begin{bmatrix} 3 & -4 \\ 4 & 3 \end{bmatrix}$$

---

### Step 3: Verify the Decomposition $A = U \Sigma V^T$
$$U \Sigma V^T = \left(\frac{1}{5}\begin{bmatrix} 3 & -4 \\ 4 & 3 \end{bmatrix}\right) \begin{bmatrix} 5 & 0 \\ 0 & 0 \end{bmatrix} \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} = \begin{bmatrix} 3 & 0 \\ 4 & 0 \end{bmatrix} = A \quad (\text{Verified})$$

---

## 📈 Quant & Data Science Bridging

### 1. The Mathematical Derivation of PCA
Principal Component Analysis (PCA) is widely used in data science and quantitative finance to compress highly dimensional features. Here is the mathematical derivation of PCA using the SVD of the covariance matrix:

#### Data Pre-processing (Centering and Scaling):
Let $X_{raw} \in \mathbb{R}^{m \times n}$ represent a raw data matrix of $m$ observations and $n$ features. We **must** pre-process the matrix to prevent features with large arbitrary scales from dominating the principal components:
1.  **Centering**: Subtract the mean vector $\boldsymbol{\mu} \in \mathbb{R}^n$ from each row:

$$X_{centered} = X_{raw} - \mathbf{1}\boldsymbol{\mu}^T \quad \text{where} \quad \mu_j = \frac{1}{m}\sum_{i=1}^{m} X_{raw, ij}$$

2.  **Scaling**: Divide each column by its standard deviation $\sigma_j$ to yield unit variance:

$$X_{ij} = \frac{X_{centered, ij}}{\sigma_j} \quad \text{where} \quad \sigma_j = \sqrt{\frac{1}{m-1}\sum_{i=1}^{m} (X_{raw, ij} - \mu_j)^2}$$

#### PCA Covariance Derivations:
1.  The $n \times n$ **empirical covariance matrix** $\Sigma_C$ of the centered and scaled matrix $X$ is:

$$\Sigma_C = \frac{1}{m-1} X^T X$$

2.  We apply SVD to the data matrix $X$:

$$X = U \Sigma V^T$$

3.  Substitute this SVD into the covariance matrix:

$$\Sigma_C = \frac{1}{m-1} (U \Sigma V^T)^T (U \Sigma V^T) = \frac{1}{m-1} (V \Sigma^T U^T U \Sigma V^T)$$

4.  Since $U$ is orthogonal, $U^T U = I_m$:

$$\Sigma_C = \frac{1}{m-1} V \Sigma^2 V^T = V \left( \frac{\Sigma^2}{m-1} \right) V^T$$

5.  This is exactly the **Spectral Decomposition** ($Q \Lambda Q^T$) of the covariance matrix $\Sigma_C$:
    *   The columns of $V$ (the right singular vectors of $X$) are the **Principal Components** (or loading vectors).
    *   The eigenvalues $\lambda_i$ of the covariance matrix are proportional to the squared singular values $\sigma_i^2$:

$$\lambda_i = \frac{\sigma_i^2}{m-1}$$

        These represent the **variance explained** along each principal component axis.

### 2. Quant Application: Term Structure of Interest Rates (Yield Curve PCA)
In quantitative fixed income trading, interest rates across different maturities are highly correlated. By running PCA on the interest rate covariance matrix, quants compress these movements into three distinct principal components that explain $98\%+$ of the yield curve's variance:

*   **PC1 (Shift - ~90% variance)**: Linear parallel movement across all maturities.
*   **PC2 (Slope - ~7% variance)**: Twist where short-term rates rise while long-term rates fall.
*   **PC3 (Curvature - ~1% variance)**: Butterfly bending of the middle of the curve.
*   **Trading Execution**: Rather than trading 15 different maturities, fixed-income arbitrage desks trade three synthetic portfolios mapped directly to **Shift, Slope, and Curvature**, hedging out unwanted interest rate exposure while isolating mispriced yield curve anomalies.