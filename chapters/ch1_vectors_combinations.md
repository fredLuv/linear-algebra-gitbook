# Chapter 1: Vectors, Linear Combinations, and Vector Spaces

## 🎯 Core Thesis
The starting point of linear algebra is the **vector**. Gilbert Strang emphasizes that linear algebra is built on two primary operations: **vector addition** ($\mathbf{u} + \mathbf{v}$) and **scalar multiplication** ($c\mathbf{u}$). By combining these two operations, we form **linear combinations**:

$$\mathbf{w} = c\mathbf{u} + d\mathbf{v}$$

Geometrically, the set of all possible linear combinations of a group of vectors defines their **Span**. If these vectors reside in a space and satisfy ten strict closure and algebraic properties, they form a **Vector Space**. The inner product (or **dot product**) introduces the geometric concepts of **length, distance, and angle** to vector spaces, providing the mathematical foundation for measuring similarity and projection.

---

## 🔑 Key Mathematical Formulations & Axioms

### 1. Vector Space Axioms
A real vector space $V$ is a set of elements (vectors) closed under addition and scalar multiplication that satisfies the following ten axioms for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and scalars $c, d \in \mathbb{R}$:

1.  **Closure under Addition**: $\mathbf{u} + \mathbf{v} \in V$
2.  **Commutativity of Addition**: $\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}$
3.  **Associativity of Addition**: $(\mathbf{u} + \mathbf{v}) + \mathbf{w} = \mathbf{u} + (\mathbf{v} + \mathbf{w})$
4.  **Additive Identity**: There exists a zero vector $\mathbf{0} \in V$ such that $\mathbf{u} + \mathbf{0} = \mathbf{u}$
5.  **Additive Inverse**: There exists $-\mathbf{u} \in V$ such that $\mathbf{u} + (-\mathbf{u}) = \mathbf{0}$
6.  **Closure under Scalar Multiplication**: $c\mathbf{u} \in V$
7.  **Distributivity of Scalar Addition**: $(c + d)\mathbf{u} = c\mathbf{u} + d\mathbf{u}$
8.  **Distributivity of Vector Addition**: $c(\mathbf{u} + \mathbf{v}) = c\mathbf{u} + c\mathbf{v}$
9.  **Associativity of Scalar Multiplication**: $c(d\mathbf{u}) = (cd)\mathbf{u}$
10. **Scalar Identity**: $1\mathbf{u} = \mathbf{u}$

### 2. The Dot Product, Lengths, and Angles
For two vectors $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$, the inner product (dot product) is defined as:

$$\mathbf{u} \cdot \mathbf{v} = \mathbf{u}^T \mathbf{v} = \sum_{i=1}^{n} u_i v_i$$

*   **Length (Euclidean $L_2$ Norm)**: The length of a vector $\mathbf{u}$, denoted as $\|\mathbf{u}\|$, is:
    $$\|\mathbf{u}\| = \sqrt{\mathbf{u} \cdot \mathbf{u}} = \sqrt{\sum_{i=1}^{n} u_i^2}$$
*   **Angle Between Vectors**: The cosine of the angle $\theta$ between two non-zero vectors is:
    $$\cos\theta = \frac{\mathbf{u} \cdot \mathbf{v}}{\|\mathbf{u}\| \|\mathbf{v}\|}$$
    *   If $\mathbf{u} \cdot \mathbf{v} = 0$, then $\cos\theta = 0$, meaning the vectors are **orthogonal** ($\theta = 90^\circ$).

### 3. Core Inequalities
*   **Cauchy-Schwarz Inequality**: The absolute value of the dot product is bounded by the product of their lengths:
    $$|\mathbf{u} \cdot \mathbf{v}| \le \|\mathbf{u}\| \|\mathbf{v}\|$$
*   **Triangle Inequality**: The length of the sum of two vectors is bounded by the sum of their individual lengths:
    $$\|\mathbf{u} + \mathbf{v}\| \le \|\mathbf{u}\| + \|\mathbf{v}\|$$

---

## ✍️ Step-by-Step Worked Example: Proving Linear Independence

### Problem:
Determine if the following three vectors in $\mathbb{R}^3$ are **linearly independent**:

$$\mathbf{v}_1 = \begin{bmatrix} 1 \\ 2 \\ 1 \end{bmatrix}, \quad \mathbf{v}_2 = \begin{bmatrix} 2 \\ 9 \\ 0 \end{bmatrix}, \quad \mathbf{v}_3 = \begin{bmatrix} 3 \\ 3 \\ 4 \end{bmatrix}$$

### Mathematical Setup:
By definition, vectors are linearly independent if the only scalar combination that yields the zero vector is the trivial solution ($c_1 = c_2 = c_3 = 0$):

$$c_1 \mathbf{v}_1 + c_2 \mathbf{v}_2 + c_3 \mathbf{v}_3 = \mathbf{0}$$

This can be written as a matrix equation $A\mathbf{c} = \mathbf{0}$:

$$\begin{bmatrix} 1 & 2 & 3 \\ 2 & 9 & 3 \\ 1 & 0 & 4 \end{bmatrix} \begin{bmatrix} c_1 \\ c_2 \\ c_3 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}$$

### Step-by-Step Row Elimination:
We apply Gaussian row elimination to the coefficient matrix $A$:

$$\begin{bmatrix} 1 & 2 & 3 \\ 2 & 9 & 3 \\ 1 & 0 & 4 \end{bmatrix}$$

1.  Subtract $2$ times Row 1 from Row 2 ($R_2 \leftarrow R_2 - 2R_1$):
    $$\begin{bmatrix} 1 & 2 & 3 \\ 0 & 5 & -3 \\ 1 & 0 & 4 \end{bmatrix}$$
2.  Subtract Row 1 from Row 3 ($R_3 \leftarrow R_3 - R_1$):
    $$\begin{bmatrix} 1 & 2 & 3 \\ 0 & 5 & -3 \\ 0 & -2 & 1 \end{bmatrix}$$
3.  Eliminate the second entry in Row 3 by adding $\frac{2}{5}$ of Row 2 to Row 3 ($R_3 \leftarrow R_3 + \frac{2}{5}R_2$):
    $$\begin{bmatrix} 1 & 2 & 3 \\ 0 & 5 & -3 \\ 0 & 0 & 1 - \frac{6}{5} \end{bmatrix} = \begin{bmatrix} 1 & 2 & 3 \\ 0 & 5 & -3 \\ 0 & 0 & -1/5 \end{bmatrix}$$

### Analysis of Pivots:
The matrix is now in upper triangular form. The diagonal entries (pivots) are:
*   Pivot 1 = $1 \neq 0$
*   Pivot 2 = $5 \neq 0$
*   Pivot 3 = $-1/5 \neq 0$

Since we have **three non-zero pivots** for a $3 \times 3$ matrix, the matrix $A$ is full rank (invertible). The only solution to $A\mathbf{c} = \mathbf{0}$ is the trivial solution $\mathbf{c} = \mathbf{0}$ ($c_1 = c_2 = c_3 = 0$).

### Conclusion:
The vectors $\mathbf{v}_1$, $\mathbf{v}_2$, and $\mathbf{v}_3$ are **linearly independent**.

---

## 📈 Quant & Data Science Bridging

### 1. Portfolios as Linear Combinations
In quantitative portfolio management, let $\mathbf{r}_1, \mathbf{r}_2, \dots, \mathbf{r}_n \in \mathbb{R}^T$ represent vectors of historical returns for $n$ distinct assets over $T$ time periods. Let $w_1, w_2, \dots, w_n$ represent the portfolio allocation weights. 

The portfolio's historical return vector $\mathbf{r}_p$ is a **linear combination** of the individual asset vectors:

$$\mathbf{r}_p = \sum_{i=1}^{n} w_i \mathbf{r}_i$$

*   The allocation weights are constrained such that $\sum_{i=1}^{n} w_i = 1$. The set of all possible long-only portfolios ($w_i \ge 0$) forms a convex hull within the span of the asset vectors.

### 2. Cosine Similarity in NLP & Recommendation Engines
In data science, text documents or user profiles are represented as high-dimensional sparse vectors of word frequencies or ratings ($\mathbf{u}, \mathbf{v}$). The similarity between two elements is measured by the **inner product** normalized by their Euclidean lengths:

$$\text{Cosine Similarity}(\mathbf{u}, \mathbf{v}) = \frac{\mathbf{u}^T \mathbf{v}}{\|\mathbf{u}\| \|\mathbf{v}\|} = \cos\theta$$

*   **Bounded Metrics**: By the Cauchy-Schwarz inequality, this similarity metric is bounded between $-1$ and $+1$. For non-negative frequency vectors, it lies between $0$ (completely orthogonal, sharing no terms) and $1$ (collinear, identical proportions).
