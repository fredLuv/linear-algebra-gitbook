# Chapter 1: Vectors, Linear Combinations, and Vector Spaces

## 🎯 Core Thesis
The starting point of linear algebra is the **vector**. Gilbert Strang emphasizes that linear algebra is built on two primary operations: **vector addition** ($\mathbf{u} + \mathbf{v}$) and **scalar multiplication** ($c\mathbf{u}$). By combining these two operations, we form **linear combinations**:

$$\mathbf{w} = c\mathbf{u} + d\mathbf{v}$$

Geometrically, the set of all possible linear combinations of a group of vectors defines their **Span**. If these vectors reside in a space and satisfy ten strict closure and algebraic properties, they form a **Vector Space**. The inner product (or **dot product**) introduces the geometric concepts of **length, distance, and angle** to vector spaces, providing the mathematical foundation for measuring similarity, projection, and distance metrics.

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

---

### 2. Inner Product Spaces
An inner product space is a vector space $V$ equipped with a map $\langle \mathbf{u}, \mathbf{v} \rangle : V \times V \to \mathbb{R}$ that satisfies the following four axioms for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and scalars $c \in \mathbb{R}$:

1.  **Symmetry**: $\langle \mathbf{u}, \mathbf{v} \rangle = \langle \mathbf{v}, \mathbf{u} \rangle$
2.  **Linearity**: $\langle c\mathbf{u} + \mathbf{v}, \mathbf{w} \rangle = c\langle \mathbf{u}, \mathbf{w} \rangle + \langle \mathbf{v}, \mathbf{w} \rangle$
3.  **Positive Definiteness**: $\langle \mathbf{u}, \mathbf{u} \rangle \ge 0$, and $\langle \mathbf{u}, \mathbf{u} \rangle = 0 \iff \mathbf{u} = \mathbf{0}$

*   **The Standard Euclidean Dot Product**: In $\mathbb{R}^n$, the standard inner product is:

$$\mathbf{u} \cdot \mathbf{v} = \mathbf{u}^T \mathbf{v} = \sum_{i=1}^{n} u_i v_i$$

*   **Angle Between Vectors**: The cosine of the angle $\theta$ between two non-zero vectors is:

$$\cos\theta = \frac{\langle \mathbf{u}, \mathbf{v} \rangle}{\|\mathbf{u}\| \|\mathbf{v}\|}$$

    *   If $\langle \mathbf{u}, \mathbf{v} \rangle = 0$, then $\cos\theta = 0$, meaning the vectors are **orthogonal** ($\theta = 90^\circ$).

---

### 3. Vector Norms: $L_1$ and $L_2$ Space
A norm $\|\mathbf{x}\|$ is a function that assigns a strictly positive length to a vector. In quantitative and data science models, we focus on the $L_p$ family of norms:

$$\|\mathbf{x}\|_p = \left( \sum_{i=1}^{n} |x_i|^p \right)^{1/p}$$

*   **$L_2$ Norm (Euclidean Norm)**: Represents the direct straight-line distance:

$$\|\mathbf{x}\|_2 = \sqrt{\mathbf{x}^T \mathbf{x}} = \sqrt{\sum_{i=1}^{n} x_i^2}$$

*   **$L_1$ Norm (Taxicab / Manhattan Norm)**: Represents the sum of absolute values:

$$\|\mathbf{x}\|_1 = \sum_{i=1}^{n} |x_i|$$

*   **$L_\infty$ Norm (Max Norm)**: Represents the maximum absolute component:

$$\|\mathbf{x}\|_\infty = \max_i |x_i|$$

---

## ✍️ Step-by-Step Worked Proof: The Cauchy-Schwarz Inequality

### Theorem:
For any vectors $\mathbf{u}, \mathbf{v}$ in an inner product space:

$$|\langle \mathbf{u}, \mathbf{v} \rangle| \le \|\mathbf{u}\| \|\mathbf{v}\|$$

### Step-by-Step Algebraic Proof:
1.  If $\mathbf{v} = \mathbf{0}$, then both sides equal $0$, and the inequality holds trivially. Assume $\mathbf{v} \neq \mathbf{0}$.
2.  Consider a scalar $t \in \mathbb{R}$. By the positive definiteness axiom of inner products, the length of the vector $t\mathbf{u} + \mathbf{v}$ squared must be non-negative for any value of $t$:

$$\langle t\mathbf{u} + \mathbf{v}, t\mathbf{u} + \mathbf{v} \rangle \ge 0$$

3.  Expand this using the linearity and symmetry axioms:

$$\langle t\mathbf{u}, t\mathbf{u} \rangle + 2\langle t\mathbf{u}, \mathbf{v} \rangle + \langle \mathbf{v}, \mathbf{v} \rangle \ge 0$$

$$t^2 \langle \mathbf{u}, \mathbf{u} \rangle + 2t \langle \mathbf{u}, \mathbf{v} \rangle + \langle \mathbf{v}, \mathbf{v} \rangle \ge 0$$

4.  This is a quadratic polynomial in $t$ of the form $a t^2 + b t + c \ge 0$, where:

$$a = \|\mathbf{u}\|^2, \quad b = 2\langle \mathbf{u}, \mathbf{v} \rangle, \quad c = \|\mathbf{v}\|^2$$

5.  For a quadratic polynomial $at^2 + bt + c$ to remain non-negative ($\ge 0$) for all real values of $t$, it must have at most one real root. Mathematically, its **discriminant** must be less than or equal to zero:

$$\Delta = b^2 - 4ac \le 0$$

6.  Substitute the expressions for $a$, $b$, and $c$:

$$(2\langle \mathbf{u}, \mathbf{v} \rangle)^2 - 4 \|\mathbf{u}\|^2 \|\mathbf{v}\|^2 \le 0$$

$$4\langle \mathbf{u}, \mathbf{v} \rangle^2 - 4 \|\mathbf{u}\|^2 \|\mathbf{v}\|^2 \le 0$$

$$\langle \mathbf{u}, \mathbf{v} \rangle^2 \le \|\mathbf{u}\|^2 \|\mathbf{v}\|^2$$

7.  Taking the square root of both sides yields the final inequality:

$$|\langle \mathbf{u}, \mathbf{v} \rangle| \le \|\mathbf{u}\| \|\mathbf{v}\|$$

The theorem is proved.

---

## 📈 Quant & Data Science Bridging

### 1. Portfolios as Linear Combinations
In quantitative portfolio management, let $\mathbf{r}_1, \mathbf{r}_2, \dots, \mathbf{r}_n \in \mathbb{R}^T$ represent vectors of historical returns for $n$ distinct assets over $T$ time periods. Let $w_1, w_2, \dots, w_n$ represent the portfolio allocation weights. 

The portfolio's historical return vector $\mathbf{r}_p$ is a **linear combination** of the individual asset vectors:

$$\mathbf{r}_p = \sum_{i=1}^{n} w_i \mathbf{r}_i$$

*   The allocation weights are constrained such that $\sum_{i=1}^{n} w_i = 1$. The set of all possible long-only portfolios ($w_i \ge 0$) forms a convex hull within the span of the asset vectors.

### 2. Norm Penalties in Regularization: Ridge vs. Lasso Regression
In statistical modeling and machine learning, we solve regression models by minimizing a loss function. To prevent overfitting, we add a norm-based penalty to the regression coefficients vector $\boldsymbol{\beta} \in \mathbb{R}^n$:

*   **Ridge Regression ($L_2$ Regularization)**:

$$\min_{\boldsymbol{\beta}} \|y - X\boldsymbol{\beta}\|_2^2 + \alpha \|\boldsymbol{\beta}\|_2^2$$

    *   *Mathematical Behavior:* The quadratic $L_2$ penalty $\|\boldsymbol{\beta}\|_2^2 = \sum \beta_i^2$ keeps the weights small but non-zero. It shrinks coefficients smoothly, which is ideal for managing high correlation (multicollinearity) in features.
*   **Lasso Regression ($L_1$ Regularization)**:

$$\min_{\boldsymbol{\beta}} \|y - X\boldsymbol{\beta}\|_2^2 + \alpha \|\boldsymbol{\beta}\|_1$$

    *   *Mathematical Behavior:* The absolute $L_1$ penalty $\|\boldsymbol{\beta}\|_1 = \sum |\beta_i|$ features sharp corners (is non-differentiable at $\beta_i = 0$). This forces negligible features to have coefficients of **exactly zero**, performing automated **feature selection** and yielding sparse models.