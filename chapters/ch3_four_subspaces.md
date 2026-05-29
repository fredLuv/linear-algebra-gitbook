# Chapter 3: The Four Fundamental Subspaces & Dimension Theorem

## 🎯 Core Thesis
The master key to understanding any matrix $A$ (size $m \times n$, rank $r$) lies in the **Four Fundamental Subspaces** defined by its linear transformations. Gilbert Strang formalized this as the **Fundamental Theorem of Linear Algebra**. 

These four subspaces describe the mapping behavior of $A$ from $\mathbb{R}^n$ to $\mathbb{R}^m$. The row space $C(A^T)$ and nullspace $N(A)$ partition the input space $\mathbb{R}^n$, while the column space $C(A)$ and left nullspace $N(A^T)$ partition the output space $\mathbb{R}^m$. 

The dimensions of these spaces are governed by the **Rank-Nullity Theorem**, ensuring a perfect mathematical symmetry between a matrix's information content (rank) and its redundancies (nullspace).

---

## 🔑 Key Mathematical Formulations

Let $A$ be an $m \times n$ matrix with rank $r$ (where $r \le \min(m, n)$).

```
          INPUT SPACE R^n                           OUTPUT SPACE R^m
   ┌───────────────────────────┐             ┌───────────────────────────┐
   │  Row Space C(A^T)         │             │  Column Space C(A)        │
   │  Dimension = r            │ ───[Ax]───► │  Dimension = r            │
   ├───────────────────────────┤             ├───────────────────────────┤
   │  Nullspace N(A)           │             │  Left Nullspace N(A^T)    │
   │  Dimension = n - r        │ ───[Ax=0]─► │  Dimension = m - r        │
   └───────────────────────────┘             └───────────────────────────┘
```

### 1. The Column Space $C(A)$
*   **Definition**: The set of all linear combinations of the columns of $A$. It resides in $\mathbb{R}^m$.
*   **Dimension**: $\dim C(A) = r$
*   **Solvability Axiom**: The linear system $Ax = b$ has a solution if and only if the vector $b$ lies in the column space of $A$ ($b \in C(A)$).

### 2. The Nullspace $N(A)$
*   **Definition**: The set of all vectors $x \in \mathbb{R}^n$ that solve the homogeneous equation $Ax = \mathbf{0}$. It resides in $\mathbb{R}^n$.
*   **Dimension**: $\dim N(A) = n - r$

### 3. The Row Space $C(A^T)$
*   **Definition**: The column space of $A^T$; the set of all linear combinations of the rows of $A$. It resides in $\mathbb{R}^n$.
*   **Dimension**: $\dim C(A^T) = r$
*   *Theorem:* The column space and the row space always have the exact same dimension, equal to the rank $r$, regardless of whether the matrix is tall, wide, or square.

### 4. The Left Nullspace $N(A^T)$
*   **Definition**: The nullspace of $A^T$; the set of all vectors $y \in \mathbb{R}^m$ such that $A^T y = \mathbf{0}$ (or $y^T A = \mathbf{0}^T$). It resides in $\mathbb{R}^m$.
*   **Dimension**: $\dim N(A^T) = m - r$

### 5. The Dimension Theorem (Rank-Nullity)
The dimensions of the input subspaces must sum to the total number of inputs:

$$\text{Rank} + \text{Nullity} = n \implies r + (n - r) = n$$

---

## ✍️ Step-by-Step Worked Example: Finding Basis and Dimension

### Problem:
Given the following $3 \times 4$ matrix $A$ of rank $r = 2$, find the **basis and dimension** for all four fundamental subspaces:

$$A = \begin{bmatrix} 1 & 2 & 0 & 1 \\ 2 & 4 & 1 & 4 \\ 1 & 2 & 1 & 3 \end{bmatrix}$$

### Step 1: Reduce to Row Echelon Form ($R$)
We apply Gaussian row elimination to find the pivot columns:

1.  Subtract $2$ times Row 1 from Row 2 ($R_2 \leftarrow R_2 - 2R_1$):
    $$\begin{bmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 2 \\ 1 & 2 & 1 & 3 \end{bmatrix}$$
2.  Subtract Row 1 from Row 3 ($R_3 \leftarrow R_3 - R_1$):
    $$\begin{bmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 2 \\ 0 & 0 & 1 & 2 \end{bmatrix}$$
3.  Subtract Row 2 from Row 3 ($R_3 \leftarrow R_3 - R_2$):
    $$R = \begin{bmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 \end{bmatrix}$$

*   **Pivots**: The pivots are in Column 1 and Column 3.
*   **Rank**: The number of non-zero rows is $r = 2$.

---

### Step 2: The Row Space $C(A^T)$
The basis for the row space is formed by the non-zero rows of the reduced matrix $R$:

$$\text{Basis of } C(A^T) = \left\{ \begin{bmatrix} 1 \\ 2 \\ 0 \\ 1 \end{bmatrix}, \begin{bmatrix} 0 \\ 0 \\ 1 \\ 2 \end{bmatrix} \right\}$$

*   **Dimension**: $\dim C(A^T) = r = 2$.

---

### Step 3: The Column Space $C(A)$
*Warning:* The column space is preserved by column operations, but **not by row operations**. We must select the pivot columns from the **original matrix $A$** (Columns 1 and 3):

$$\text{Basis of } C(A) = \left\{ \begin{bmatrix} 1 \\ 2 \\ 1 \end{bmatrix}, \begin{bmatrix} 0 \\ 1 \\ 1 \end{bmatrix} \right\}$$

*   **Dimension**: $\dim C(A) = r = 2$.

---

### Step 4: The Nullspace $N(A)$
To find the nullspace, we solve $Rx = \mathbf{0}$. We express the pivot variables ($x_1, x_3$) in terms of the free variables ($x_2, x_4$):

*   Row 2: $x_3 + 2x_4 = 0 \implies x_3 = -2x_4$
*   Row 1: $x_1 + 2x_2 + x_4 = 0 \implies x_1 = -2x_2 - x_4$

We write the general vector $x \in N(A)$:

$$x = \begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{bmatrix} = \begin{bmatrix} -2x_2 - x_4 \\ x_2 \\ -2x_4 \\ x_4 \end{bmatrix} = x_2 \begin{bmatrix} -2 \\ 1 \\ 0 \\ 0 \end{bmatrix} + x_4 \begin{bmatrix} -1 \\ 0 \\ -2 \\ 1 \end{bmatrix}$$

The basis vectors are the coefficients of the free variables:

$$\text{Basis of } N(A) = \left\{ \begin{bmatrix} -2 \\ 1 \\ 0 \\ 0 \end{bmatrix}, \begin{bmatrix} -1 \\ 0 \\ -2 \\ 1 \end{bmatrix} \right\}$$

*   **Dimension**: $\dim N(A) = n - r = 4 - 2 = 2$.

---

### Step 5: The Left Nullspace $N(A^T)$
We solve $A^T y = \mathbf{0}$, which is equivalent to finding the linear combination of rows of $A$ that yields the zero vector. We can do this by tracking row operations on an augmented matrix $[A \mid I]$ reduced to $[R \mid E]$:

$$\left[ \begin{array}{cccc|ccc} 1 & 2 & 0 & 1 & 1 & 0 & 0 \\ 2 & 4 & 1 & 4 & 0 & 1 & 0 \\ 1 & 2 & 1 & 3 & 0 & 0 & 1 \end{array} \right]$$

1.  Row operations:
    *   $R_2 \leftarrow R_2 - 2R_1$
    *   $R_3 \leftarrow R_3 - R_1$
    $$\left[ \begin{array}{cccc|ccc} 1 & 2 & 0 & 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & 2 & -2 & 1 & 0 \\ 0 & 0 & 1 & 2 & -1 & 0 & 1 \end{array} \right]$$
2.  Subtract Row 2 from Row 3 ($R_3 \leftarrow R_3 - R_2$):
    $$\left[ \begin{array}{cccc|ccc} 1 & 2 & 0 & 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & 2 & -2 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 & -1 & 1 \end{array} \right]$$

The row in the right-hand identity matrix corresponding to the zero row in $R$ represents the vector in the left nullspace:

$$\text{Basis of } N(A^T) = \left\{ \begin{bmatrix} 1 \\ -1 \\ 1 \end{bmatrix} \right\}$$

*   **Verification**: Check $A^T y = \mathbf{0}$:
    $$1\begin{bmatrix} 1 \\ 2 \\ 1 \end{bmatrix} - 1\begin{bmatrix} 2 \\ 4 \\ 2 \end{bmatrix} + 1\begin{bmatrix} 1 & 2 & 1 \end{bmatrix}^T = \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}$$
*   **Dimension**: $\dim N(A^T) = m - r = 3 - 2 = 1$.

---

## 📈 Quant & Data Science Bridging

### 1. Column Space as Achievable Portfolio Returns
Let $A$ be a $T \times n$ matrix representing the historical returns of $n$ distinct trading strategies over $T$ days. The columns of $A$ are vectors $\mathbf{a}_i \in \mathbb{R}^T$. 

*   The **Column Space $C(A)$** represents the space of all **achievable portfolio returns** that can be constructed by allocating capital weights $x$ to these strategies.
*   If a client requests a specific daily return profile vector $b \in \mathbb{R}^T$ (e.g., matching a benchmark index), this target return is **only achievable** if $b$ lies in the column space of $A$ ($b \in C(A)$). If $b$ is outside $C(A)$, we must project $b$ onto $C(A)$ to find the closest achievable portfolio (yielding the Ordinary Least Squares solution).

### 2. Nullspace as Market-Neutral Hedge Portfolios
In multi-factor quantitative risk models, let $A$ be an $m \times n$ exposure matrix, where $A_{ij}$ is the exposure of asset $j$ to systemic risk factor $i$ (e.g., beta to growth, momentum, or value).

*   A portfolio allocation vector $x \in \mathbb{R}^n$ represents a **factor-hedged portfolio** if and only if it lies in the **Nullspace $N(A)$**:
    $$Ax = \mathbf{0}$$
*   Every vector $x$ in the nullspace represents a market-neutral portfolio with **absolute zero net exposure** to all $m$ systemic risk factors. Finding the basis of $N(A)$ allows risk managers to construct arbitrage portfolios that are completely insulated from macroeconomic factor movements.
