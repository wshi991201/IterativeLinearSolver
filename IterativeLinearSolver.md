---
Name: Wangzhuo Shi
Topic: [40.md]
Title: What are Iterative linear solvers and where do they apply?
----
# IterativeLinearSolver
## Introduction
Iterative linear solvers are numerical methods used to approximate the solution to a system of linear equations. As opposed to direct linear solvers which attempt to obtain the exact solution via a finite number of operatiors, these solvers iteratively refine an initial guess until a sufficiently accurate solution is obtained. 
## Why
Why would we use iterative linear solver to approximate the solution of $Ax = b$ when we have learned to use LU, QR or Cholesky factorization to directly solve it?
  1. All the methods mentioned above take $O(n^3)$ time! If our iterative solver can arrive at a close-enough approximation with a short amount of iterations, we would be able to save so much time! 
  2. Somtimes if $A$ is sparse, e.g. $A$ only contains $O(n)$ nonzero entries, intermediate operation in direct methods such as matrix inversion may degenerate the complexity of our problem to $O(n^2)$ by introducing additional nonzero entries.
  3. Up until this point, iterative methods are often the only approach to nonlinear equations. Remeber in Newton's methods, we are actually solving the linear system equtaion $f(x_n) - hf'(x_n) = 0$ for $h$ per iteration. Since the overall scheme to solving nonlinear equations is iterative, we may as well use iterative methods per iteration if that means time reduction.
## Stationary iterative methods
To solve $Ax = b$, the iteration scheme for stationary iterative methods can be represented in the following form: 
> ### $x^{(k+1)} = Bx^{(k)} + c$ 

where $B, c$ do not depend on w.r.t iteration count. 
Although the stationary iterative methods are no longer popular in today's application, it is worth discussing the Jacobi Method as it serves as a prototype for other more advanced Iteratvie Methods.
### Jacobi Method
Consider one of the $n$ equations represented by some matrix $A$. We can write the following,
> ### $\sum_{j} a_{ij}x_j = b_i$

Given some $x_{(k)}$. Assume the other entries are correct, how would one caculate entry $x_i$ for $x_{(k+1)}$? 
> ### $x_i^{(k+1)} = \frac{b_i - \sum_{j\neq i} a_{ij}x_j^{(k)}}{a_ii}$

We can also rewrite this in matrix form. Decompose matrix $A = L + D + U$ where $A, L, D, U$ are $N*N$ matrices, $L$ consists of the lower triangular part of $A$ and zeors, $D$ consists of the diagonal of A and zeros, $D$ consists of the upper triangular part of $A$ and zeors. Then the formula for the update rule becomes,

> ### $x^{(x+1)} = D^{-1}(L + U)x^{(k)} + D^{-1}b$

or if we substitue $L+U$ with $A-D$,

> ### $x^{(x+1)} = (I - D^{-1}A)x^{(k)} + D^{-1}b$

There, we can map the above equation to the definition $x^{(k+1)} = Bx^{(k)} + c$ where $B = I - D^{-1}A$, $c = D^{-1}b$!

#### Convergence
The standard convergence condition (for any iterative method) is when the spectral radius of the iteration matrix is less than 1(Beyond the scope of this project). In the case of the Jacobi Method, it translates to $\rho(D^{-1}(L + U)x^{(k)}) < 1$. Thus, a digonally dominant matrix would be a sufficient condition for convergence. For ill-conditioned problems, Jacobi Method converges very slowly or diverge in some cases.
Note that the Jacobi Method does not converge even for some symmetric, positive-definite matrix. 

#### Pseudocode
```
Input: initial guess x0, matrix A, vector b, convergence criterion tol, max iteration K
Output: solution when convergence is reached or none if iteration does not converge/max iteration is met
Comments: pseudocode based on the element-based formula above

k = 0
n = length of a row in A
while (convergence not reached & max iteration not reached) do
    for i = 1 to n:
        sum = 0
        for j = 1 to n:
            if j ≠ i then
                sum = sum + aij * xj
            end
        end
        xi = (bi − sum) / aii
    end
    k = k + 1
end
```
Done! There we finished the Jacobi Method. Of course there are many variants of it. To name a few, Gauss-Seidel Method, Successive Over-Relaxation (SOR)... The Gauss-Seidel Method substitutes the new $x_i$ to the solution vector and immediately put to use DURING each iteration. SOR bascially introduces a wegiht parameter to our update rule ,
> ### $x^{(k+1)} = (I - wD^{-1}A)x^{(k)} + wD^{-1}b$

Now that we have gained some knowledge of stationary iterative linear solvers, it's time to move on to their modern-day counterparts, the Krylov subspace methods.

## Krylov subspace methods
Before we officially delve into any specific methods, let's first understand what Krylov subspaces are and hopefully we can gain some intuition along the way. 
First consider how we discover convergence in iterative methods. One would immediately relate convergence to the magnitude of error at each iteration. However, $e_{(k+1)} = x^* - x^{(k)}$ is impossible to caculate unless we know the exact solution $x^{*}$. Thus, alternatively we use $k_{th}$ residual, $r^{(k)} = b - Ax^{(k)}$ to check for convergence. Then,
> ### $r^{(k)} =  b - Ax^{(k)} = -A(x^{(k)} - x^*) = - Ae^{(k)}$

Recall the generic form of stationary iterative methods we just went through, $x^{(x+1)} = (I - D^{-1}A)x^{(k)} + D^{-1}b$. Now, assume $D  = I$ and let $B = I - A$. Then, we can rewrite $k_{th}$ residual as the following:
> ### $r^{(k)} =  b - Ax^{(k)} = b + Bx^{(k)} - x^{(k)}  = x^{(k+1)} - x^{(k)}$

And then Jocabi Iteration becomes $x^{(k+1)} = x^{(k)} + r^{(k)}$. The new residual $r^{(k+1)} = b - Ax^{(k+1)} = b - A(x^{(k)} + r^{(k)}) =
r^{(k)} -  Ar^{(k)} = Br^{(k)}$. Hey, we just find ourselves a recursive way to compute $r^{(k)}$!
So $k_{th}$ residual $r^{(k)}$ corresponds to $p_{k}(I-A)r^{(0)} \in$ span $\[ r^{(0)}, Ar^{(0)}, A^2r^{(0)},...,A^kr^{(0)} \]$, where $p_k$ is a polynomial of exact degree $k$.
Moreover, we can also represent $x^{(k)}$ in the following way,
> ### $x^{(k)} =  x^{(0)} + r^{(0)} + r^{(1)} + ... + r^{(k-1)} = x^{(0)} + q_{k-1}(B)r^{(0)}$

,where $q_{k-1}$ is a polynomial of exact degree $k-1$.
Therefore, $x_k$ lies in the affine space of $x^{(0)} + $ span $\[ r^{(0)}, Ar^{(0)}, A^2r^{(0)},...,A^{k-1}r^{(0)} \]$, which translates to a shift in the subspace $r^{(k-1)}$. 
Now, we have finally gatherd everything we need to understand what is Krylov subspace. Let's define it!

![alt text][definition] [1]

[definition]: https://github.com/wshi991201/IterativeLinearSolver/blob/main/Screen%20Shot%202023-12-07%20at%204.01.38%20PM.png 



Nice! Now, we are ready to tackle Krylov Space Methods. Krylov subspace methods operate by constructing a basis from a series of matrix powers applied to the initial residual, forming the "Krylov sequence." The solution approximations are subsequently generated by minimizing the residual within this subspace. 

### Conjugate Gradient(CG) Method
#### Why
Consider symmetric, positive-definite $n*n$ matrix $A$(Context for the entire CG section). Recall Gradient Descent is an optimization algorithm used to minimize some objective function by iteratively updating parameters in the direction opposite to the gradient of the objective function. In the case of our $A$, one can easily prove solving $Ax = b$ is equivalent to minimize function $\frac{1}{2} x^TAx - b^Tx + c$ for $c \in \mathbb{R}$(This is quite literally an objective of minimizing the residual, described by Krylov subspace methods). Now that we see how we can use Gradient Descent to solve linear system of equations, why invent other methods? Well, it takes $O(n^2)$ time per iteration($Ax$, Matrix-vector Multiplication) but we don't know hom many iterations we need before convergence. That means if the number of iterations is equal to or greater than O(n), we may as well use Gaussian Elimination, which is $O(n^3)$, to solve for exact solutions! In poor-conditioned cases, Gradient Descent sometimes takes a huge number of iterations to converge. Thus, in order for our method to converge definitely in $O(n^3)$ time, we introduce the Conjugate Gradient(CG) Method.
#### How
##### Conjugate Vector Driect Method 
First, let's review what conjugate vectors are. Two nonzero vectors $u,v$ are said to be conjugate w.r.t matrix $A$ if $u^TAv = 0$. This is also a symmetric relation. If there is a set $P$ of $n$ mutually conjugate vectors with respect to $A$,  then $P$ forms a basis for $\mathbb{R}^n$. Accordingly, we can reprsent the exact solution as a linear combination of these basis vectors in $P$, $x^{*} = \sum_j a_jp_j$. And,
> ### $Ax^{*} = \sum_j a_jAp_j$.

Where $a_j$ is some scalar. Now try multiply $Ax=b$ by some conjugate vector $P_i$, we get $p_i^TAx^{*} = p_i^Tb = \sum_j a_jp_i^TAp_j = 
a_ip_j^TAp_i$. Thus, after we find the $n$ conjugate vectors, we can directly solve for these $a_j$ and compute the exact solution.
##### Iteration Scheme
We ask ourselves the question: is it really necessary to find all $n$ A-conjugate vectors to directly solve the method even when $n$ is very large? It would be nice if we can pick the next conjugate vector as we approximate the solution. Here comes the iterative approach.
Still gradient descent, still minimizing the residual but we are using conjugate vectors as gradient and the scalar $a$ in front of the conjugate vector as line search parameter. We first choose the negative gradient of our objective function $\frac{1}{2} x^TAx - b^Tx + c$, or $Ax-b$ evaluated at $x^{(0)}$ to be the first conjugate vector $p_0$. This happens to be the residual at the first iteration as well. 
Remeber we need to find conjugate vectors along the way. In future iterations, we apply Gram-Schmidt orthonormalization to compute $p_k$, that is 
> ### $p_k = r^{(k)} - \sum_{j} \frac{p_j^TAr^{(k)}}{p_j^TAp_j}p_j, \forall j < k$

And recall our equation for scalar parameter caculation from the direct method above. Since we don't have all the conjugate vectors, we aim to have $a_k$ account for the current residual $r^{(k)}$. We have,
> ### $a_k = \frac{p_k^TAr^{(k)}}{p_k^TAp_k}p_k$

Fianlly, we update $x^{(k)} = x^{(k+1)} + a_kp_k$. We are done!

##### Convergence
By design, CG is guaranteed to converge within $n$ interations and thus performs better than normal Gradient Descent. "[Given cond $A=k$], the number of iterations needed for conjugate gradient to reach a given error value usually can be bounded by a function of $\sqrt{k}$, whereas bounds for convergence of gradient descent are proportional to $k$." [2]


#### Code
The following is a numpy implementation

```python
import numpy as np

def conjugate_gradient(A, b, x0, tol=1e-6, max_iterations=None):
    """
    Conjugate Gradient algorithm for solving Ax = b.

    Parameters:
    - A: Symmetric positive definite matrix
    - b: Right-hand side vector
    - x0: Initial guess for the solution
    - tol: Tolerance for convergence
    - max_iterations: Maximum number of iterations (optional)

    Returns:
    - x: Approximate solution vector
    - iterations: Number of iterations performed
    """

    # Initialization
    x = x0
    r = b - A @ x
    p = r
    iterations = 0

    while np.linalg.norm(r) > tol:
        Ap = A @ p
        alpha = np.dot(r, r) / np.dot(p, Ap)
        x = x + alpha * p
        r_new = r - alpha * Ap
        beta = np.dot(r_new, r_new) / np.dot(r, r)
        p = r_new + beta * p
        r = r_new

        iterations += 1

        # Check for maximum iterations
        if max_iterations is not None and iterations >= max_iterations:
            break

    return x, iterations

```

### Generalized Minimal Residual(GMRES) Method
#### Why
At this point, you are probably wondering how we are gonna solve $Ax=b$ iteratively if $A$ does not have all the pretty properties like before. There you go. Generalized Minimal Residual(GMRES) Method will give you answer! Note that from this point on, we use subscripts to represent iteration number. 
#### How
Take a look at [Arnoldi Iteration](https://en.wikipedia.org/wiki/Arnoldi_iteration). For our purposes, Arnoldi iteration find $n$ orthonormal vectors to form a basis for $n_th$ Krylov subspace. From the definition of the Krylov subspace, we know that $x_n \in x_0 + K_n \implies x_n = x_0 + Q_ny_n$, where $Q_n$ is the orthonormal basis given by the Arnoldi process and $y_n$ is the vector of parameters we use to minimize the residue. In addition, the Arnoldi process gives $n+1 * n$ upper Heisenberg matrix $H_n$ that satisfies $AQ_n = Q_{n+1}H_n$. The combination of these matrices essentially reduce the work we do per iteration to a least-square problem that finds best $y_n$.

##### Iteration Scheme
Apply Arnoldi Iteration to find $Q_n$. Then, find $||yn_n||$ that minimizes the residual $r_n$. Compute $x_{n+1} = x_0 + Q_ny_n$. 

##### Convergence
The method converges regardlessly within n iterations. GMERS does not care about symmetry, positive-definiteness, etc. For spares matrix, the time it takes to converge may even reduce to linear. However, if $A$ is nearly singular or $b$ is nearly orthogonal to the Krylov subspace, the convergency rate will be slow.

##### Code
``` python
import numpy as np
from scipy.linalg import solve, norm, qr

def gmres(A, b, x0=None, tol=1e-6, max_iterations=None, restart=None):
    """
    Generalized Minimal Residual (GMRES) algorithm for solving Ax = b.

    Parameters:
    - A: Coefficient matrix
    - b: Right-hand side vector
    - x0: Initial guess for the solution
    - tol: Tolerance for convergence
    - max_iterations: Maximum number of iterations (optional)
    - restart: Number of iterations before restarting (optional)

    Returns:
    - x: Approximate solution vector
    - iterations: Number of iterations performed
    """

    n = len(b)
    m = restart if restart is not None else n

    # Initialization
    x0 = np.zeros(n) if x0 is None else x0
    r0 = b - A @ x0
    q = [r0 / norm(r0)]

    H = np.zeros((m + 1, m))
    beta = norm(r0)
    H[0, 0] = beta

    x = x0
    for k in range(m):
        # Arnoldi process
        v = A @ q[k]
        for j in range(k + 1):
            H[j, k] = np.dot(q[j], v)
            v = v - H[j, k] * q[j]

        H[k + 1, k] = norm(v)
        q.append(v / H[k + 1, k])

        # Solve the least squares problem
        y, _, _, _ = np.linalg.lstsq(H[: k + 1, : k + 1], beta * np.eye(k + 1))

        # Update solution
        x = x0 + q[0] * y[0]
        for j in range(1, k + 1):
            x = x + q[j] * y[j]

        # Check for convergence
        if norm(H[k + 1, k]) < tol:
            break

    return x, k + 1
```
### Other Methods
Krylov subspace methods is a big, fast-growing family of methods. For example, derived from methods we just saw, GC, GMRES, minimal residual method(MINRES) can be used to solve linear systems with the symmetric but positive-indefinite coefficient matrix; biconjugate gradient method(BiCG) is applied to solve with the non-symmetric matrix. In addition, technique such as preconditioning are developed for iterative numerical methods to improve the convergence and efficiency of the iterative process by transforming the original linear system into an equivalent, better-conditioned one. And, there are many, many more for us to discover and even to invent ourselves. 

## Discussion
Now that we have a basic understanding of iterative linear solvers and what each method is good for, it's time that we talk about the problems one may encounter if he/she does not first analyze the problem. First is obviously result. Using the wrong method on the given question will lead to slow convergence or even divergence. Imagine you are training a large-scale dataset for machine learning. It could mean waste of time in the unit of months or years. Meanwhile, the computation cost such as memory usage, GPU/CPU usage and physical cost such as eletricity bill, inevitably become lareger. In conclusion, we need to conduct a thorough analysis before deciding which iteration method to use. Otherwise, it is going to be a waste of time, space and money.

# Reference
1. Gutknecht, M.H. (2007). A Brief Introduction to Krylov Space Methods for Solving Linear Systems.
2. Iterative Linear Solvers - Stanford University. (n.d.). https://graphics.stanford.edu/courses/cs205a-13-fall/assets/notes/chapter10.pdf 
3. Hestenes, Magnus R.; Stiefel, Eduard (December 1952). "Methods of Conjugate Gradients for Solving Linear Systems" (PDF). Journal of Research of the National Bureau of Standards. 49 (6): 409. doi:10.6028/jres.049.044.
