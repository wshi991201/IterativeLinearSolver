---
Name: Wangzhuo Shi
Topic: 40
Title: What are Iterative Linear Solvers and where do they apply?
----
# IterativeLinearSolver
## Introduction
Iterative linear solvers are numerical methods used to approximate the solution to a system of linear equations. As opposed to direct linear solvers, which attempt to obtain the exact solution via a finite number of operatiors, these solvers iteratively refine an initial guess until a sufficiently accurate solution is obtained. 
## Why
Why would we use iterative linear solver to approximate the solution of $Ax = b$ when we have learned to use LU, QR or Cholesky factorization to directly solve it?
  1. All the methods mentioned above take $O(n^3)$ time! If our iterative solver can arrive at a close-enough approximation with a short amount of iterations, we would be able to save so much time! 
  2. Somtimes if $A$ is sparse, e.g. $A$ only contains $O(n)$ nonzero entries, intermediate operation in direct methods such as matrix inversion may degenerate the complexity of our problem to $O(n^2)$ by introducing additional nonzero entries.
  3. Up until this point, iterative methods are often the only approach to nonlinear equations. Remeber in Newton's methods, we are actually solving the linear system equtaion $f(x_n) - hf'(x_n) = 0$ for $h$ per iteration. Since the overall scheme to solving nonlinear equations is iterative, we may as well use iterative methods per iteration if that means time reduction.
## Stationary Iterative Method
To solve $Ax = b$, the iteration scheme for Stationary Iterative Methods can be represented in the following form: 
> ### $x^{(k+1)} = Bx^{(k)} + c$ 

where $B, c$ do not depend on w.r.t iteration count. 
Although the Stationary Iterative Methods are no longer popular in today's application, it is worth discussing the Jacobi Method as it serves as a prototype for other more advanced Iteratvie Methods.
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

### Krylov subspace methods
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

![alt text][definition]

[definition]: https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Logo Title Text 2"
