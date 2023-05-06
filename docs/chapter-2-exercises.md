# Exercises: Basics

## Matrix Algebra

$$\begin{align}
3x + 2y + 5z &= 0\\
7x + 6y + 4z &= -2\\
 x + 3y + 2z &= -6\\
\end{align}$$

``` {.julia file=src/exercises/ch2.jl #exercise-ch2}
function exercise_1()
    A = [ 3 2 5 ;
          7 6 4
          1 3 2 ]
    b = [ 0, -2, -6 ]
    println("Exercise 1: A \\ b = $(A \ b)")
end

exercise_1()
```

The answer can be checked using Gaussian elimination.

## The First Numerical Model
Equation 2.26. We have a 1D diffusion model of a Newtonian fluid adjacent to a wall. The wall (at $y = 0$) moves with velocity $V_0$, and the fluid is dragged along with a velocity $V(y)$

$$\partial_t V - \nu \partial_y^2 V = 0$$

Solve using FTCS (Forward-in-time, Centered-in-space).

$$V_j^{n+1} = sV_{j-1}^n + (1 - 2s)V_j^n + sV^n_{j+1},$$

where $s = \nu \Delta t / \Delta y^2$.

We are given an analytic solution

$$V = V_0 \left\{\sum_{n=0}^{\infty} {\rm erfc}(2n\eta_1 + \eta) - \sum_{n=0}^{\infty} {\rm erfc}(2(n+1)\eta_1 - \eta)\right\},$$

where $\eta_1 = h / 2\sqrt{\nu} t$ and $\eta = y/2\sqrt{\nu}t$.

``` {.julia #exercise-ch2}

```
