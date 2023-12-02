---
layout: post
title: approximation of discrete delays with continuous-time dynamics
date: 2023-11-27 21:01:00
description: some pointless and basic math noodling
tags: math code
categories: dynamical-systems
---

In discrete-time dynamical systems, it’s easy to create a state variable that is a time-delayed copy of another. For example, if we have a linear system 

$$ 
\begin{bmatrix}
    x_1(t + 1) \\ 
    x_2(t + 1) \\ 
    x_3(t + 1) 
\end{bmatrix} = A \cdot \begin{bmatrix}
    x_1(t) \\ 
    x_2(t) \\ 
    x_3(t) 
\end{bmatrix} 
$$

and we want $$ x_3 $$ to be equal to $$ x_1 $$, but delayed by one timestep, we can design A as

$$ 
A = \begin{bmatrix}
    a_{11} & a_{12} & a_{13} \\
    a_{21} & a_{22} & a_{23} \\
    1 & 0 & 0
\end{bmatrix}
$$

so that $$ x_3(t + 1) = x_1(t) $$. Of course, we can only ever create delays of exactly one timestep, and to achieve greater delays we'd need to chain multiple variables each delayed by a timestep. More generally, in discrete-time systems, state variables can only affect each other with a one timestep delay.

In continuous-time systems, state variables affect each other instantaneously. It is therefore not obvious how to recreate this sort of time-delayed relationship between variables, though this occurs often in real-world systems, e.g. in interacting brain regions. This problem has been extensively studied in a variety of fields, but because I am too lazy to look into it, I just messed around a bit myself and found a solution that is probably pretty useless.

## Problem statement

Let me try to state the problem a little more clearly. We have a continuous-time dynamical system

$$ \dot{x} = f(x,u) $$

where $$ x \in \mathbb{R}^{N} $$. We want to make some sort of coupled dynamical system with a new state variable $$ y $$ so that $$ y(t) = x_{i}(t - \tau) $$ for some index $$ i \in [1, N] $$ and some delay $$ \tau > 0 $$. This means that we want

$$ \dot{y}(t) = \dot{x}_{i}(t - \tau). $$

In theory, we could maintain an entire delayed copy of $$ x $$ and compute $$ \dot{x}_{i}(t - \tau) $$ just from $$ f $$, but this is generally impractical because
1. we'd need to also store all external inputs $$ u(t) $$ for $$ [t - \tau, t] $$
2. we'd need to exactly initialize this copied system to exactly $$ x(t - \tau) $$
Instead, if we could estimate $$ \dot{x}_{i}(t - \tau) $$ purely from observables at time $$ t $$ and allow for some margin of error in initializing $$ y $$, that would be potentially less accurate but more useful. 

## Prior work: Finite differences

Though I didn't do much research on prior work, I did look at the first search result[^fn1]. The method proposed there is very simple: approximate 

$$ \dot{x}_{i}(t - \tau) \approx (x_{i}(t) - y(t)) / \tau. $$

If $$ y(t) = x_{i}(t - \tau) $$ exactly, then the estimate is equal to 

$$ (x_{i}(t) - x_{i}(t - \tau)) / \tau, $$

which converges to $$ \dot{x}_{i}(t - \tau) $$ as $$ \tau \to 0 $$.

As you can imagine, this becomes less accurate for larger $$ \tau $$, so instead of going straight for $$ \tau $$, you can break up the delay into many smaller intervals where each derivative estimate is likely to be more accurate. The resulting system would look something like this:

$$ 
\dot{x} = f(x,u) 
$$

$$ 
\dot{y}_1 = (x_i - y_1) / \Delta t 
$$

$$ 
\dot{y}_2 = (y_1 - y_2) / \Delta t 
$$

$$ 
... 
$$

$$ 
\dot{y}_m = (y_{m-1} - y_m) / \Delta t 
$$

which gives a total delay of $$ m \cdot \Delta t $$. If it helps to see it in matrix form, the above could also be written like this:

$$
\begin{bmatrix} 
    \dot{x} \\
    \dot{y}_1 \\
    \\
    \vdots \\
    \\
    \dot{y}_m
\end{bmatrix} = 
\begin{bmatrix} 
    f_{lin} & 0 &  & \dots &  & 0 \\
    \frac{1}{\Delta t} & -\frac{1}{\Delta t} & 0 &  & \dots & \\
    0 & \frac{1}{\Delta t} & -\frac{1}{\Delta t} &  & \dots & \\
     &  & \ddots & \ddots &  & \\
     &  &  & \ddots & \ddots & \\
    0 &  & \dots &  & \frac{1}{\Delta t} & -\frac{1}{\Delta t}
\end{bmatrix}
\begin{bmatrix}
    x \\
    y_1 \\
    \\
    \vdots \\
    \\
    y_m
\end{bmatrix} + 
\begin{bmatrix}
    1 \\
    0 \\
    \\
    \vdots \\
    \\
    0
\end{bmatrix}
f_{nonlin}(x,u)
$$

where $$ f_{lin} $$ is whatever component of $$ f $$ can be represented as a linear function of $$ x $$ and $$ f_{nonlin} $$ is everything else. (My weird notation of $$ f_{lin} $$ and $$ f_{nonlin} $$ is irrelevant though - the main idea is the rest of the matrix!)

There are a few nice things about this method:

1. The dynamics of the delayed variables do not depend explicitly on the inputs at all. Instead, they account for the effect of the inputs only by looking at $$ x_{i} $$ directly.
2. The system sort of self-corrects for inaccurate initializations. If $$ \Delta t $$ is small enough, we can assume that all gradient estimates are accurate when $$ y $$ values are accurate. Then, for example, if $$ y $$ is initialized too small, the gradient estimate from the finite difference will be bigger than the true gradient estimate, which compensates for the initialization error.

The method can be extended/improved with central differences or unevenly-spaced $$ y $$ points, which may give better accuracy[^fn1][^fn2].

Lastly, notice that if our integration timestep is also $$ \Delta t $$, then $$ \pm 1 / \Delta t = \pm 1 $$, and we have recreated the original discrete system. The -1's on the diagonal would be canceled out by adding the identity matrix, since discrete systems output the next state, not the difference between the current and next state.

## Finding an alternate solution by inverting the continuous-to-discrete transformation

Instead of doing something principled like the paper above, I decided to just try to directly convert a discrete-time linear system with delays into a continuous-time one. Thanks to Stack Exchange, I found that the bilinear transform is easily invertible^[fn3]. For some background, the bilinear transform is essentially a trapezoidal approximation of a continuous-time linear system. 

Taking a discrete-time linear system with a delay, like the one shown in the very first section, we can apply the inverse bilinear transform to get a continuous-time system that should display the same behavior. Though I did not expect it to, it worked:

[insert image here]

Doing this a few times, there was a clear pattern to the resulting continuous-time systems. The dynamics of the delayed variable were always given by

$$ \dot{y}(t) = -\dot{x}_{i}(t) + \frac{2}{\tau}(x_{i}(t) - y(t)). $$

This equation did not immediately make sense to me, but ultimately it's quite straightforward if we rearrange it. 

$$ \dot{y}(t) + \dot{x}_{i}(t) = \frac{2}{\tau}(x_{i}(t) - y(t)) $$

$$ (\dot{y}(t) + \dot{x}_{i}(t)) / 2 = (x_{i}(t) - y(t)) / \tau $$

If we assume the ideal of $$ y(t) = x_{i}(t - \tau) $$, we get

$$ (\dot{x}_{i}(t - \tau) + \dot{x}_{i}(t)) / 2 = (x_{i}(t) - x_{i}(t - \tau)) / \tau. $$

Basically, the assumption of this solution is that the average slope between $$ x_{i}(t - \tau) $$ and $$ x_{i}(t) $$ ( $$ (x_{i}(t) - x_{i}(t - \tau)) / \tau $$ ) is right in the middle of the slope at the left edge ( $$ \dot{x}_{i}(t - \tau) $$ ) and the slope at the right edge ( $$ \dot{x}(t - \tau) $$ ). This would be true of the second derivative is constant.

The good thing about this approximation is that it retains all the benefits of the finite difference method, relying only on current observables and compensating for inaccurate initializations, while taking into account the additional information we know about $$ \dot{x}_{i}(t) $$ that was ignored by the finite difference method.

## Generalizing the alternate solution with Taylor expansions

TODO

## References

[^fn1]: Sun, Jian-Qiao. A method of continuous time approximation of delayed dynamical systems. In *Communications in Nonlinear Science and Numerical Simulation*, 2009. [URL](https://www.sciencedirect.com/science/article/pii/S1007570408000610)


[^fn2]: Butcher, Eric A. and Bobrenkov, Oleg A. On the Chebyshev spectral continuous time approximation for constant and periodic delay differential equations. In *Communications in Nonlinear Science and Numerical Simulation*, 2011. [URL](https://www.sciencedirect.com/science/article/pii/S1007570410003539)


[^fn3]: van der Veen, Kwin. discrete-time to continuous-time state space. On *Mathematics Stack Exchange*, 2020. [URL](https://math.stackexchange.com/q/3820405)