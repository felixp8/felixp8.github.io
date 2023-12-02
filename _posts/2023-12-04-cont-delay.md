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
\begin{bmatrix}
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
1) we'd need to also store all external inputs $$ u(t) $$ for $$ [t - \tau, t] $$
2) we'd need to exactly initialize this copied system to exactly $$ x(t - \tau) $$
Instead, if we could estimate $$ \dot{x}_{i}(t - \tau) $$ purely from observables at time $$ t $$ and allow for some margin of error in initializing $$ y $$ where $$ y $$ still converges to $$ x_{i}(t - \tau) $$, that would be potentially less accurate but more useful. 

## Prior work: Finite differences

Though I didn't do much research on prior work, I did look at the first search result, which was this paper[^1]. The method proposed there is very simple: approximate the time derivative with the difference between $$ x_{i}(t) $$ and $$ y(t) $$ divided by $$ \tau $$. If $$ y(t) $$ is equal to $$ x_{i}(t - \tau) $$, then the estimate is equal to $$ (x_{i}(t) - x_{i}(t - \tau)) / \tau $$, which converges to $$ \dot{x}_{i}(t - \tau) $$ as $$ \tau \to 0 $$. In general, we can expect that the estimate is pretty decent as long as $$ \tau $$ is small relative to the timescale of $$ x_{i} $$.

As you can imagine, this is not very accurate for larger $$ \tau $$, so instead of going straight for $$ \tau $$, you can instead break up the delay into many smaller intervals where each derivative estimate is likely to be more accurate. The resulting system would look something like this:

$$ \dot{x} = f(x,u) $$
$$ \dot{y}_1 = (x_i - y_1) / \Delta t $$
$$ \dot{y}_2 = (y_1 - y_2) / \Delta t $$
$$ ... $$
$$ \dot{y}_m = (y_{m-1} - y_m) / \Delta t $$

where $$ m = \tau / \Delta t $$. If it helps to see it in matrix form, the above could also be written like this:

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

where $$ f_{lin} $$ is whatever component of $$ f $$ can be represented as a linear function of $$ x $$ and $$ f_{nonlin} $$ is everything else. (The weird notation of $$ f_{lin} $$ and $$ f_{nonlin} $$ is irrelevant though - the main idea is the rest of the matrix!)

Notice that if our integration timestep is also $$ \Delta t $$, then we have recreated the original discrete system. (The -1's on the diagonal would be canceled out by adding the identity matrix, since discrete systems output the next state, not the difference between the current and next state.)

Notice also that the system sort of self-corrects for inaccurate initializations. If $$ \Delta t $$ is small enough, we can suppose that all gradient estimates are accurate when $$ y $$ values are accurate, e.g. $$ \dot{y}_1(t) = (x_i(t) - y_1(t)) / \Delta t \approx \dot{x}_{i}(t - \Delta t) $$ when $$ y_1(t) \approx x_{i}(t - \Delta t) $$. If initialization is bad, so $$ y_1(t) < x_{i}(t - \Delta t) $$, then $$ \dot{y}_1(t) = (x_i(t) - y_1(t)) / \Delta t > \dot{x}_{i}(t - \Delta t) $$, which compensates for the inaccurate initialization. The same is true if $$ y_1(t) > x_{i}(t - \Delta t) $$.

This method can also be extended with central differences or unevenly-spaced $$ y $$ points, which may give better accuracy[^1][^2].

[^1]: Sun, Jian-Qiao. A method of continuous time approximation of delayed dynamical systems. In *Communications in Nonlinear Science and Numerical Simulation*, 2009.

[^2]: Butcher, Eric A. and Bobrenkov, Oleg A. On the Chebyshev spectral continuous time approximation for constant and periodic delay differential equations. In *Communications in Nonlinear Science and Numerical Simulation*, 2011.
