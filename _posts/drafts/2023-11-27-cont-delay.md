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

so that $$ x_3(t + 1) = x_1(t) $$. Of course, we can only ever create delays of exactly one timestep, and more generally in discrete-time systems state variables can only affect each other with a one timestep delay.

In continuous-time systems, state variables affect each other instantaneously. It is therefore not obvious how to recreate this sort of time-delayed relationship between variables, though this occurs often in real-world systems, e.g. in interacting brain regions. This problem has been extensively studied in a variety of fields, but because I am too lazy to look into it, I just messed around a bit myself and found a solution that is probably pretty useless.

## Problem statement

Let me try to formulate the problem a little more clearly. We have a continuous-time dynamical system

$$ \dot{x} = f(x,u) $$

where $$ x \in \mathbb{R}^{N} $$. We want to introduce a new state variable $$ y $$ so that $$ y(t) = x_{i}(t - \tau) $$ for some index $$ i \in [1, N] $$ and some delay $$ \tau > 0 $$. This means that we want

$$ y(t) = \dot{x}_{i}(t - \tau). $$

In theory, we could maintain an entire delayed copy of $$ x $$ and compute $$ \dot{x}_{i}(t - \tau) $$ just from $f$, but this is generally impractical because
1) we'd need to also store all external inputs $$ u(t) $$ for $$ [t - \tau, t] $$
2) we'd need to exactly initialize this copied system to exactly $x(t - \tau)$
Instead, if we could estimate $$ \dot{x}_{i}(t - \tau) $$ purely from observables at time $$ t $$ and allow for some margin of error in initializing $$ y $$, that would be less accurate but more manageable.

## Prior work: Finite differences

Though I didn't do much research on prior work, I did look at the first search result, which was [1]. The method proposed there is very simple: approximate the time derivative with the difference between $$ y $$ and $$ x_{i} $$ divided by $$ \tau $$. As you can imagine, this is not very accurate for larger $$ \tau $$, so instead of going straight for $$ \tau $$, you can instead break up the delay into many smaller intervals where each derivative estimate is likely to be more accurate. The resulting system would look something like this:

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

where $$ f_{lin} $$ is whatever component of $$ f $$ can be represented as a linear function of $$ x $$ and $$ f_{nonlin} $$ is everything else.

There are a couple interesting things to notice about this. 