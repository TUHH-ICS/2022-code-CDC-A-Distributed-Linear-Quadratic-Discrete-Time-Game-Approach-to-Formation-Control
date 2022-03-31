# 2022-code-CDC-A-Distributed-Linear-Quadratic-Discrete-Time-Game-Approach-to-Multi-Agent-Consensus

## General
This repository contains an implementation of the algorithms and simulations described in 
> Prima Aditya and Herbert Werner, "A Distributed Linear Quadratic Discrete-Time Game Approach to Multi-Agent Consensus" submitted for CDC 2022

## Simulation
The main code `decoupling_lqdtg.m` contains everything in one file. This distributed framework depends on the predefined communication graph structure. In this case, we took four vertices and six edges as an example. Settings that can be adjusted are the finite time (tf) and the number of prediction horizons (Np)

## Supplementary files
Algorithms of Nash strategy (Problem 1) and distributed framework in receding horizon fashion (Problem 2) that can't be presented in the paper are stored in the supplementary files folder

## Prerequisites
The simulation code in this repository was tested in the following environment:
- *Windows 11* Pro Version 21H2
- *Matlab* 2020a
