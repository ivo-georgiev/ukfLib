# ukfLib: Unscented Kalman filter C library
What is included:
- mtxLib.c, mtxLib.h : Library with some basic matrix functions required for UKF implementation(choleski decomposition, matrix inversion, matrix operation, i.e).
- ukfLib.c, ukfLib.h : Unscented Kalman filter implementation.
- ukfCfg.c,ukfCfg.h : These files contain UKF matrix and mathematical state space models/equation(state and measurement equations).They are specific for every different nonlinear filtering problem. In current project these files adress the issue from Matthew B Rhudy tutorial(Understanding Nonlinear Kalman Filters, Part II:An Implementation Guide):
https://web.statler.wvu.edu/~irl/IRL_WVU_Online_UKF_Implementation_V1.0_06_28_2013.pdf
- test.c : This file is used for both simulation and unit testing. It also show how filter should be initialized and used in practice.
- Matlab model with s-fun for simulation and verification in Simulink enviroment.This model use files listed above to build s-function.

Why this project is useful:
- Implementation could be easy adapted for every different nonlinear filterin problem only by modification on
  ukfCfg.c,ukfCfg.h files.
- Implementation follow the prety good Matthew Rhudy tutorial step by step and solve the same example. Users could compare difference between matlab and C implementation.
- The project contain simulink model which shown how s-function could be created and simulated from given source files. First step is to load data from m-script ukf.m . The second step is to build s-fun and run simulink model. Filter behavior could be checked at every simulation step and compared with m-script implementation from tutorial. 
