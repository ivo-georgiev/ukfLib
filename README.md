# ukfLib: Unscented Kalman filter C library
Посвещавам тази работа на моите любими деца! Приятно пътуване във времето на динозаврите, когато бяхте съвсем мънички! Нека искрата на живота разпали във вас огън от любознателност, старателност, решителност и доброта. Учете се от всеки който може и иска да сподели с вас знание. Обичайте се и се пазете от алчност и завист. Вярвайте и се доверявайте най-вече на себе си!
Обичам ви толкова много, Тати!

The idea of the library is to deliver free open source C implementation on UKF with different examples, documentation and GUI for easy filter configuration. Currently there are two available examples. They show how to use code and how to import and test it as s-function in Simulink(optional). My wish is to extend and test UKF library with more examples from different areas in order to be useful for more people. Any ideas, feedback and help are welcome.

What is included:
- mtxLib.c, mtxLib.h : Library with some basic matrix functions required for UKF implementation(cholesky decomposition, matrix inversion, matrix operation, i.e).
- ukfLib.c, ukfLib.h : Unscented Kalman filter implementation.
- ukfCfg.c,ukfCfg.h : These files contain UKF matrix and mathematical state space models/equation(state and measurement equations).They are specific for every different nonlinear filtering problem. In current project these files adress the issue from Matthew B Rhudy tutorial(Understanding Nonlinear Kalman Filters, Part II:An Implementation Guide):
https://web.statler.wvu.edu/~irl/IRL_WVU_Online_UKF_Implementation_V1.0_06_28_2013.pdf
- ukfCfg1.c,ukfCfg1.h : This configuration is used for estimation of angle and angular speed of free pendulum.
- ukfCfg2.c,ukfCfg2.h : nonlinear state estimation for the Van der Pol oscillator.
- test2.c : file is used for both simulation and unit testing. It also show how filter should be initialized and used in practice.
- Matlab model(R2013b) with s-fun for simulation and verification of each configuration(cfg1,cfg2,..,cfgN) in Simulink enviroment.This model use files listed above to build and simulate s-function.
- Matlab GUI tool for easy filter configuration and automatic file generation(ukfCfgN.c ukfCfgN.h) 

Why this project is useful:
- There is not any specific external dependency. All needed sources are available in repository.
- Implementation could be easy adapted for every different nonlinear filtering problem only by creation of own
  ukfCfgN.c,ukfCfgN.h files.
- Implementation follow the prety good Matthew Rhudy tutorial step by step and solve the same example. Users could easy compare the difference between matlab and C implementation.
- The project contain simulink model(R2013b) which show how s-function could be created(s-function builder) and simulated from given set of source files. There is two steps needed to simulate UKF filter in Simulink. First step is to load data in WS by running of ukf.m. The second step is to build s-fun and run simulink model. Filter behavior could be checked at every simulation step and compared with expected data generated from m-script implementation(applied also in tutorial). 

For more information please see UKF wiki page: https://github.com/ivo-georgiev/ukfLib/wiki

[![View ivo-georgiev/ukfLib on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://nl.mathworks.com/matlabcentral/fileexchange/65092-ivo-georgiev-ukflib)
