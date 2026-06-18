# ukfLib: Unscented Kalman filter C library

The idea of the library is to deliver free open source C implementation on UKF with different examples, documentation. Currently there are two available examples. My wish is to extend and test UKF library with more examples from different areas in order to be useful for more people. Any ideas, feedback and help are welcome.

What is included:
- ukfLib.c, ukfLib.h : Unscented Kalman filter implementation.
- ukfCfg.c,ukfCfg.h : These files contain UKF matrix and mathematical state space models/equation(state and measurement equations).They are specific for every different nonlinear filtering problem. In current project these files adress the issue from Matthew B Rhudy tutorial(Understanding Nonlinear Kalman Filters, Part II:An Implementation Guide):
https://yugu.faculty.wvu.edu/files/d/2cbb566f-9936-4033-bb1c-6d887c30d45a/irl_wvu_online_ukf_implementation_v1-0_06_28_2013.pdf
- ukfCfg1.c,ukfCfg1.h : This configuration is used for estimation of angle and angular speed of free pendulum.
- ukfCfg2.c,ukfCfg2.h : nonlinear state estimation for the Van der Pol oscillator.
- ukfTest.c : file is used for both simulation and unit testing. It also show how filter should be initialized and used in practice.

Why this project is useful:
- Implementation could be easy adapted for every different nonlinear filtering problem only by creation of own
  ukfCfgN.c,ukfCfgN.h files.
- Implementation follow the prety good Matthew Rhudy tutorial step by step and solve the same example. Users could easy compare the difference between matlab and C implementation.

For more information please see UKF wiki page: https://github.com/ivo-georgiev/ukfLib/wiki

## Contact

For technical questions, architectural discussions regarding `ukfLib`, please reach out to the author at:
* **Email:** ivo.pl.georgiev@gmail.com

[![View ivo-georgiev/ukfLib on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://nl.mathworks.com/matlabcentral/fileexchange/65092-ivo-georgiev-ukflib)
