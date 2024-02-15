Matlab code for reproducing results presented in the manuscript "**Shaping Robust Dynamic Inversion Control of Neural Cell Dynamics**" (**doi.org/10.21203/rs.3.rs-2986266/v1**). Please contact rongting.yue@uconn.edu for any questions.


File "IDNI_gains_11_27_no_EKF_burst_traj.m" uses the control method in the manuscript to activate a spiking train.

![spike_train](https://github.com/yuerongting/Neuron-Control/assets/42655883/b7719528-5694-4380-9b96-c63cc6998713)


If the initial condition is not accurate and an Extended Kalman Filter is needed, the performance can be checked in "IDNI_gains_EKF.m"

![EKF](https://github.com/yuerongting/Neuron-Control/assets/42655883/bf62c683-8d64-4e46-ad86-de306e4c306d)

Explore different combinations of controller gains "IDNI_gains_11_27_no_EKF_diff_gains.m"

![diff_gains_combination](https://github.com/yuerongting/Neuron-Control/assets/42655883/8a4e69cf-80a9-46f8-9d29-bc72f013c9da)

Comparison of robustness performanc run "IDNI_robust_12_30.m"

![Comparison_all](https://github.com/yuerongting/Neuron-Control/assets/42655883/fa19fdd4-bdd8-43ae-9a1e-501c540b896f)
