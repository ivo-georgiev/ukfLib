clear all
%% Before filter execution
% System properties
T = 0.1; % Sampling time
N = 600; % Number of time steps for filter
N1 = 20; % Station 1 North coordinate
E1 = 0; % Station 1 East coordinate
N2 = 0; % Station 2 North coordinate
E2 = 20; % Station 2 East coordinate
% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;
% Step 2: Define noise assumptions
Q = diag([0 0 4 4]);
R = diag([1 1]);
% Step 3: Initialize state and covariance
x = zeros(4, N); % Initialize size of state estimate for all k
x(:,1) = [0; 0; 50; 50]; % Set initial state estimate
P0 = eye(4,4); % Set initial error covariance
% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
w = sqrt(Q)*randn(4, N); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(2, N); % Generate random measurement noise (from assumed R)
xt = zeros(4, N); % Initialize size of true state for all k
xt(:,1) = [0; 0; 50; 50] + sqrt(P0)*randn(4,1); % Set true initial state
yt = zeros(2, N); % Initialize size of output vector for all k
for k = 2:N
xt(:,k) = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1]*xt(:,k-1) + w(:,k-1);
yt(:,k) = [sqrt((xt(1,k)-N1)^2 + (xt(2,k)-E1)^2); ...
sqrt((xt(1,k)-N2)^2 + (xt(2,k)-E2)^2)] + v(:,k);
end
time=0:0.1:length(yt)*0.1-0.2;
exp =[time' xt(:,2:end)'] ;
meas=[time' yt(:,2:end)'];
%load('simdata.mat')

%% Initialize and run EKF for comparison
xe = zeros(4,N);
xe(:,1) = x(:,1);
P = P0;
F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1]; % Linear prediction
for k = 2:N
% Prediction
x_m = F*xe(:,k-1);
P_m = F*P*F' + Q;
% Observation
y_m = [sqrt((x_m(1)-N1).^2 + (x_m(2)-E1).^2); ...
sqrt((x_m(1)-N2).^2 + (x_m(2)-E2).^2)];

H = [(x_m(1)-N1)/sqrt((x_m(1)-N1)^2 + (x_m(2)-E1)^2), ...
(x_m(2)-E1)/sqrt((x_m(1)-N1)^2 + (x_m(2)-E1)^2), 0, 0; ...
(x_m(1)-N2)/sqrt((x_m(1)-N2)^2 + (x_m(2)-E2)^2), ...
(x_m(2)-E2)/sqrt((x_m(1)-N2)^2 + (x_m(2)-E2)^2), 0, 0];
% Measurement Update
K = P_m*H'/(H*P_m*H' + R); % Calculate Kalman gain
xe(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = (eye(4)-K*H)*P_m; % Update covariance estimate
end
clc
%% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0
for k = 2:N
% Step 1: Generate the sigma-points
sP = chol(P,'lower'); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1), x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];
% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
% chi_m = "chi minus" = chi(k|k-1)
chi_m = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1]*chi_p;
x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end
% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [sqrt((chi_m(1,:)-N1).^2 + (chi_m(2,:)-E1).^2); ...
sqrt((chi_m(1,:)-N2).^2 + (chi_m(2,:)-E2).^2)];
y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end
% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
x(:,k);
P = P_m - K*Pyy*K'; % Update covariance estimate
end
%% Display results
% figure(1);
% t = T*(1:N);
% for i = 1:4
% subplot(4,2,2*i-1); plot(t,x(i,:),'b-', t,xe(i,:),'g-.', t,xt(i,:),'r--', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel(['x_',num2str(i)]); grid on; legend('UKF','EKF','True');
% subplot(4,2,2*i); plot(t,x(i,:)-xt(i,:),'b-', t,xe(i,:)-xt(i,:),'g-.', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel(['\Deltax_',num2str(i)]); grid on; legend('UKF','EKF');
% end
exp_ukf =[time' x(:,2:end)'] ; 

