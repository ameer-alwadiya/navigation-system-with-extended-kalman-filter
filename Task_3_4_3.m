% clear all;
% close all;
clc; 
% Estimation of the data in the nominal state, to get angle of attack
% measurements
% Template for Kalman Filter (KF), Extended Kalman Filter (EKF) and Iterated Extended Kalman Filter (IEKF)

% Read data files
load('Data/dataTask2.mat');
z_k = d_k;
u_k = c_k;

% Initialization
Ts=dt;     %time step (already provided by the data)
N=length(t);      %total number of steps

% omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r
stdw = [0.01, 0.01, 0.01, 0.01*pi()/180, 0.01*pi()/180, 0.01*pi()/180, 5, 5, 5, 1.5*pi()/180, 1.5*pi()/180, 1.5*pi()/180];  %standard deviation of w

% nu_x nu_y nu_z nu_u nu_v nu_w nu_phi nu_theta nu_psi nu_V nu_alpha nu_beta
stdv = [5, 5, 10, 0.1, 0.1, 0.1, 0.1*pi()/180, 0.1*pi()/180, 0.1*pi()/180, 0.1, 0.1*pi()/180, 0.1*pi()/180];      %standard deviation of v

% x y z u v w phi theta psi v_wxE v_wyE v_wzE
Ex_0 = [d_k(1,1) d_k(1,2) d_k(1,3) d_k(1,10) 0 0 d_k(1,7) d_k(1,8) d_k(1,9) 0 0 0 0 0 0 0 0 0];   %expected value of x_0
stdx_0 = [1 1 1 1 100 100 0.1*pi()/180 0.1*pi()/180 0.1*pi()/180 100 100 100 200 200 200 200 200 200];  %standard deviation of x_0

xhat_km1_km1 = Ex_0; % x(0|0) = E{x_0}
P_km1_km1 = diag(stdx_0.^2);  % P(0|0) = P(0)
Q=diag(stdw.^2);
R=diag(stdv.^2);

n = length(xhat_km1_km1); % n: state dimension
m = size(u_k, 2);     % m: observation dimension
p = size(z_k, 2);     % m: observation dimension
u_km1 = [zeros(1,m); u_k]; % shifted to have the right indices

% Preallocate storage
stdx_cor  = zeros(N, n);  % \sigma(k-1|k-1), standard deviation of state estimation error (hint: diagonal elements of P(k-1|k-1))
x_cor     = zeros(N, n);  % \hat{x}(k-1|k-1), previous estimation
K       = cell(N, 1);   % K(k) Kalman Gain
innov     = zeros(N, p);  % y(k)-y(k|k-1), innovation, with y(k|k-1)=h(\hat{x}(k|k-1),u(k|k-1),k);

for k=1:N
    % Step 1: Prediction
    [t_nonlin, x_nonlin] = ode45(@(t,x) funcf(x, u_km1(k,:), t), [0 Ts], xhat_km1_km1);
    xhat_k_km1 = x_nonlin(end,:); % x(k|k-1) (prediction)

    % Step 2: Covariance matrix of state prediction error / Minimum
    % prediction MSE
    [Phi_km1, Gamma_km1] = funcLinDisDyn(xhat_km1_km1, u_km1(k,:), Ts); % Phi(k,k-1), Gamma(k,k-1)
    P_k_km1 = Phi_km1 * P_km1_km1 * Phi_km1' + Gamma_km1 * Q * Gamma_km1'; % P(k|k-1) (prediction)

    % Step 3: Kalman Gain
    H_k = funcLinDisObs(xhat_k_km1, u_km1(k,:), []);
    Ve = (H_k * P_k_km1 * H_k' + R); % Pz(k|k-1) (prediction)
    K_k = P_k_km1 * H_k' / Ve; % K(k) (gain)
    
    % Step 4: Measurement Update (Correction)
    z_k_km1 = funch(xhat_k_km1,u_km1(k,:),[]); % z(k|k-1) (prediction of output)
    xhat_k_k = xhat_k_km1 + (z_k(k,:) - z_k_km1)*K_k'; % x(k|k) (correction)

    % Step 5: Correction for Covariance matrix of state Estimate error /
    % Minimum MSE
    I_KH = eye(n) - K_k * H_k;
    P_k_k = I_KH * P_k_km1 * I_KH' + K_k * R * K_k'; % P(k|k) (correction)

    % Save data: State estimate and std dev
    stdx_cor(k,:) = sqrt(diag(P_km1_km1)); % \sigma(k-1|k-1) Standard deviation of state estimation error
    x_cor(k,:) = xhat_km1_km1; % \hat{x}(k-1|k-1), estimated state
    K{k,1} = K_k; % K(k) (gain)
    innov(k,:)= z_k(k,:) - z_k_km1;
    
    % Recursive step
    xhat_km1_km1 = xhat_k_k; 
    P_km1_km1 = P_k_k;
end


%% Plots

titles = {'x(1)=x', 'x(2)=y', 'x(3)=z', 'x(4)=u', 'x(5)=v', 'x(6)=w', 'x(7)=phi', 'x(8)=theta', 'x(9)=psi', 'x(10)=v_wxE', 'x(11)=v_wyE', 'x(12)=v_wzE', 'x(13)=b_A_x', 'x(14)=b_A_y', 'x(15)=b_A_z', 'x(16)=b_p', 'x(17)=b_q', 'x(18)=b_r'};

figure;
for i = 1:18
    if ismember(i, 1:12)
        subplot(3, 6, i)
        plot(t(1:N), x_cor_free_bias(1:N,i), 'b', t(1:N), x_cor(1:N,i), 'r')
        title(['Estimation of the state ' titles{i}])
        legend('Nominal', 'Nominal faults')
        grid on
    else
        subplot(3, 6, i)
        plot(t(1:N), x_cor(1:N,i), 'b')
        title(['Estimation of the state ' titles{i}])
        grid on
    end
end
sgtitle('Estimations of the states with fault measurements');

figure;
for i = 1:12
    subplot(2, 6, i)
    plot(t(1:N), innov(1:N,i), 'b')
    title(['Innovation of the state ' titles{i}])
    grid on
end
sgtitle('Innovations of the states');

Innov_nominal = innov;

figure;
plot(t(1:N), innov(1:N,11), 'b')
title('Innovation of Alpha - Angle of Attack')
grid on

%% dynamics in continous time
function x_dot = funcf(xin, uin, t)
    % Extract values
    x = xin(1);
    y = xin(2);
    z = xin(3);
    u = xin(4);
    v = xin(5);
    w = xin(6);
    phi = xin(7);
    theta = xin(8);
    psi = xin(9);
    v_wxE = xin(10);
    v_wyE = xin(11);
    v_wzE = xin(12);

    b_A_x = xin(13);
    b_A_y = xin(14);
    b_A_z = xin(15);
    b_p = xin(16);
    b_q = xin(17);
    b_r = xin(18);

    A_x = uin(1);
    A_y = uin(2);
    A_z = uin(3);
    p = uin(4);
    q = uin(5);
    r = uin(6);
    
    g = 9.81;
    Ts = 0.01;

    omega_b_A_x = 0;
    omega_b_A_y = 0;
    omega_b_A_z = 0;
    omega_b_p = 0;
    omega_b_q = 0;
    omega_b_r = 0;

    % original continuous-time nonlinear system equations
    x_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
    y_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
    z_dot =-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
    
    u_dot =(A_x-b_A_x)-g*sin(theta)+(r-b_r)*v-(q-b_q)*w;
    v_dot =(A_y-b_A_y)+g*cos(theta)*sin(phi)+(p-b_p)*w-(r-b_r)*u;
    w_dot =(A_z-b_A_z)+g*cos(theta)*cos(phi)+(q-b_q)*u-(p-b_p)*v;
    
    phi_dot =(p-b_p)+(q-b_q)*sin(phi)*tan(theta)+(r-b_r)*cos(phi)*tan(theta);
    theta_dot =(q-b_q)*cos(phi)-(r-b_r)*sin(phi);
    psi_dot =(q-b_q)*sin(phi)/cos(theta)+(r-b_r)*cos(phi)/cos(theta);
     
    v_wxE_dot =0;
    v_wyE_dot =0;
    v_wzE_dot =0;  
    
    b_A_x_dot = omega_b_A_x/Ts;
    b_A_y_dot = omega_b_A_y/Ts;
    b_A_z_dot = omega_b_A_z/Ts;
    b_p_dot = omega_b_p/Ts;
    b_q_dot = omega_b_q/Ts;
    b_r_dot = omega_b_r/Ts;

    % original, continuous-time, noiseless nonlinear dynamics
    x_dot = [x_dot; y_dot; z_dot; u_dot; v_dot; w_dot; phi_dot; theta_dot; psi_dot; v_wxE_dot; v_wyE_dot; v_wzE_dot; b_A_x_dot; b_A_y_dot; b_A_z_dot; b_p_dot; b_q_dot; b_r_dot];
end

function d = funch(xin,uin,t)
    % Extract values
    x = xin(1);
    y = xin(2);
    z = xin(3);
    u = xin(4);
    v = xin(5);
    w = xin(6);
    phi = xin(7);
    theta = xin(8);
    psi = xin(9);
    v_wxE = xin(10);
    v_wyE = xin(11);
    v_wzE = xin(12);

    b_A_x = xin(13);
    b_A_y = xin(14);
    b_A_z = xin(15);
    b_p = xin(16);
    b_q = xin(17);
    b_r = xin(18);

    A_x = uin(1);
    A_y = uin(2);
    A_z = uin(3);
    p = uin(4);
    q = uin(5);
    r = uin(6);

    g = 9.81;

    % original continuous-time nonlinear system equations
    x_GPS = x;
    y_GPS = y;
    z_GPS = z;
        
    u_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
    v_GPS = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
    w_GPS = -u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
        
    phi_GPS = phi;
    theta_GPS = theta;
    psi_GPS = psi;
    
    V = sqrt(u^2+v^2+w^2);
    alpha = atan(w/u);
    beta = atan(v/sqrt(u^2+w^2));
    
    % Final observation dynamics
    d = [x_GPS, y_GPS, z_GPS, u_GPS, v_GPS, w_GPS, phi_GPS, theta_GPS, psi_GPS, V, alpha, beta];
end

%% lineaized dynamics in discrete time
function [Phi,Gamma] = funcLinDisDyn(x_linpt,u_linpt,Ts)
    % Extract values
    x = x_linpt(1);
    y = x_linpt(2);
    z = x_linpt(3);
    u = x_linpt(4);
    v = x_linpt(5);
    w = x_linpt(6);
    phi = x_linpt(7);
    theta = x_linpt(8);
    psi = x_linpt(9);
    v_wxE = x_linpt(10);
    v_wyE = x_linpt(11);
    v_wzE = x_linpt(12);

    b_A_x = x_linpt(13);
    b_A_y = x_linpt(14);
    b_A_z = x_linpt(15);
    b_p = x_linpt(16);
    b_q = x_linpt(17);
    b_r = x_linpt(18);

    A_x_m = u_linpt(1);
    A_y_m = u_linpt(2);
    A_z_m = u_linpt(3);
    p_m = u_linpt(4);
    q_m = u_linpt(5);
    r_m = u_linpt(6);

    omega_A_x = 0;
    omega_A_y = 0;
    omega_A_z = 0;
    omega_p = 0;
    omega_q = 0;
    omega_r = 0;

    % F/A Matrix
    F = [[0, 0, 0, cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)),                                   -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 1, 0, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0, cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)),                                   -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 1, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,         -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),                                             - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 1,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,                   0,                                        r_m - b_r,                                        b_q - q_m,                                                                                  0,                                                                             -(981*cos(theta))/100,                                                                                                     0, 0, 0, 0, -1,  0,  0,  0,                    w,                   -v]
        [0, 0, 0,           b_r - r_m,                                                0,                                        p_m - b_p,                                                      (981*cos(phi)*cos(theta))/100,                                                                    -(981*sin(phi)*sin(theta))/100,                                                                                                     0, 0, 0, 0,  0, -1,  0, -w,                    0,                    u]
        [0, 0, 0,           q_m - b_q,                                        b_p - p_m,                                                0,                                                     -(981*cos(theta)*sin(phi))/100,                                                                    -(981*cos(phi)*sin(theta))/100,                                                                                                     0, 0, 0, 0,  0,  0, -1,  v,                   -u,                    0]
        [0, 0, 0,                   0,                                                0,                                                0,                  sin(phi)*tan(theta)*(b_r - r_m) - cos(phi)*tan(theta)*(b_q - q_m),               - cos(phi)*(b_r - r_m)*(tan(theta)^2 + 1) - sin(phi)*(b_q - q_m)*(tan(theta)^2 + 1),                                                                                                     0, 0, 0, 0,  0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta)]
        [0, 0, 0,                   0,                                                0,                                                0,                                        cos(phi)*(b_r - r_m) + sin(phi)*(b_q - q_m),                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,            -cos(phi),             sin(phi)]
        [0, 0, 0,                   0,                                                0,                                                0,              (sin(phi)*(b_r - r_m))/cos(theta) - (cos(phi)*(b_q - q_m))/cos(theta), - (cos(phi)*sin(theta)*(b_r - r_m))/cos(theta)^2 - (sin(phi)*sin(theta)*(b_q - q_m))/cos(theta)^2,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta)]
        [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]
        [0, 0, 0,                   0,                                                0,                                                0,                                                                                  0,                                                                                                 0,                                                                                                     0, 0, 0, 0,  0,  0,  0,  0,                    0,                    0]];
     

    % G Matrix
    G = [[ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
        [-1,  0,  0,  0,                    w,                   -v,   0,   0,   0,   0,   0,   0]
        [ 0, -1,  0, -w,                    0,                    u,   0,   0,   0,   0,   0,   0]
        [ 0,  0, -1,  v,                   -u,                    0,   0,   0,   0,   0,   0,   0]
        [ 0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta),   0,   0,   0,   0,   0,   0]
        [ 0,  0,  0,  0,            -cos(phi),             sin(phi),   0,   0,   0,   0,   0,   0]
        [ 0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta),   0,   0,   0,   0,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0, 100,   0,   0,   0,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0, 100,   0,   0,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0,   0, 100,   0,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0, 100,   0,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0, 100,   0]
        [ 0,  0,  0,  0,                    0,                    0,   0,   0,   0,   0,   0, 100]];
         
    [Phi,Gamma]=c2d(F,G,Ts);
end

function H = funcLinDisObs(x_linpt,u_linpt,t)
    % Extract values
    x = x_linpt(1);
    y = x_linpt(2);
    z = x_linpt(3);
    u = x_linpt(4);
    v = x_linpt(5);
    w = x_linpt(6);
    phi = x_linpt(7);
    theta = x_linpt(8);
    psi = x_linpt(9);
    v_wxE = x_linpt(10);
    v_wyE = x_linpt(11);
    v_wzE = x_linpt(12);

    b_A_x = x_linpt(13);
    b_A_y = x_linpt(14);
    b_A_z = x_linpt(15);
    b_p = x_linpt(16);
    b_q = x_linpt(17);
    b_r = x_linpt(18);

    A_x = u_linpt(1);
    A_y = u_linpt(2);
    A_z = u_linpt(3);
    p = u_linpt(4);
    q = u_linpt(5);
    r = u_linpt(6);

    H = [[1, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        [0, 1, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        [0, 0, 1,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        [0, 0, 0,                              cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), sin(psi)*(w*cos(phi) + v*sin(phi)) + cos(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)), -cos(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psi)*(v*cos(phi) - w*sin(phi)), 1, 0, 0, 0, 0, 0, 0, 0, 0]
        [0, 0, 0,                              cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), sin(psi)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psi)*(w*cos(phi) + v*sin(phi)), -sin(psi)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psi)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psi)*(v*cos(phi) - w*sin(phi)), 0, 1, 0, 0, 0, 0, 0, 0, 0]
        [0, 0, 0,                                      -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta),                                               cos(theta)*(v*cos(phi) - w*sin(phi)),           - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                     0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
        [0, 0, 0,                                                0,                                                0,                                                0,                                                                                  1,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        [0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               1,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        [0, 0, 0,                                                0,                                                0,                                                0,                                                                                  0,                                                               0,                                                                                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        [0, 0, 0,                        u/(u^2 + v^2 + w^2)^(1/2),                        v/(u^2 + v^2 + w^2)^(1/2),                        w/(u^2 + v^2 + w^2)^(1/2),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        [0, 0, 0,                           -w/(u^2*(w^2/u^2 + 1)),                                                0,                              1/(u*(w^2/u^2 + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        [0, 0, 0, -(u*v)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),      1/((u^2 + w^2)^(1/2)*(v^2/(u^2 + w^2) + 1)), -(v*w)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),                                                                                  0,                                                               0,                                                                                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
end
