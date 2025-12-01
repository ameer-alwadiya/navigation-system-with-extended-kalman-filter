clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                          CAUTION                                            %%
%% This test will only help you spot syntex errors and check for dimension mismatches.         %%
%% Passing the test does not mean your script will run error free with another flight dataset. %%
%% It is YOUR responsibility to design your script to ROBUSTLY handle different scenarios.     %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load data to test your design
% load('dataTask3_1.mat') %with only input measurement sensor faults
% load('dataTask3_2.mat') %with only angle of attack sensor fault
load('Data/dataTask4.mat') %with both input measurement sensor and angle of attack sensor faults. 

%Test your script, remember to change the function name according to your student ID (as on Blackboard)
tic;
[x_est,b_est,Ax_f_instance,Ay_f_instance,Az_f_instance,p_f_instance,q_f_instance,r_f_instance,AoA_f_instance, Innov] = SID230255137(c_k, d_k, t, dt);
tc=toc;

% Please make sure your code finishes running within 5 minutes. 
disp(['Total run time is ',num2str(tc),' seconds'])

if size(x_est,1)==length(t) && size(x_est,2)==12 && ~any(isnan(x_est),'all')
    disp('Dimesion of x_est returned by the function is in good order')
else
    error('Dimesion of x_est returned by the function is not correct or x_est contains NaN elements')
end

if size(b_est,1)==length(t) && size(b_est,2)==6 && ~any(isnan(b_est),'all')
    disp('Dimesion of b_est returned by the function is in good order')
else
    error('Dimesion of b_est returned by the function is not correct or b_est contains NaN elements')
end


N = 8001;
figure;
for i = 1:12
    subplot(2, 6, i)
    plot(t(1:N), x_est(1:N,i), 'b')
    grid on
end

figure;
for i = 1:6
    subplot(2, 3, i)
    plot(t(1:N), b_est(1:N,i), 'b')
    grid on
end

%% fault detection
figure;
subplot(2, 3, 1)
plot(t(1:8001), b_est(1:8001,1), 'b')
hold on
grid on
plot(t(Ax_f_instance), b_est(Ax_f_instance, 1), 'rx', 'MarkerSize', 20) % Plotting red 'x' markers at fault instances
xlabel('Time')
ylabel('b_est')
legend('b_est', 'Fault Instances')
hold off

subplot(2, 3, 2)
plot(t(1:8001), b_est(1:8001,2), 'b')
hold on
grid on
plot(t(Ay_f_instance), b_est(Ay_f_instance, 2), 'rx', 'MarkerSize', 20) % Plotting red 'x' markers at fault instances
xlabel('Time')
ylabel('b_est')
legend('b_est', 'Fault Instances')
hold off

subplot(2, 3, 3)
plot(t(1:8001), b_est(1:8001,3), 'b')
hold on
grid on
plot(t(Az_f_instance), b_est(Az_f_instance, 3), 'rx', 'MarkerSize', 20) % Plotting red 'x' markers at fault instances
xlabel('Time')
ylabel('b_est')
legend('b_est', 'Fault Instances')
hold off

subplot(2, 3, 4)
plot(t(1:8001), b_est(1:8001,4), 'b')
hold on
grid on
plot(t(p_f_instance), b_est(p_f_instance, 4), 'rx', 'MarkerSize', 20) % Plotting red 'x' markers at fault instances
xlabel('Time')
ylabel('b_est')
legend('b_est', 'Fault Instances')
hold off

subplot(2, 3, 5)
plot(t(1:8001), b_est(1:8001,5), 'b')
hold on
grid on
plot(t(q_f_instance), b_est(q_f_instance, 5), 'rx', 'MarkerSize', 20) % Plotting red 'x' markers at fault instances
xlabel('Time')
ylabel('b_est')
legend('b_est', 'Fault Instances')
hold off

subplot(2, 3, 6)
plot(t(1:8001), b_est(1:8001,6), 'b')
hold on
grid on
if numel(r_f_instance) == 1 && r_f_instance(1) == 0
    disp('no faults');
else
    plot(t(r_f_instance), b_est(r_f_instance, 6), 'rx', 'MarkerSize', 20) % Plotting red 'x' markers at fault instances
end
xlabel('Time')
ylabel('b_est')
legend('b_est', 'Fault Instances')
hold off

%% 
angle_of_attack = Innov(1:N,11);
figure;
plot(t(1:N), angle_of_attack, 'b')
title('Innovation of Alpha - Angle of Attack')
hold on
grid on
if numel(AoA_f_instance) == 1 && AoA_f_instance(1) == 0
    disp('no faults');
else
    plot(t(AoA_f_instance), angle_of_attack(AoA_f_instance), 'rx', 'MarkerSize', 20) % Plotting red 'x' markers at fault instances
end
xlabel('Time')
ylabel('b_est')
legend('b_est', 'Fault Instances')
hold off
grid on