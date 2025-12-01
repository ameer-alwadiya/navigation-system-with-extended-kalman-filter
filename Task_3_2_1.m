clear variables

%% Kinematic Model
% Read data files
load('Data/dataTask3_1.mat');

% define symbolic variables
syms x y z u v w phi theta psi v_wxE v_wyE v_wzE b_A_x b_A_y b_A_z b_p b_q b_r
syms A_x A_y A_z p q r

% including the additive measurement noise dynamics
syms omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r
syms A_x_m A_y_m A_z_m p_m q_m r_m 

syms omega_b_A_x omega_b_A_y omega_b_A_z omega_b_p omega_b_q omega_b_r
     
g = 9.81;
Ts=dt;

% original continuous-time nonlinear system equations
x_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
y_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
z_dot =-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
    
u_dot =A_x-g*sin(theta)+r*v-q*w;
v_dot =A_y+g*cos(theta)*sin(phi)+p*w-r*u;
w_dot =A_z+g*cos(theta)*cos(phi)+q*u-p*v;
    
phi_dot =p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);
theta_dot =q*cos(phi)-r*sin(phi);
psi_dot =q*sin(phi)/cos(theta)+r*cos(phi)/cos(theta);
    
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
x_dot = [x_dot y_dot z_dot u_dot v_dot w_dot phi_dot theta_dot psi_dot v_wxE_dot v_wyE_dot v_wzE_dot b_A_x_dot b_A_y_dot b_A_z_dot b_p_dot b_q_dot b_r_dot];
     
% substitute the input measurement model into dynamics 
xm_dot = subs(x_dot, [A_x A_y A_z p q r], [A_x_m-omega_A_x-b_A_x, A_y_m-omega_A_y-b_A_y, A_z_m-omega_A_z-b_A_z, p_m-omega_p-b_p, q_m-omega_q-b_q, r_m-omega_r-b_r]);

%% Observation Model

x_GPS=x;
y_GPS=y;
z_GPS=z;
    
u_GPS=(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+v_wxE;
v_GPS=(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+v_wyE;
w_GPS=-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
    
phi_GPS=phi;
theta_GPS=theta;
psi_GPS=psi;

V=sqrt(u^2+v^2+w^2);
alpha=atan(w/u);
beta=atan(v/sqrt(u^2+w^2));

% Final observation dynamics
d = [x_GPS y_GPS z_GPS u_GPS v_GPS w_GPS phi_GPS theta_GPS psi_GPS V alpha beta];

%% Linearisation of the system 
f_expression = subs(xm_dot, [omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r omega_b_A_x omega_b_A_y omega_b_A_z omega_b_p omega_b_q omega_b_r], [0 0 0 0 0 0 0 0 0 0 0 0]);
h_expression = d;

% Jacobian linearization to obtain the expression for F | A
F = jacobian (f_expression, [x y z u v w phi theta psi v_wxE v_wyE v_wzE b_A_x b_A_y b_A_z b_p b_q b_r]);

% Jacobian linearization to obtain the expression for B (not required)
% B = jacobian (f_expression, [A_x_m A_y_m A_z_m p_m q_m r_m]);

% Jacobian linearization to obtain the expression for G 
G = jacobian (xm_dot, [omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r omega_b_A_x omega_b_A_y omega_b_A_z omega_b_p omega_b_q omega_b_r]);

% Jacobian linearization to obtain the expression for H
H = jacobian (h_expression, [x y z u v w phi theta psi v_wxE v_wyE v_wzE b_A_x b_A_y b_A_z b_p b_q b_r]);



     



