clear variables

%% Kinematic Model (state equations)

% define symbolic variables
syms x y z u v w phi theta psi v_wxE v_wyE v_wzE 
syms A_x A_y A_z p q r

% including the additive measurement noise dynamics
syms omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r
syms A_x_m A_y_m A_z_m p_m q_m r_m 
g = 9.81;

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

% original, continuous-time, noiseless nonlinear dynamics
X_dot = [x_dot y_dot z_dot u_dot v_dot w_dot phi_dot theta_dot psi_dot v_wxE_dot v_wyE_dot v_wzE_dot];

% substitute the input measurement model into dynamics 
Xm_dot = subs(X_dot, [A_x A_y A_z p q r], [A_x_m-omega_A_x A_y_m-omega_A_y A_z_m-omega_A_z p_m-omega_p q_m-omega_q r_m-omega_r]);

%% Observation  (measurement equations)

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
X = [x y z u v w phi theta psi v_wxE v_wyE v_wzE];
U = [A_x_m A_y_m A_z_m p_m q_m r_m];
W = [omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r];

f_expression = subs(Xm_dot, [omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r], [0 0 0 0 0 0]);
h_expression = d;

% Jacobian linearization to obtain the expression for F | A
F = jacobian (f_expression, X);

% Jacobian linearization to obtain the expression for B (not required)
B = jacobian (f_expression, U);

% Jacobian linearization to obtain the expression for G 
G = jacobian (Xm_dot, W);

% Jacobian linearization to obtain the expression for H
H = jacobian (h_expression, X);

%% task 1.1: Extracting the requested elements

% The 5th row, 4th column of the F matrix
F_5_4 = F(5, 4);

% The 6th row, 7th column of the F matrix
F_6_7 = F(6, 7);

% The 4th row, 6th column of the G matrix
G_4_6 = G(4, 6);

% The 9th row, 5th column of the G matrix
G_9_5 = G(9, 5);

% The 7th row, 7th column of the H matrix
H_7_7 = H(7, 7);

% The 11th row, 4th column of the H matrix
H_11_4 = H(11, 4);

% Displaying the extracted elements
disp(['F(5, 4) = ', char(F_5_4)]);
disp(['F(6, 7) = ', char(F_6_7)]);
disp(['G(4, 6) = ', char(G_4_6)]);
disp(['G(9, 5) = ', char(G_9_5)]);
disp(['H(7, 7) = ', char(H_7_7)]);
disp(['H(11, 4) = ', char(H_11_4)]);