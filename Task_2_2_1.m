clear variables

%% Kinematic Model

% define symbolic variables
syms x y z u v w phi theta psi v_wxE v_wyE v_wzE b_A_x b_A_y b_A_z b_p b_q b_r
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

b_A_x_dot=0;
b_A_y_dot=0;
b_A_z_dot=0;
b_p_dot=0;
b_q_dot=0;
b_r_dot=0;

% original, continuous-time, noiseless nonlinear dynamics
X_dot = [x_dot y_dot z_dot u_dot v_dot w_dot phi_dot theta_dot psi_dot v_wxE_dot v_wyE_dot v_wzE_dot b_A_x_dot b_A_y_dot b_A_z_dot b_p_dot b_q_dot b_r_dot];
     
% substitute the input measurement model into dynamics 
Xm_dot = subs(X_dot, [A_x A_y A_z p q r], [A_x_m-omega_A_x-b_A_x, A_y_m-omega_A_y-b_A_y, A_z_m-omega_A_z-b_A_z, p_m-omega_p-b_p, q_m-omega_q-b_q, r_m-omega_r-b_r]);

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
f_expression = subs(Xm_dot, [omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r], [0 0 0 0 0 0]);
h_expression = d;

X = [x y z u v w phi theta psi v_wxE v_wyE v_wzE b_A_x b_A_y b_A_z b_p b_q b_r];
U = [A_x_m A_y_m A_z_m p_m q_m r_m];
W = [omega_A_x omega_A_y omega_A_z omega_p omega_q omega_r];

% Jacobian linearization to obtain the expression for F | A
F = jacobian (f_expression, X);

% Jacobian linearization to obtain the expression for B (not required)
% B = jacobian (f_expression, [A_x_m A_y_m A_z_m p_m q_m r_m]);

% Jacobian linearization to obtain the expression for G 
G = jacobian (Xm_dot, W);

% Jacobian linearization to obtain the expression for H
H = jacobian (h_expression, X);



     



