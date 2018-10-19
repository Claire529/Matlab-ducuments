% backstepping with fuzzy based on quaternion representation
function [dx,w,q,angle,tao]=Quaternion_fuzzy_backstepping(t,x,k1,k2,k3,k4,k5,k6,c_q1,ce_q1,w_c_q1,w_ce_q1,c_q2,ce_q2,w_c_q2,w_ce_q2,c_q3,ce_q3,w_c_q3,w_ce_q3)

w_1 = x(1,:);
w_2 = x(2,:);
w_3 = x(3,:);
w = [w_1;w_2;w_3];
ww=size(w);
q_w=[0;w];

q_0 = x(4,:);
q_1 = x(5,:);
q_2 = x(6,:);
q_3 = x(7,:);
q = [q_0;q_1;q_2;q_3];

w_1_hat=x(8:16);
w_2_hat=x(17:25);
w_3_hat=x(26:34);

[theta phi gamma]=quat2angle(q');
angle=[theta phi gamma];
q_dot=1/2*quatmultiply(q',q_w');
q_dot=q_dot';
q_0_dot = q_dot(1,1);
q_dot(1,1)  = q_dot(2,1);
q_dot(2,1)  = q_dot(3,1);
q_dot(3,1)  = q_dot(4,1);

rho=1.2;
V=5;
Sb=0.005;
l1=0.07;
l2=0.0125;
alpha=5*pi/180;
mx=0.0117e-7;
my=0.0118e-7;
mz=0.0181e-7;
Jx=4e-7;
Jy=7e-7;
Jz=1.1e-7;

Mbx = rho/2*V*V*Sb*l1*(mx*cos(alpha)+my*sin(alpha));
Mby = rho/2*V*V*Sb*l1*(-mx*sin(alpha)+my*cos(alpha));
Mbz = rho/2*V*V*mz*Sb*l2;

q_1_d=0;
q_2_d=0;
q_3_d=0;
q_0_d=sqrt(1-q_1_d^2-q_2_d^2-q_3_d^2);

q_d_inv=[q_0_d;-q_1_d;-q_2_d;-q_3_d];
q_e=quatmultiply(q',q_d_inv');
q_e=q_e';
q_0_e=q_e(1,1);
q_1_e=q_e(2,1);
q_2_e=q_e(3,1);
q_3_e=q_e(4,1);

q_0_e_dot=0.5*(-w_1*q_1_e-w_2*q_2_e-w_3*q_3_e);
q_1_e_dot=0.5*(w_1*q_0-w_2*q_3_e+w_3*q_2_e);
q_2_e_dot=0.5*(w_2*q_0-w_3*q_2_e+w_1*q_3_e);
q_3_e_dot=0.5*(w_3*q_0-w_1*q_2_e+w_2*q_1_e);


w_1_star=-(k1*q_1_e);
w_2_star=-(k2*q_2_e);
w_3_star=-(k3*q_3_e);


w_1_star_dot=-(k1*q_1_e_dot);
w_2_star_dot=-(k2*q_2_e_dot);
w_3_star_dot=-(k3*q_3_e_dot);

s_1=Quaternion_calCtrl_fuz_adap(w_2,w_3,c_q1,ce_q1,w_c_q1,w_ce_q1);
s_2=Quaternion_calCtrl_fuz_adap(w_1,w_3,c_q2,ce_q2,w_c_q2,w_ce_q2);
s_3=Quaternion_calCtrl_fuz_adap(w_1,w_2,c_q3,ce_q3,w_c_q3,w_ce_q3);

Gamma_1=10;
Gamma_2=10;
Gamma_3=10;
delta_1=10;
delta_2=10;
delta_3=10;
w_1_hat_dot=s_1'*(w_1-w_1_star)*Gamma_1-w_1_hat*delta_1;
w_2_hat_dot=s_2'*(w_2-w_2_star)*Gamma_2-w_2_hat*delta_2;
w_3_hat_dot=s_3'*(w_3-w_3_star)*Gamma_3-w_3_hat*delta_3;

tao_1=-q_1_e-w_1_hat'*s_1'-Mbx/Jx+((Jz-Jy)/Jx)*w_3*w_2+w_1_star_dot-(k4+0.5)*(w_1-w_1_star);
tao_2=-q_2_e-w_2_hat'*s_2'-Mby/Jy+((Jx-Jz)/Jy)*w_1*w_2+w_2_star_dot-(k5+0.5)*(w_2-w_2_star);
tao_3=-q_3_e-w_3_hat'*s_3'-Mbz/Jz+((Jy-Jx)/Jz)*w_2*w_1+w_3_star_dot-(k6+0.5)*(w_3-w_3_star);
tao=[tao_1,tao_2,tao_3];

DELTA_1=1;
DELTA_2=1;
DELTA_3=1;
w_dot=[Mbx/Jx-((Jz-Jy)/Jx)*w_3*w_2+s_1*w_1_hat+DELTA_1+tao_1;
       Mby/Jy-((Jx-Jz)/Jy)*w_1*w_3+s_2*w_2_hat+DELTA_2+tao_2;
       Mbz/Jz-((Jy-Jx)/Jz)*w_2*w_1+s_3*w_3_hat+DELTA_3+tao_3];

f1 = w_dot(1,1);
f2 = w_dot(2,1);
f3 = w_dot(3,1);
f4 = q_0_dot;
f5 = q_dot(1,1);
f6 = q_dot(2,1);
f7 = q_dot(3,1);
dx = [f1;f2;f3;f4;f5;f6;f7;w_1_hat_dot;w_2_hat_dot;w_3_hat_dot];
