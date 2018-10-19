function [Mxyz,Pxyz,dx]=attitude_angle_lxp(t,x,tau)
theta = x(1,:);%pitch
phi = x(2,:);%yaw
gamma = x(3,:);%roll
status_x1=[theta;phi;gamma];
theta_dot = x(4,:);
phi_dot = x(5,:);
gamma_dot = x(6,:);
% p = x(4,:);
% q = x(5,:);
% r = x(6,:);
status_x2=[theta_dot;phi_dot;gamma_dot];

mx=0.0117e-7;
my=0.0118e-7;
mz=0.0181e-7;
Jx=4e-7;
Jy=7e-7;
Jz=1.1e-7;

% mx=0.0117;
% my=0.0118;
% mz=0.0181;
% Jx=4;
% Jy=7;
% Jz=1.1;

rho=1.2;
R=0.07;
c_bar=0.05;
Sb=0.003;
l1=0.07;
l2=0.0125;
yita=0.75;
V=5;
c=0.05;
R=0.07;
xi=0.25*c;
zeta=0.7*R;
alpha1=5*pi/180;
beita=0;

eta_l=0;
eta_r=0;
beta_l=0;
beta_r=0;

% pexi_l=(30*pi/180)*sin(50*pi*t-10*pi/180);
% pexi_r=((-30*pi/180)*sin(50*pi*t-10*pi/180));
% phi_l=((-60*pi/180)*cos(50*pi*t)-20*pi/180);
% phi_r=(60*pi/180)*cos(50*pi*t)+20*pi/180;
% pexi_l_dot=(30*pi/180)*50*pi*cos(50*pi*t-10*pi/180);
% pexi_r_dot=((-30*pi/180)*50*pi*cos(50*pi*t-10*pi/180));
% phi_l_dot=(-60*pi/180)*50*pi*sin(50*pi*t);
% phi_r_dot=(60*pi/180)*50*pi*sin(50*pi*t);
% 
% Mbx = rho/2*V*V*Sb*l1*(mx*cos(alpha1)+my*sin(alpha1));
% Mby = rho/2*V*V*Sb*l1*(-mx*sin(alpha1)+my*cos(alpha1));
% Mbz = rho/2*V*V*mz*Sb*l2;


f=25;

alpha_w=5*pi/180;%*sin(2*pi*f*t+0*pi/2);

omega=2*pi*f;
Psi_l=30*pi/180;
Psi_r=-30*pi/180;
Phi_l=60*pi/180;
Phi_r=-60*pi/180;
d_phi_l=10*pi/180;
d_phi_r=-10*pi/180;
lambda=10*pi/180;
pexi_l=Psi_l*sin(omega*t-lambda);
pexi_r=Psi_r*sin(omega*t-lambda);
phi_l=Phi_l*cos(omega*t)+d_phi_l;
phi_r=Phi_r*cos(omega*t)+d_phi_r;

pexi_l_dot=omega*Psi_l*cos(omega*t-lambda);
pexi_r_dot=omega*Psi_r*cos(omega*t-lambda);
phi_l_dot=-omega*Phi_l*sin(omega*t);
phi_r_dot=-omega*Phi_r*sin(omega*t);

Mbx = (rho/2)*V*V*Sb*l1*(mx*cos(alpha1)+my*sin(alpha1));
Mby = (rho/2)*V*V*Sb*l1*(-mx*sin(alpha1)+my*cos(alpha1));
Mbz = (rho/2)*V*V*Sb*l2*mz;%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mbx = (rho/2)*V*V*Sb*l1*mx;
% Mby = (rho/2)*V*V*Sb*l1*my;
Mbz = (rho/2)*V*V*Sb*l1*mz;

theta_bolang=0.3*pi/180*sin(pi*t/3);
phi_bolang=0.5*pi/180*sin(pi*t/6);
gamma_bolang=0.2*pi/180*sin(pi*t/4);

Mdx=0.1*Jz*pi/180;
Mdy=0.2*Jy*pi/180;
Mdz=0.3*Jx*pi/180;
Mux=0.2*Jz*pi/180;
Muy=0.15*Jy*pi/180;
Muz=0.1*Jx*pi/180;

C_l=0.225+1.58*sin(2.13*alpha_w-7.2*pi/180);
C_d=1.92-1.55*cos(2.04*alpha_w-9.82*pi/180);
C_rot=pi*(0.75-0.25);

F_L_l=0.5*rho*C_l*phi_l_dot*phi_l_dot*c_bar*R*R*R/3;
No=floor(omega*t/(2*pi));
remainder=omega*t-No*2*pi;
if (remainder<pi)
    F_L_l=F_L_l;
else
    F_L_l=-F_L_l;
end
F_D_l=0.5*rho*C_d*phi_l_dot*phi_l_dot*c_bar*R*R*R/3;

F_L_r=0.5*rho*C_l*phi_r_dot*phi_r_dot*c_bar*R*R*R/3;
if (remainder<pi)
    F_L_r=F_L_r;
else
    F_L_r=-F_L_r;
end
F_D_r=0.5*rho*C_d*phi_r_dot*phi_r_dot*c_bar*R*R*R/3;

Frot_l=rho*C_rot*phi_l_dot*pexi_l_dot*R*c_bar*c_bar*c_bar*c_bar*R*R*R*R/4;
Frot_r=rho*C_rot*phi_r_dot*pexi_r_dot*R*c_bar*c_bar*c_bar*c_bar*R*R*R*R/4;  

F_x_l= F_D_l;
F_y_l=-F_L_l*sin(phi_l);
F_z_l= F_L_l*cos(phi_l);
Frot_x_l=Frot_l;
Frot_y_l=-Frot_l*sin(phi_l);
Frot_z_l=Frot_l*cos(phi_l);
Frot_xyz_l=[Frot_x_l;Frot_y_l;Frot_z_l];
F_xyz_l=[F_x_l;F_y_l;F_z_l];
F_x_r= F_D_r;
F_y_r= F_L_r*sin(phi_r);
F_z_r=-F_L_r*cos(phi_r);
Frot_x_r=Frot_r;
Frot_y_r=Frot_r*sin(phi_r);
Frot_z_r=-Frot_r*cos(phi_r);
Frot_xyz_r=[Frot_x_r;Frot_y_r;Frot_z_r];
F_xyz_r=[F_x_r;F_y_r;F_z_r];
L_l=[ cos(beta_l)*cos(eta_l)        sin(eta_l)   sin(beta_l)*cos(eta_l)
     -cos(beta_l)*sin(eta_l)        cos(eta_l)  -sin(beta_l)*sin(eta_l)
     -sin(beta_l)                   0            cos(beta_l)];
L_r=[ cos(beta_r)*cos(eta_r)        sin(eta_r)  -sin(beta_r)*cos(eta_r)
      cos(beta_r)*sin(eta_r)       -cos(eta_r)  -sin(beta_r)*sin(eta_r)
     -sin(beta_r)                   0           -cos(beta_r)];

Pxyz=L_l*F_xyz_l+L_l*Frot_xyz_l+L_r*F_xyz_r+L_r*Frot_xyz_r;
% Pxyz=L_l*F_xyz_l+L_r*F_xyz_r;

v_b=[V*cos(alpha1)*cos(beita)
    -V*sin(alpha1)*cos(beita)
     V*sin(beita)               ];
v_l=[-xi*pexi_l_dot*sin(pexi_l)
     -xi*pexi_l_dot*cos(pexi_l)*cos(phi_l)+xi*phi_l_dot*sin(pexi_l)*sin(phi_l)-zeta*phi_l_dot*cos(phi_l)
      xi*pexi_l_dot*cos(pexi_l)*sin(phi_l)+xi*phi_l_dot*sin(pexi_l)*cos(phi_l)+zeta*phi_l_dot*sin(phi_l)];

v_r=[-xi*pexi_r_dot*sin(pexi_r)
      xi*pexi_r_dot*cos(pexi_r)*cos(phi_r)-xi*phi_r_dot*sin(pexi_r)*sin(phi_r)+zeta*phi_r_dot*cos(phi_r)
     -xi*pexi_r_dot*cos(pexi_r)*sin(phi_r)-xi*phi_r_dot*sin(pexi_r)*cos(phi_r)-zeta*phi_r_dot*sin(phi_r)];
u_l=v_l+v_b;
u_r=v_r+v_b;

% v_l_x=;
% v_l_y=;
% v_l_z=+;
% 
% v_r_x=-a*pexi_r_dot*sin(pexi_r)+V*cos(alpha1)*cos(beita);
% v_r_y=-V*sin(alpha1)*cos(beita)+a*pexi_r_dot*cos(pexi_r)*cos(phi_r)-a*phi_r_dot*sin(pexi_r)*sin(phi_r)+b*phi_r_dot*cos(phi_r);
% v_r_z=V*sin(beita)-a*pexi_r_dot*cos(pexi_r)*sin(phi_r)-a*phi_r_dot*sin(pexi_r)*cos(phi_r)-b*phi_r_dot*sin(phi_r);
% 
% FLlx=rho*C_l*v_l_x*v_l_x*0.5*c*R;
% FDlx=rho*C_d*v_l_x*v_l_x*0.5*c*R;
% Frotlx=rho*C_rot*v_l_x*w_l_x*c_bar*c_bar*R*0.5*c_bar*c_bar*R*R;
% FLrx=rho*C_l*v_r_x*v_r_x*0.5*c*R;
% FDrx=rho*C_d*v_r_x*v_r_x*0.5*c*R;
% Frotrx=rho*C_rot*v_r_x*w_r_x*c_bar*c_bar*R*4*c_bar*c_bar*0.5*R*R;
% FLly=rho*C_l*v_l_y*v_l_y*0.5*c*R;
% FDly=rho*C_d*v_l_y*v_l_y*0.5*c*R;
% Frotly=rho*C_rot*v_l_y*w_l_y*c_bar*c_bar*R*0.5*c_bar*c_bar*R*R;
% FLry=rho*C_l*v_r_y*v_r_y*0.5*c*R;
% FDry=rho*C_d*v_r_y*v_r_y*0.5*c*R;
% Frotry=rho*C_rot*v_r_y*w_r_y*c_bar*c_bar*R*0.5*c_bar*c_bar*R*R;
% FLlz=rho*C_l*v_l_z*v_l_z*0.5*c*R;
% FDlz=rho*C_d*v_l_z*v_l_z*0.5*c*R;
% Frotlz=rho*C_rot*v_l_z*w_l_z*c_bar*c_bar*R*0.5*c_bar*c_bar*R*R;
% FLrz=rho*C_l*v_r_z*v_r_z*0.5*c*R;
% FDrz=rho*C_d*v_r_z*v_r_z*0.5*c*R;
% Frotrz=rho*C_rot*v_r_z*w_r_z*c_bar*c_bar*R*0.5*c_bar*c_bar*R*R;
% Frx=FLrx-FDrx+Frotrx;
% Fry=FLry+FDry+Frotry;
% Frz=FLrz+FDrz+Frotrz;
% Flx=FLlx-FDlx+Frotlx;
% Fly=FLly+FDly+Frotly;
% Flz=FLlz+FDlz+Frotlz;
% Px=Frx+Flx;
% Py=Fry+Fly;
% Pz=Frz+Flz;
% Mwx=Fry*(a*sin(pexi_r)*sin(phi_r)-b*cos(phi_r)-yita)-Fly*(a*sin(pexi_l)*sin(phi_l)-b*cos(phi_l)-yita)+Frz*(a*sin(pexi_r)*cos(phi_r)+b*sin(phi_r))-Flz*(a*sin(pexi_l)*cos(phi_l)+b*sin(phi_l));
% Mwy=Frx*(-a*sin(pexi_r)*sin(phi_r)+b*cos(phi_r)+yita)+Flx*(a*sin(pexi_l)*sin(phi_l)-b*cos(phi_l)-yita)-Frz*a*cos(pexi_r)+Flz*a*cos(pexi_l);
% Mwz=-Frx*(a*sin(pexi_r)*cos(pexi_r)+b*sin(phi_r))+Fry*a*cos(pexi_r)+Fly*a*cos(pexi_l)+Flx*(a*sin(pexi_l)*cos(phi_l)+b*sin(phi_l));
% Mwx=0;
% Mwy=0;
% Mwz=2*Fry*a*cos(pexi_r);

A_l=[0            -sin(phi_l)      cos(phi_l);
    -sin(phi_l)    0               0;
    -cos(phi_l)    0               0];
A_r=[0              sin(phi_r)     -cos(phi_r);
     sin(phi_r)     0               0;
     cos(phi_r)     0               0];

d1_l=0.5*rho*C_d*phi_l_dot*phi_l_dot+0.5*rho*C_rot*phi_l_dot*phi_l_dot;
d2_l=0.5*rho*C_l*phi_l_dot*phi_l_dot*sin(phi_l)+0.5*rho*C_rot*phi_l_dot*phi_l_dot*sin(phi_l);
d1_l=0.5*rho*C_d*phi_l_dot*phi_l_dot;
d2_l=0.5*rho*C_l*phi_l_dot*phi_l_dot*sin(phi_l);
d1_rot_l=0.5*rho*C_rot*phi_l_dot*phi_l_dot;
d2_rot_l=0.5*rho*C_rot*phi_l_dot*phi_l_dot*sin(phi_l);
if (remainder<pi)
    d2_l=d2_l;
else
    d2_l=-d2_l;
end
d3_l=0.5*rho*C_rot*phi_l_dot*phi_l_dot*cos(phi_l)+0.5*rho*C_l*phi_l_dot*phi_l_dot*cos(phi_l);
d3_l=0.5*rho*C_l*phi_l_dot*phi_l_dot*cos(phi_l);
d3_rot_l=0.5*rho*C_rot*phi_l_dot*phi_l_dot*cos(phi_l);
if (remainder<pi)
    d3_l=d3_l;
else
    d3_l=-d3_l;
end
d_l=[d1_l;d2_l;d3_l];
d_rot_l=[d1_rot_l;d2_rot_l;d3_rot_l];
d1_r=0.5*rho*C_d*phi_r_dot*phi_r_dot+0.5*rho*C_rot*phi_r_dot*phi_r_dot;
d2_r=-0.5*rho*C_l*phi_r_dot*phi_r_dot*sin(phi_r)-0.5*rho*C_rot*phi_r_dot*phi_r_dot*sin(phi_r);
d1_r=0.5*rho*C_d*phi_r_dot*phi_r_dot;
d2_r=-0.5*rho*C_l*phi_r_dot*phi_r_dot*sin(phi_r);
d1_rot_r=0.5*rho*C_rot*phi_r_dot*phi_r_dot;
d2_rot_r=-0.5*rho*C_rot*phi_r_dot*phi_r_dot*sin(phi_r);
if (remainder<pi)
    d2_r=d2_r;
else
    d2_r=-d2_r;
end
d3_r=-0.5*rho*C_rot*phi_r_dot*phi_r_dot*cos(phi_r)-0.5*rho*C_l*phi_r_dot*phi_r_dot*cos(phi_r);
d3_r=-0.5*rho*C_l*phi_r_dot*phi_r_dot*cos(phi_r);
d3_rot_r=-0.5*rho*C_rot*phi_r_dot*phi_r_dot*cos(phi_r);
if (remainder<pi)
    d3_r=d3_r;
else
    d3_r=-d3_r;
end
d_r=[d1_r;d2_r;d3_r];
d_rot_r=[d1_rot_r;d2_rot_r;d3_rot_r];
Mxyz=(A_l*L_l*d_l+A_r*L_r*d_r)*c_bar*0.25*R*R*R*R;
Mxyz_rot=(A_l*L_l*d_rot_l+A_r*L_r*d_rot_r)*c_bar*c_bar*c_bar*c_bar*R*R*R*R*R*R/6;
Mxyz=Mxyz+Mxyz_rot;
Mwx=Mxyz(1);
Mwy=Mxyz(2);
Mwz=Mxyz(3);


T=[0            sin(theta)                  1
   sin(gamma)   cos(theta)*cos(gamma)       0
   cos(gamma)  -cos(theta)*sin(gamma)       0];
   
p = phi_dot*sin(theta)+gamma_dot;
q = phi_dot*cos(theta)*cos(gamma)+theta_dot*sin(gamma);
r = -phi_dot*cos(theta)*sin(gamma)+theta_dot*cos(gamma);
P=[p;q;r];

f1=q*sin(gamma+gamma_bolang)+r*cos(gamma+gamma_bolang);
f2=(1/(cos(theta+theta_bolang)))*(q*cos(gamma+gamma_bolang)-r*sin(gamma+gamma_bolang));
f3=p-tan(theta+theta_bolang)*(q*cos(gamma+gamma_bolang)-r*sin(gamma+gamma_bolang));

f4=(Mbx+Mux+Mwx+Mdx-(Jz-Jy)*r*q)/Jx;
f5=(Mby+Muy+Mwy+Mdy-(Jx-Jz)*r*p)/Jy;
f6=(Mbz+Muz+Mwz+Mdz-(Jy-Jx)*p*q)/Jz;
P_dot=[f4;f5;f6];
T_dot=[0                       cos(theta)*theta_dot                                                 0
       cos(gamma)*gamma_dot   -cos(gamma)*sin(theta)*theta_dot-cos(theta)*sin(gamma)*gamma_dot      0
      -sin(gamma)*gamma_dot     sin(gamma)*sin(theta)*theta_dot-cos(theta)*cos(gamma)*gamma_dot     0];
f456=inv(T)*(P_dot-T_dot*status_x2);
f456=T_dot*P+T_dot*P_dot;
dx=[theta_dot;phi_dot;gamma_dot;f456];