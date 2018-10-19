clc
close all
clear all 
t=0:1/25/100:1/25;
 x=[0;0;0;0;0;0]; 
 tau=0.2;
 f=25;
 step1=0.2;
  step2=0.2;
   step3=0.2;
    step4=0.2;
 v1n=0;
for v1=-1:step1:1
    v1n=v1n+1;
     v2n=0;
    for v2=-1:step2:1
        v2n=v2n+1;
         v3n=0;
        for v3=-1:step3:1
            v3n=v3n+1;
             v4n=0;
            for v4=-1:step4:1   
                v4n=v4n+1;
 for k=1:length(t)
     [f_b(k,:),tau_b(k,:),F_D(k,:),F_L(k,:),pexi_l(k,:),pexi_r(k,:),phi_l(k,:),phi_r(k,:)]=attitude_deng(t(k),x,v1,v2,v3,v4,f,tau);
 end
%  px= f_b(:,1);
% pxm=mean(px)
% py= f_b(:,2);
% pym=mean(py)
% pz= f_b(:,3);
% pzm=mean(pz)
mx= tau_b(:,1);
mxm=mean(mx);
tau_x(v1n,v2n,v3n,v4n)=mxm;
my= tau_b(:,2);
mym=mean(my);
tau_y(v1n,v2n,v3n,v4n)=mym;
mz= tau_b(:,3);
mzm=mean(mz);
tau_z(v1n,v2n,v3n,v4n)=mzm;
            end
        end
    end
end
v1=[-1:step1:1];
v2=[-1:step2:1];
v3=[-1:step3:1];
v4=[-1:step4:1];
tau_xxx=[];
tau_yyy=[];
tau_zzz=[];
V1=[];
V2=[];
V3=[];
V4=[];
v1n=0;
for v1n=1:length(v1)
    for v2n=1:length(v2)
        for v3n=1:length(v3)
            for v4n=1:length(v4) 
                V1((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n)=v1(v1n);   
                V2((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n)=v2(v2n);   
                V3((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n)=v3(v3n);   
                V4((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n)=v4(v4n);   
%                 tau_z((v1n)*length(v1),v2n,v3n,v4n)=mzm;
                tau_xxx((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n)=tau_x(v1n,v2n,v3n,v4n);
                tau_yyy((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n)=tau_y(v1n,v2n,v3n,v4n);
                tau_zzz((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n)=tau_z(v1n,v2n,v3n,v4n);
            end
        end
    end
end

x=[V1;V2;V3;V4];
t=[tau_xxx];
% [x,t] = simplefit_dataset;
tmin=min(t);
tmax=max(t);
t=(t-min(t))/(max(t)-min(t));
     net = feedforwardnet(10);
%      net=configure(net,x,t);
%      net=configure(net,'outputs',t,1);
%      net=configure(net,'inputs',x,2);
     net = train(net,x,t);
%      view(net)
     y = net(x);
          perf = perform(net,t,y')
for v1n=1:length(v1)
    for v2n=1:length(v2)
        for v3n=1:length(v3)
            for v4n=1:length(v4) 
                tau_xxxx(v1n,v2n,v3n,v4n) = y((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n);
%                 tau_yyy((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n)=tau_y(v1n,v2n,v3n,v4n);
%                 tau_zzz((v1n-1)*length(v2)*length(v3)*length(v4)+(v2n-1)*length(v3)*length(v4)+(v3n-1)*length(v4)+v4n)=tau_z(v1n,v2n,v3n,v4n);
            end
        end
    end
end

     plot(V1,tau_xxx,'r')
     figure
     plot(V2,tau_xxx,'b')
     figure
     plot(V1,y(1,:))
    tau_xxxxx =tau_xxxx(:,:,1,1)*(tmax-tmin)+tmin;
     figure
mesh(v1,v2,tau_x(:,:,1,1))
% figure
% plot(t,f_b)
% xlabel('Time(s)')
% ylabel('Angle(rad)')
% legend('1','2','3')
% hold on; 
% grid on;
figure
mesh(v1,v2,tau_xxxxx(:,:,1,1))
subtau_xxxxx=tau_x(:,:,1,1)-tau_xxxxx(:,:,1,1);
figure
mesh(v1,v2,subtau_xxxxx(:,:,1,1))

figure
plot(t,tau_b)
xlabel('Time(s)')
ylabel('Angle(rad)')
legend('1','2','3')
hold on; 
grid on;
% figure(1)
% plot(t,f_b(:,1))
% xlabel('Time(s)')
% ylabel('Angle(rad)')
% legend('1','2','3')
% hold on; 
% grid on;
% figure(2)
% plot(t,f_b(:,2))
% xlabel('Time(s)')
% ylabel('Angle(rad)')
% legend('1','2','3')
% hold on; 
% grid on;
% figure(3)
% plot(t,f_b(:,3))
% xlabel('Time(s)')
% ylabel('Angle(rad)')
% legend('1','2','3')
% hold on; 
% grid on;
% 
% figure(4)
% plot(t,tau_b(:,1))
% xlabel('Time(s)')
% ylabel('Angle(rad)')
% legend('1','2','3')
% hold on; 
% grid on;
% figure(5)
% plot(t,tau_b(:,2))
% xlabel('Time(s)')
% ylabel('Angle(rad)')
% legend('1','2','3')
% hold on; 
% grid on;
% figure(6)
% plot(t,tau_b(:,3))
% xlabel('Time(s)')
% ylabel('Angle(rad)')
% legend('1','2','3')
% hold on; 
% grid on;

% figure
% plot(t,F_D)
% xlabel('Time(s)')
% ylabel('Force(N)') 
% legend('DRAG')
% figure
% plot(t,F_L)
% xlabel('Time(s)')
% ylabel('Force(N)')
% legend('LIFT')
% figure
% plot(t,pexi_l)
% xlabel('Time(s)')
% ylabel('Force(N)')
% legend('Psi_l')
% figure
% plot(t,pexi_r)
% xlabel('Time(s)')
% ylabel('Force(N)')
% legend('Psi_r')
% figure
% plot(t,phi_l)
% xlabel('Time(s)')
% ylabel('Force(N)')
% legend('Phi_l')
% figure
% plot(t,phi_r)
% xlabel('Time(s)')
% ylabel('Force(N)')
% legend('Phi_r')
