load tau_x;
load tau_y;
load tau_z;
load f_x;
load f_y;
load f_z;
v1=[-1:step1:1];
v2=[-1:step2:1];
v3=[-1:step3:1];
v4=[-1:step4:1];
v5=[-1:step5:1];
v6=[-1:step6:1];
tau_xxx=[];
tau_yyy=[];
tau_zzz=[];
V1=[];
V2=[];
V3=[];
V4=[];
V5=[];
V6=[];
for v1n=1:length(v1)
    for v2n=1:length(v2)
        for v3n=1:length(v3)
            for v4n=1:length(v4) 
                for v5n=1:length(v5)
                    for v6n=1:length(v6) 
                        nnn=(v1n-1)*length(v2)*length(v3)*length(v4)*length(v5)*length(v6)+(v2n-1)*length(v3)*length(v4)*length(v5)*length(v6)+(v3n-1)*length(v4)*length(v5)*length(v6)+(v4n-1)*length(v5)*length(v6)+(v5n-1)*length(v6)+v6n;
                        V1(nnn)=v1(v1n);   
                        V2(nnn)=v2(v2n);   
                        V3(nnn)=v3(v3n);   
                        V4(nnn)=v4(v4n);   
                        V5(nnn)=v5(v5n);   
                        V6(nnn)=v6(v6n);   
                        f_xxx(nnn)=f_x(v1n,v2n,v3n,v4n,v5n,v6n);
                        f_yyy(nnn)=f_y(v1n,v2n,v3n,v4n,v5n,v6n);
                        f_zzz(nnn)=f_y(v1n,v2n,v3n,v4n,v5n,v6n);
                        tau_xxx(nnn)=tau_x(v1n,v2n,v3n,v4n,v5n,v6n);
                        tau_yyy(nnn)=tau_y(v1n,v2n,v3n,v4n,v5n,v6n);
                        tau_zzz(nnn)=tau_z(v1n,v2n,v3n,v4n,v5n,v6n);
                    end
                end
            end
        end
    end
end

x=[V1;V2;V3;V4;V5;V6];
t=[f_xxx];
tmin=min(t);
tmax=max(t);
t=(t-min(t))/(max(t)-min(t));
net = feedforwardnet(10);
net.trainParam.epochs =200; % Max number of iterations
net.trainParam.goal = 1e-5; % Error tolerance; stopping criterion
net = train(net,x,t);
z_x = net(x);
save z_x

t=[f_yyy];
i=1
tmin=min(t);
tmax=max(t);
t=(t-min(t))/(max(t)-min(t));
net = feedforwardnet(10);
net.trainParam.epochs =200; % Max number of iterations
net.trainParam.goal = 1e-5; % Error tolerance; stopping criterion
net = train(net,x,t);
z_y = net(x);
save z_y

t=[f_zzz];
i=2
% [x,t] = simplefit_dataset;
tmin=min(t);
tmax=max(t);
t=(t-min(t))/(max(t)-min(t));
net = feedforwardnet(10);
net.trainParam.epochs =200; % Max number of iterations
net.trainParam.goal = 1e-5; % Error tolerance; stopping criterion
net = train(net,x,t);
z_z = net(x);
save z_z

i=3
t=[tau_xxx];
tmin=min(t);
tmax=max(t);
t=(t-min(t))/(max(t)-min(t));
net = feedforwardnet(10);
net.trainParam.epochs =200; % Max number of iterations
net.trainParam.goal = 1e-5; % Error tolerance; stopping criterion
net = train(net,x,t);
y_x = net(x);
save y_x

i=4
t=[tau_yyy];
tmin=min(t);
tmax=max(t);
t=(t-min(t))/(max(t)-min(t));
net = feedforwardnet(10);
net.trainParam.epochs =200; % Max number of iterations
net.trainParam.goal = 1e-5; % Error tolerance; stopping criterion
net = train(net,x,t);
y_y = net(x);
save y_y

i=5
t=[tau_zzz];
tmin=min(t);
tmax=max(t);
t=(t-min(t))/(max(t)-min(t));
net = feedforwardnet(10);
net.trainParam.epochs =200; % Max number of iterations
net.trainParam.goal = 1e-5; % Error tolerance; stopping criterion
net = train(net,x,t);
y_z = net(x);
save y_z

i=6
for v1n=1:length(v1)
    for v2n=1:length(v2)
        for v3n=1:length(v3)
            for v4n=1:length(v4) 
                for v5n=1:length(v5)
                    for v6n=1:length(v6) 
                        nnn=(v1n-1)*length(v2)*length(v3)*length(v4)*length(v5)*length(v6)+(v2n-1)*length(v3)*length(v4)*length(v5)*length(v6)+(v3n-1)*length(v4)*length(v5)*length(v6)+(v4n-1)*length(v5)*length(v6)+(v5n-1)*length(v6)+v6n;
                        f_xxxx(v1n,v2n,v3n,v4n,v5n,v6n) = z_x(nnn);
                        f_yyyy(v1n,v2n,v3n,v4n,v5n,v6n) = z_y(nnn);
                        f_zzzz(v1n,v2n,v3n,v4n,v5n,v6n) = z_z(nnn);
                        tau_xxxx(v1n,v2n,v3n,v4n,v5n,v6n) = y_x(nnn);
                        tau_yyyy(v1n,v2n,v3n,v4n,v5n,v6n) = y_y(nnn);
                        tau_zzzz(v1n,v2n,v3n,v4n,v5n,v6n) = y_z(nnn);
                    end
                end
            end
        end
    end
end
save f_xxxx;
save f_yyyy;
save f_zzzz;
save tau_xxxx;
save tau_yyyy;
save tau_zzzz;

load tau_xxxx;
load tau_yyyy;
load tau_zzzz;

load tau_x;
load tau_y;
load tau_z;
v1=[-1:step1:1];
v2=[-1:step2:1];
v3=[-1:step3:1];
v4=[-1:step4:1];
v5=[-1:step5:1];
v6=[-1:step6:1];

f_xxxx_min=min(min(min(min(min(min(f_x))))));
f_xxxx_max=max(max(max(max(max(max(f_x))))));

f_xxxxx =f_xxxx(:,:,1,1)*(f_xxxx_max-f_xxxx_min)+f_xxxx_min;

f_yyyy_min=min(min(min(min(min(min(f_y))))));
f_yyyy_max=max(max(max(max(max(max(f_y))))));

f_yyyyy =f_yyyy(:,:,1,1)*(f_yyyy_max-f_yyyy_min)+f_yyyy_min;

f_zzzz_min=min(min(min(min(min(min(f_y))))));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% y-->z
f_zzzz_max=max(max(max(max(max(max(f_y))))));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% y-->z

f_zzzzz =f_zzzz(:,:,1,1)*(f_zzzz_max-f_zzzz_min)+f_zzzz_min;

tau_xxxx_min=min(min(min(min(min(min(tau_x))))));
tau_xxxx_max=max(max(max(max(max(max(tau_x))))));

tau_xxxxx =tau_xxxx(:,:,1,1)*(tau_xxxx_max-tau_xxxx_min)+tau_xxxx_min;

tau_yyyy_min=min(min(min(min(min(min(tau_y))))));
tau_yyyy_max=max(max(max(max(max(max(tau_y))))));

tau_yyyyy =tau_yyyy(:,:,1,1)*(tau_yyyy_max-tau_yyyy_min)+tau_yyyy_min;

tau_zzzz_min=min(min(min(min(min(min(tau_z))))));
tau_zzzz_max=max(max(max(max(max(max(tau_z))))));

tau_zzzzz =tau_zzzz(:,:,1,1)*(tau_zzzz_max-tau_zzzz_min)+tau_zzzz_min;

figure
mesh(v1,v2,f_x(:,:,1,1))
figure
mesh(v1,v2,f_xxxxx(:,:,1,1))
subtau_xxxxx=f_x(:,:,1,1)-f_xxxxx(:,:,1,1);
figure
mesh(v1,v2,subtau_xxxxx(:,:,1,1))

figure
mesh(v1,v2,f_y(:,:,1,1))
figure
mesh(v1,v2,f_yyyyy(:,:,1,1))
subtau_yyyyy=f_y(:,:,1,1)-f_yyyyy(:,:,1,1);
figure
mesh(v1,v2,subtau_yyyyy(:,:,1,1))


figure
mesh(v1,v2,f_y(:,:,1,1))%%%%%%%%%%%%%%%%%%%%% y-->z
figure
mesh(v1,v2,f_zzzzz(:,:,1,1))
subtau_zzzzz=f_y(:,:,1,1)-f_zzzzz(:,:,1,1);%%%%%%%%%%%%%%%%%%%%% y-->z
figure
mesh(v1,v2,subtau_zzzzz(:,:,1,1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
mesh(v1,v2,tau_x(:,:,1,1))
figure
mesh(v1,v2,tau_xxxxx(:,:,1,1))
subtau_xxxxx=tau_x(:,:,1,1)-tau_xxxxx(:,:,1,1);
figure
mesh(v1,v2,subtau_xxxxx(:,:,1,1))

figure
mesh(v1,v2,tau_y(:,:,1,1))
figure
mesh(v1,v2,tau_yyyyy(:,:,1,1))
subtau_yyyyy=tau_y(:,:,1,1)-tau_yyyyy(:,:,1,1);
figure
mesh(v1,v2,subtau_yyyyy(:,:,1,1))


figure
mesh(v1,v2,tau_z(:,:,1,1))
figure
mesh(v1,v2,tau_zzzzz(:,:,1,1))
subtau_zzzzz=tau_z(:,:,1,1)-tau_zzzzz(:,:,1,1);
figure
mesh(v1,v2,subtau_zzzzz(:,:,1,1))