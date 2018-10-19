% the relationship between force and parameters
load tau_x;
load tau_y;
load tau_z;
load f_x;
load f_y;
load f_z;
tau_x_min=min(min(min(min(min(min(tau_x))))));
tau_x_max=max(max(max(max(max(max(tau_x))))));
tau_y_min=min(min(min(min(min(min(tau_y))))));
tau_y_max=max(max(max(max(max(max(tau_y))))));
tau_z_min=min(min(min(min(min(min(tau_z))))));
tau_z_max=max(max(max(max(max(max(tau_z))))));

f_x_min=min(min(min(min(min(min(f_x))))));
f_x_max=max(max(max(max(max(max(f_x))))));
f_y_min=min(min(min(min(min(min(f_y))))));
f_y_max=max(max(max(max(max(max(f_y))))));
f_z_min=min(min(min(min(min(min(f_y))))));%%%%  y--.z
f_z_max=max(max(max(max(max(max(f_y))))));%%%%  y--.z

NOP=5;
step_tau_x=(tau_x_max-tau_x_min)/NOP;
step_tau_y=(tau_y_max-tau_y_min)/NOP;
step_tau_z=(tau_z_max-tau_z_min)/NOP;
step_f_x=(f_x_max-f_x_min)/NOP;
step_f_y=(f_y_max-f_y_min)/NOP;
step_f_z=(f_z_max-f_z_min)/NOP;
tau_xx=[tau_x_min:step_tau_x:tau_x_max];
tau_yy=[tau_y_min:step_tau_y:tau_y_max];
tau_zz=[tau_z_min:step_tau_z:tau_z_max];
f_xx=[f_x_min:step_f_x:f_x_max];
f_yy=[f_y_min:step_f_y:f_y_max];
f_zz=[f_z_min:step_f_z:f_z_max];
v1=tau_xx;
v2=tau_yy;
v3=tau_zz;
v1=f_xx;
v2=f_yy;
v3=f_zz;
while(1)
                        f=[v1(1);v2(1);v3(1);v4(1);v5(1);v6(1)];
                        [V,FVAL,EXITFLAG]=fsolve(@(V) FW_MAV_Deng_V2f(V,f),[1 1 1 1 1 1]*rand());
                        if EXITFLAG==1
                            V1(v1n,v2n,v3n,v4n,v5n,v6n)=V(1);
                            V2(v1n,v2n,v3n,v4n,v5n,v6n)=V(2);
                            V3(v1n,v2n,v3n,v4n,v5n,v6n)=V(3);
                            V4(v1n,v2n,v3n,v4n,v5n,v6n)=V(4);
                            V5(v1n,v2n,v3n,v4n,v5n,v6n)=V(5);
                            V6(v1n,v2n,v3n,v4n,v5n,v6n)=V(6);
                        else
                            f=[v1(1);v2(1);v3(1);v4(1);v5(1);v6(1)];
                        end
                        
end

for v1n=1:length(v1)
    for v2n=1:length(v2)
        for v3n=1:length(v3)
            for v4n=1:length(v4) 
                for v5n=1:length(v5)
                    for v6n=1:length(v6) 
                        f=[v1(v1n);v2(v2n);v3(v3n);v4(v4n);v5(v5n);v6(v6n)];
                        [V,FVAL,EXITFLAG]=fsolve(@(V) FW_MAV_Deng_V2f(V,f),[1 1 1 1 1 1]*rand());
                        if EXITFLAG==1
                            V1(v1n,v2n,v3n,v4n,v5n,v6n)=V(1);
                            V2(v1n,v2n,v3n,v4n,v5n,v6n)=V(2);
                            V3(v1n,v2n,v3n,v4n,v5n,v6n)=V(3);
                            V4(v1n,v2n,v3n,v4n,v5n,v6n)=V(4);
                            V5(v1n,v2n,v3n,v4n,v5n,v6n)=V(5);
                            V6(v1n,v2n,v3n,v4n,v5n,v6n)=V(6);
                        end
                    end
                end
            end
        end
    end
end
