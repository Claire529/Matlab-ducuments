
% loop for fuzzy system
function s_fun=Quaternion_calCtrl_fuz_adap(e,ce,e_center,ce_center,e_width,ce_width)
maxNum=50.0;
e_n=length(e_center);
ce_n=length(ce_center);
for j=1:ce_n
    if (j==1)
        for i=1:e_n
              if (i==1)
                  s_fun((i+e_n*(j-1)))=trapezoid(e,-maxNum,-maxNum,e_center(i),e_center(i+1))*trapezoid(ce,-maxNum,-maxNum,ce_center(j),ce_center(j+1));
              elseif (i==e_n)
                  s_fun((i+e_n*(j-1)))=trapezoid(e,e_center(i-1),e_center(i),maxNum,maxNum)*trapezoid(ce,-maxNum,-maxNum,ce_center(j),ce_center(j+1));
              else
                  s_fun((i+e_n*(j-1)))=trapezoid(e,e_center(i-1),e_center(i),e_center(i),e_center(i+1))*trapezoid(ce,-maxNum,-maxNum,ce_center(j),ce_center(j+1));
              end
         end
    elseif (j==ce_n)
        for i=1:e_n
              if (i==1)
                  s_fun((i+e_n*(j-1)))=trapezoid(e,-maxNum,-maxNum,e_center(i),e_center(i+1))*trapezoid(ce,ce_center(j-1),ce_center(j),maxNum,maxNum);
              elseif (i==e_n)
                  s_fun((i+e_n*(j-1)))=trapezoid(e,e_center(i-1),e_center(i),maxNum,maxNum)*trapezoid(ce,ce_center(j-1),ce_center(j),maxNum,maxNum);
              else
                  s_fun((i+e_n*(j-1)))=trapezoid(e,e_center(i-1),e_center(i),e_center(i),e_center(i+1))*trapezoid(ce,ce_center(j-1),ce_center(j),maxNum,maxNum);
              end
         end
       
    else
        for i=1:e_n
            if (i==1)
                s_fun((i+e_n*(j-1)))=trapezoid(e,-maxNum,-maxNum,e_center(i),e_center(i+1))*trapezoid(ce,ce_center(j-1),ce_center(j),ce_center(j),ce_center(j+1));
            elseif (i==e_n)
                s_fun((i+e_n*(j-1)))=trapezoid(e,e_center(i-1),e_center(i),maxNum,maxNum)*trapezoid(ce,ce_center(j-1),ce_center(j),ce_center(j),ce_center(j+1));
            else
               s_fun((i+e_n*(j-1)))=trapezoid(e,e_center(i-1),e_center(i),e_center(i),e_center(i+1))*trapezoid(ce,ce_center(j-1),ce_center(j),ce_center(j),ce_center(j+1));
            end
        end
    end
end
sum = 0;
for k=1:e_n*ce_n
    sum = sum+s_fun(k);
end
for k=1:e_n*ce_n
   s_fun(k) = s_fun(k)/sum;
end


% two rules

% function u=more_calCtrl_fuzzy_adaptive_new_integral(e,ce,sigma1,sigma2,c1,c2)
% 
% 
% uf11 = 1/(1+exp(((e-c1)/sigma1)));
% uf21 = 1/(1+exp(((ce-c2)/sigma2)));
% uf12 = 1/(1+exp(-((e-c1)/sigma1)));
% uf22 = 1/(1+exp(-((ce-c2)/sigma2)));
% 
% sum = uf11*uf21+uf11*uf22+uf12*uf21+uf12*uf22;
% s11 = uf11*uf21/sum;
% s12 = uf11*uf22/sum;
% s13 = uf12*uf21/sum;
% s14 = uf12*uf22/sum;
% s1 = [s11;s12;s13;s14];
% u = s1;
