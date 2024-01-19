function [dh_dx, dg_dx, H_,L_Z, U_W]=coe(Y,bus,gen,branch,X,A2)
num_node=size(bus,1);
num_gen=size(gen,1);  
num_branch=size(branch,1);
num_equa=2*num_node;
num_inequa=2*num_gen+num_node+num_branch;
state=2*(num_gen+num_node); 

z = X(1:num_inequa); 
l = X(num_inequa+1:2*num_inequa); 
w = X(2*num_inequa+1:3*num_inequa);
u = X(3*num_inequa+1:4*num_inequa); 
x = X(4*num_inequa+1:5*num_inequa); 
y = X(5*num_inequa+1:5*num_inequa+num_equa);
Xtilde=x(2*num_gen+1:state);
%% 计算等式约束Jacobian矩阵
dh_dPG=zeros(num_gen,2*num_node);
dh_dQR=zeros(num_gen,2*num_node);
for i=1:num_gen
   dh_dPG(i,gen(i,1)*2-1)=1;
   dh_dQR(i,gen(i,1)*2)=1;
end
dh_dXtilde=zeros(2*num_node,2*num_node);  %潮流计算中的雅可比矩阵，详见《电力系统分析》
for i=1:num_node
    for j=1:num_node       
        if(i~=j)
            theta = Xtilde(i*2-1)-Xtilde(j*2-1);
            dh_dXtilde(j*2-1,i*2-1)=-Xtilde(i*2)*Xtilde(j*2)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta));     %Hij
            dh_dXtilde(j*2-1,i*2)=Xtilde(i*2)*Xtilde(j*2)*(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta));        %Jij
            dh_dXtilde(j*2,i*2-1)=-Xtilde(i*2)*(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta));                  %Nij
            dh_dXtilde(j*2,i*2)=-Xtilde(i*2)*Xtilde(j*2)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta));       %Lij
        
            dh_dXtilde(i*2-1,i*2-1)=dh_dXtilde(i*2-1,i*2-1)+Xtilde(i*2)*Xtilde(j*2)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta)); %Hii
            dh_dXtilde(i*2-1,i*2)=dh_dXtilde(i*2-1,i*2)-Xtilde(i*2)*Xtilde(j*2)*(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta));       %Jii
            dh_dXtilde(i*2,i*2-1)=dh_dXtilde(i*2,i*2-1)-Xtilde(j*2)*(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta));       %Nii
            dh_dXtilde(i*2,i*2)=dh_dXtilde(i*2,i*2)-Xtilde(i*2)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta));         %Lii
        end
    end
    dh_dXtilde(i*2,i*2-1)=dh_dXtilde(i*2,i*2-1)-2*Xtilde(i*2)*real(Y(i,i));     %附加项Nii
    dh_dXtilde(i*2,i*2)=dh_dXtilde(i*2,i*2)+2*Xtilde(i*2)*Xtilde(i*2)*imag(Y(i,i))/Xtilde(j*2);         %附加项Lii
end
dh_dx = [dh_dPG; dh_dQR; dh_dXtilde];
%% 计算不等式约束Jacobian矩阵
dg1_dPG = eye(num_gen); 
dg1_dQR = zeros(num_gen,num_gen); 
dg1_dXtilde = zeros(2*num_node,num_gen);
dg2_dPG = zeros(num_gen,num_gen); 
dg2_dQR = eye(num_gen); 
dg2_dXtilde = zeros(2*num_node,num_gen);
dg3_dPG = zeros(num_gen,num_node); 
dg3_dQR = zeros(num_gen,num_node);
dg3_dXtilde = zeros(2*num_node,num_node); 
dg4_dPG = zeros(num_gen,num_branch); 
dg4_dQR = zeros(num_gen,num_branch);
dg4_dXtilde = zeros(2*num_node,num_branch);

for i = 1:num_node
    dg3_dXtilde(2*i,i) = 1;
end

for k=1:num_branch
   i=branch(k,1); 
   j=branch(k,2);
   theta=Xtilde(i*2-1)-Xtilde(j*2-1);
   dg4_dXtilde(2*i-1,k)=-Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta));
   dg4_dXtilde(2*j-1,k)=Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta));
   dg4_dXtilde(2*i,k)=Xtilde(2*j)*(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta))-2*Xtilde(2*i)*real(Y(i,j));
   dg4_dXtilde(2*j,k)=Xtilde(2*i)*(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta));
end
dg_dx = [dg1_dPG       dg2_dPG       dg3_dPG        dg4_dPG;
         dg1_dQR       dg2_dQR       dg3_dQR        dg4_dQR;
         dg1_dXtilde   dg2_dXtilde   dg3_dXtilde    dg4_dXtilde];
%% 计算对角矩阵
L_Z = diag(z./l); 
U_W = diag(w./u);
%% 计算目标函数的Hessian矩阵
d2f_dx=zeros(state,state);
d2f_dx(1:num_gen,1:num_gen)=2*A2;
%% 计算等式约束的Hessian矩阵算子与Lagrange乘子y乘积
d2h_dx_y=zeros(state,state);
a=zeros(num_equa,num_equa);
for i=1:num_node
    for j=1:num_node
    theta=Xtilde(i*2-1)-Xtilde(j*2-1);  
        if(j~=i)
a(2*i-1,2*i-1)=a(2*i-1,2*i-1)+Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*(-cos(theta)*y(2*j-1)+sin(theta)*y(2*j))+imag(Y(i,j))*(-sin(theta)*y(2*j-1)-cos(theta)*y(2*j)));
a(2*i-1,2*i)=a(2*i-1,2*i)+Xtilde(2*j)*(real(Y(i,j))*(-sin(theta)*y(2*j-1)-cos(theta)*y(2*j))+imag(Y(i,j))*(cos(theta)*y(2*j-1)-sin(theta)*y(2*j)));
a(2*i,2*i-1)=a(2*i,2*i-1)+Xtilde(2*j)*(real(Y(i,j))*(-sin(theta)*y(2*j-1)-cos(theta)*y(2*j))+imag(Y(i,j))*(cos(theta)*y(2*j-1)-sin(theta)*y(2*j)));  
a(2*i,2*i)=-2*(real(Y(i,i))*y(2*j-1)-2*imag(Y(i,i))*y(2*j));            
a(2*i-1,2*j-1)=a(2*i-1,2*j-1)+Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*(cos(theta)*y(2*j-1)-sin(theta)*y(2*j))+imag(Y(i,j))*(sin(theta)*y(2*j-1)+cos(theta)*y(2*j)));
a(2*i-1,2*j)=a(2*i-1,2*j)+Xtilde(2*i)*(real(Y(i,j))*(-sin(theta)*y(2*j-1)-cos(theta)*y(2*j))+imag(Y(i,j))*(cos(theta)*y(2*j-1)-sin(theta)*y(2*j)));        
a(2*i,2*j-1)=a(2*i,2*j-1)+Xtilde(2*j)*(real(Y(i,j))*(sin(theta)*y(2*j-1)+cos(theta)*y(2*j))+imag(Y(i,j))*(-cos(theta)*y(2*j-1)+sin(theta)*y(2*j)));
a(2*i,2*j)=a(2*i,2*j)+(real(Y(i,j))*(cos(theta)*y(2*j-1)-sin(theta)*y(2*j))+imag(Y(i,j))*(sin(theta)*y(2*j-1)+cos(theta)*y(2*j)));
        end
    end
end
d2h_dx_y(2*num_gen+1:state,2*num_gen+1:state)=a;
%% 计算不等式约束的Hessian矩阵算子与Lagrange乘子z+w乘积
d2g_dx_c=zeros(state,state);
d2g4_d2Xtilde=zeros(length(Xtilde),length(Xtilde));    
c=z+w;
for k=1:num_branch                                 %细细品
   i=branch(k,1); 
   j=branch(k,2);
   theta=Xtilde(i*2-1)-Xtilde(j*2-1);

   d2g4_d2Xtilde(2*i-1,2*i-1) = d2g4_d2Xtilde(2*i-1,2*i-1)-Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta))*c(2*num_gen+num_node+k);
   d2g4_d2Xtilde(2*i-1,2*i) = d2g4_d2Xtilde(2*i-1,2*i)+(-Xtilde(2*j)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta)))*c(2*num_gen+num_node+k);
   d2g4_d2Xtilde(2*i-1,2*j-1) = d2g4_d2Xtilde(2*i-1,2*j-1)+(Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta)))*c(2*num_gen+num_node+k);
   d2g4_d2Xtilde(2*i-1,2*j) = d2g4_d2Xtilde(2*i-1,2*j)-Xtilde(2*i)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta))*c(2*num_gen+num_node+k);
   
   d2g4_d2Xtilde(2*j-1,2*i-1) = d2g4_d2Xtilde(2*i-1,2*j-1);
   d2g4_d2Xtilde(2*j-1,2*i) = d2g4_d2Xtilde(2*j-1,2*i)+(Xtilde(2*j)*(real(Y(i,j))*sin(theta)-imag(Y(i,j))*cos(theta)))*c(2*num_gen+num_node+k);
   d2g4_d2Xtilde(2*j-1,2*j-1) =  d2g4_d2Xtilde(2*i-1,2*i-1);
   d2g4_d2Xtilde(2*j-1,2*j) = d2g4_d2Xtilde(2*j-1,2*j)+(Xtilde(2*i)*(real(Y(i,j))*sin(theta))-imag(Y(i,j))*cos(theta))*c(2*num_gen+num_node+k);
   
   d2g4_d2Xtilde(2*i,2*i-1) = d2g4_d2Xtilde(2*i-1,2*i);
   d2g4_d2Xtilde(2*i,2*i) = d2g4_d2Xtilde(2*i,2*i)+2*real(Y(i,j))*c(2*num_gen+num_node+k);
   d2g4_d2Xtilde(2*i,2*j-1) = d2g4_d2Xtilde(2*j-1,2*i);
   d2g4_d2Xtilde(2*i,2*j) = d2g4_d2Xtilde(2*i,2*j)+(real(Y(i,j))*cos(theta)+imag(Y(i,j))*sin(theta))*c(2*num_gen+num_node+k);
   
   d2g4_d2Xtilde(2*j,2*i-1) = d2g4_d2Xtilde(2*i-1,2*j);
   d2g4_d2Xtilde(2*j,2*i) = d2g4_d2Xtilde(2*i,2*j);
   d2g4_d2Xtilde(2*j,2*j-1) = d2g4_d2Xtilde(2*j-1,2*j);
end
d2g_dx_c(2*num_gen+1:state,2*num_gen+1:state)=d2g4_d2Xtilde;  %右下角补缺
mix=dg_dx*(L_Z-U_W)*dg_dx';
H_=-d2f_dx+d2h_dx_y+d2g_dx_c-mix;           %海塞矩阵
end