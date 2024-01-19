clc;clear all;warning off;
opfcase;
%% 基本参数设置
num_node=size(bus,1);                              %节点个数
num_branch=size(branch,1);                         %支路数
num_gen=size(gen,1);                               %发电机个数
bus_3=find(bus(:,2)==3);                           %平衡节点
bus_2=find(bus(:,2)==2);                           %PV节点
bus_1=find(bus(:,2)==1);                           %PQ节点
num_PQ=size(bus_1,1);
num_PV=size(bus_2,1);

num_equa = 2*num_node;                                 %等式约束个数
num_inequa = 2*num_gen+num_node+num_branch;            %不等式约束个数

a2=gencost(:,5);  a1=gencost(:,6);  a0=gencost(:,7);    %耗量特性多项式系数
A2=diag(a2);  A1=diag(a1);

g_u = [gen(:,9); gen(:,4); bus(:,12); branch(:,6)];
%发电机有功上界 发电机无功上界 节点电压上限 支路上限功率
g_l = [gen(:,10); gen(:,5); bus(:,13); -branch(:,6)];
%发电机有功下界 发电机无功下界 节点电压下限 支路下限功率

Xtilde(1:2:num_node*2-1) = (bus(:,9))';               %节点初始相角
Xtilde(2:2:num_node*2) = (bus(:,8))';                 %节点初始幅值

u=1*ones(num_inequa,1);                               %拉格朗日乘子初始赋值
l=1*ones(num_inequa,1);
z=1*ones(1,num_inequa); 
w=-0.5*ones(1,num_inequa); 
y=1E-10*ones(1,num_equa); 
y(2:2:num_equa)=-1*y(2:2:num_equa);            

epsi=1E-6;
sigma=0.1;
max_iteration=50; 
gap_record=zeros(max_iteration,1);
%% 求节点导纳矩阵
Y=solveY(branch,bus,num_node,num_branch);
%% 开始计算
for num_iteration = 1:max_iteration  
gap=l'*z'-u'*w';
gap_record(num_iteration) = gap;
    if (gap<=epsi)
        disp('iterations completed');
        break; 
    end

Pg=(gen(:,9)+gen(:,10))/2;                                     % 发电机有功出力
Qr=(gen(:,4)+gen(:,5))/2;                                      % 发电机无功出力
g1=Pg;                                                         % 发电机有功出力g1
g2=Qr;                                                         % 发电机无功出力g2
g3=(Xtilde(2:2:num_node*2))';                                  % 所有节点电压幅值g3   
    for k=1:num_branch                                         % 线路潮流g4
        i=branch(k,1);
        j=branch(k,2);
        cita=Xtilde(i*2-1)-Xtilde(j*2-1);
     g4(k,1) = Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*cos(cita)+imag(Y(i,j))*sin(cita))-...
               Xtilde(2*i)*Xtilde(2*i)*real(Y(i,j));     
    end
g=[g1; g2; g3; g4];                                         %  num_inequa*1
    
for i=1:num_node                                              %各节点净注入
P_ejec(i,1)=-bus(i,3); 
Q_ejec(i,1)=-bus(i,4);
  if bus(i,2)~=1
P_ejec(i,1)=P_ejec(i,1)+gen(gen(:,1)==i,2); 
Q_ejec(i,1)=Q_ejec(i,1)+gen(gen(:,1)==i,3);
  end
end  

h=zeros(num_node*2,1);                                        %对y求偏导，等式约束，num_equa*1
for i=1:num_node
     for j=1:num_node
            cita=Xtilde(i*2-1)-Xtilde(j*2-1);
            h(2*i-1)=h(2*i-1)-Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*cos(cita)+imag(Y(i,j))*sin(cita));
            h(2*i)=h(2*i)-Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*sin(cita)-imag(Y(i,j))*cos(cita));
     end
            h(2*i-1)=h(2*i-1)+P_ejec(i);
            h(2*i)=h(2*i)+Q_ejec(i);
end
    %% 求系数矩阵
    x=[Pg;Qr;Xtilde'];                                             % 状态变量集
    X=[z';l;w';u;x;y'];                                       % 迭代集
    [dh_dx,dg_dx,H_,L_Z,U_W]=coe(Y,bus,gen,branch,X,A2);
    
    Ly=h;                                                          % Lagrange函数偏导
    Lz=g-l-g_l;
    Lw=g+u-g_u;
    L=diag(l); U=diag(u); Z=diag(z); W=diag(w);
    miu=gap*sigma/2/num_inequa;
    L_miu_l=L*Z*ones(num_inequa,1)-miu*ones(num_inequa,1);
    L_miu_u=U*W*ones(num_inequa,1)+miu*ones(num_inequa,1);
    df_dx=[2*A2*Pg+a1;zeros(num_gen,1);zeros(length(Xtilde),1)];               %目标函数梯度矢量num_inequa*1
    Lx_=df_dx-dh_dx*y'-dg_dx*(z+w)'+dg_dx*(L\(L_miu_l+Z*Lz)+U\(L_miu_u-W*Lw));    
    b=[-L\L_miu_l;Lz;-U\L_miu_u;-Lw;Lx_;-Ly];   
    %% 求解矩阵方程
    delta_X=delta(H_,dg_dx,dh_dx,L_Z,U_W,b,z,l',u',w,x',y); 
    X=X+delta_X;    
    % 更新迭代量
    z=(X(1:num_inequa))';
    l=X(num_inequa+1:2*num_inequa);
    w=(X(2*num_inequa+1:3*num_inequa))';
    u=X(3*num_inequa+1:4*num_inequa);
    x=X(4*num_inequa+1:5*num_inequa);
    y=(X(5*num_inequa+1:5*num_inequa+num_equa))';
    % 更新状态量
    Pg=(x(1:num_gen));
    Qr=(x(num_gen+1:2*num_gen));
    Xtilde=(x(2*num_gen+1:2*num_gen+2*num_node))';   
end
if (num_iteration>=max_iteration)
  disp('OPF is not convergent');  
end
%% 结果分析
semilogy(gap_record(1:num_iteration));
grid on;
title('收敛特性');
xlabel('迭代次数'); 
ylabel('Gap');

%相角变换
for k=1:num_node
    Xtilde(2*k-1)=Xtilde(2*k-1)-Xtilde(2*num_node-1);
end

source=[gen(:,1) Pg Qr A2*(Pg.*Pg)+A1*Pg+a0];
voltage=[bus(:,1)  (Xtilde(2:2:2*num_node))'  (Xtilde(1:2:2*num_node-1))'];
branch=[branch(:,1:2)  g4];

disp( '最优潮流计算结果：');
disp('====================================================');
disp('                  对偶因子Gap                 ');
disp('====================================================');
disp('   迭代次数         Gap');
for i=1:num_iteration
    fprintf('       ');
     fprintf('%2d',i);
fprintf('          ');
     fprintf('%2d',gap_record(i,:));
      fprintf('\n');
end


disp('====================================================');
disp('                  有功无功电源出力                 ');
disp('====================================================');
disp('   母线序号    有功出力    无功出力     燃料费用');
disp(source);


disp('====================================================');
disp('                     节点电压相量                  ');
disp('====================================================');
disp('母线序号        电压幅值           电压相角');
disp(voltage);
 
disp('====================================================');
disp('                    支路有功功率                     ');
disp('====================================================');
disp('    从           到            功率'); 
disp(branch);