clc;clear all;warning off;
opfcase;
%% ������������
num_node=size(bus,1);                              %�ڵ����
num_branch=size(branch,1);                         %֧·��
num_gen=size(gen,1);                               %���������
bus_3=find(bus(:,2)==3);                           %ƽ��ڵ�
bus_2=find(bus(:,2)==2);                           %PV�ڵ�
bus_1=find(bus(:,2)==1);                           %PQ�ڵ�
num_PQ=size(bus_1,1);
num_PV=size(bus_2,1);

num_equa = 2*num_node;                                 %��ʽԼ������
num_inequa = 2*num_gen+num_node+num_branch;            %����ʽԼ������

a2=gencost(:,5);  a1=gencost(:,6);  a0=gencost(:,7);    %�������Զ���ʽϵ��
A2=diag(a2);  A1=diag(a1);

g_u = [gen(:,9); gen(:,4); bus(:,12); branch(:,6)];
%������й��Ͻ� ������޹��Ͻ� �ڵ��ѹ���� ֧·���޹���
g_l = [gen(:,10); gen(:,5); bus(:,13); -branch(:,6)];
%������й��½� ������޹��½� �ڵ��ѹ���� ֧·���޹���

Xtilde(1:2:num_node*2-1) = (bus(:,9))';               %�ڵ��ʼ���
Xtilde(2:2:num_node*2) = (bus(:,8))';                 %�ڵ��ʼ��ֵ

u=1*ones(num_inequa,1);                               %�������ճ��ӳ�ʼ��ֵ
l=1*ones(num_inequa,1);
z=1*ones(1,num_inequa); 
w=-0.5*ones(1,num_inequa); 
y=1E-10*ones(1,num_equa); 
y(2:2:num_equa)=-1*y(2:2:num_equa);            

epsi=1E-6;
sigma=0.1;
max_iteration=50; 
gap_record=zeros(max_iteration,1);
%% ��ڵ㵼�ɾ���
Y=solveY(branch,bus,num_node,num_branch);
%% ��ʼ����
for num_iteration = 1:max_iteration  
gap=l'*z'-u'*w';
gap_record(num_iteration) = gap;
    if (gap<=epsi)
        disp('iterations completed');
        break; 
    end

Pg=(gen(:,9)+gen(:,10))/2;                                     % ������й�����
Qr=(gen(:,4)+gen(:,5))/2;                                      % ������޹�����
g1=Pg;                                                         % ������й�����g1
g2=Qr;                                                         % ������޹�����g2
g3=(Xtilde(2:2:num_node*2))';                                  % ���нڵ��ѹ��ֵg3   
    for k=1:num_branch                                         % ��·����g4
        i=branch(k,1);
        j=branch(k,2);
        cita=Xtilde(i*2-1)-Xtilde(j*2-1);
     g4(k,1) = Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*cos(cita)+imag(Y(i,j))*sin(cita))-...
               Xtilde(2*i)*Xtilde(2*i)*real(Y(i,j));     
    end
g=[g1; g2; g3; g4];                                         %  num_inequa*1
    
for i=1:num_node                                              %���ڵ㾻ע��
P_ejec(i,1)=-bus(i,3); 
Q_ejec(i,1)=-bus(i,4);
  if bus(i,2)~=1
P_ejec(i,1)=P_ejec(i,1)+gen(gen(:,1)==i,2); 
Q_ejec(i,1)=Q_ejec(i,1)+gen(gen(:,1)==i,3);
  end
end  

h=zeros(num_node*2,1);                                        %��y��ƫ������ʽԼ����num_equa*1
for i=1:num_node
     for j=1:num_node
            cita=Xtilde(i*2-1)-Xtilde(j*2-1);
            h(2*i-1)=h(2*i-1)-Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*cos(cita)+imag(Y(i,j))*sin(cita));
            h(2*i)=h(2*i)-Xtilde(2*i)*Xtilde(2*j)*(real(Y(i,j))*sin(cita)-imag(Y(i,j))*cos(cita));
     end
            h(2*i-1)=h(2*i-1)+P_ejec(i);
            h(2*i)=h(2*i)+Q_ejec(i);
end
    %% ��ϵ������
    x=[Pg;Qr;Xtilde'];                                             % ״̬������
    X=[z';l;w';u;x;y'];                                       % ������
    [dh_dx,dg_dx,H_,L_Z,U_W]=coe(Y,bus,gen,branch,X,A2);
    
    Ly=h;                                                          % Lagrange����ƫ��
    Lz=g-l-g_l;
    Lw=g+u-g_u;
    L=diag(l); U=diag(u); Z=diag(z); W=diag(w);
    miu=gap*sigma/2/num_inequa;
    L_miu_l=L*Z*ones(num_inequa,1)-miu*ones(num_inequa,1);
    L_miu_u=U*W*ones(num_inequa,1)+miu*ones(num_inequa,1);
    df_dx=[2*A2*Pg+a1;zeros(num_gen,1);zeros(length(Xtilde),1)];               %Ŀ�꺯���ݶ�ʸ��num_inequa*1
    Lx_=df_dx-dh_dx*y'-dg_dx*(z+w)'+dg_dx*(L\(L_miu_l+Z*Lz)+U\(L_miu_u-W*Lw));    
    b=[-L\L_miu_l;Lz;-U\L_miu_u;-Lw;Lx_;-Ly];   
    %% �����󷽳�
    delta_X=delta(H_,dg_dx,dh_dx,L_Z,U_W,b,z,l',u',w,x',y); 
    X=X+delta_X;    
    % ���µ�����
    z=(X(1:num_inequa))';
    l=X(num_inequa+1:2*num_inequa);
    w=(X(2*num_inequa+1:3*num_inequa))';
    u=X(3*num_inequa+1:4*num_inequa);
    x=X(4*num_inequa+1:5*num_inequa);
    y=(X(5*num_inequa+1:5*num_inequa+num_equa))';
    % ����״̬��
    Pg=(x(1:num_gen));
    Qr=(x(num_gen+1:2*num_gen));
    Xtilde=(x(2*num_gen+1:2*num_gen+2*num_node))';   
end
if (num_iteration>=max_iteration)
  disp('OPF is not convergent');  
end
%% �������
semilogy(gap_record(1:num_iteration));
grid on;
title('��������');
xlabel('��������'); 
ylabel('Gap');

%��Ǳ任
for k=1:num_node
    Xtilde(2*k-1)=Xtilde(2*k-1)-Xtilde(2*num_node-1);
end

source=[gen(:,1) Pg Qr A2*(Pg.*Pg)+A1*Pg+a0];
voltage=[bus(:,1)  (Xtilde(2:2:2*num_node))'  (Xtilde(1:2:2*num_node-1))'];
branch=[branch(:,1:2)  g4];

disp( '���ų�����������');
disp('====================================================');
disp('                  ��ż����Gap                 ');
disp('====================================================');
disp('   ��������         Gap');
for i=1:num_iteration
    fprintf('       ');
     fprintf('%2d',i);
fprintf('          ');
     fprintf('%2d',gap_record(i,:));
      fprintf('\n');
end


disp('====================================================');
disp('                  �й��޹���Դ����                 ');
disp('====================================================');
disp('   ĸ�����    �й�����    �޹�����     ȼ�Ϸ���');
disp(source);


disp('====================================================');
disp('                     �ڵ��ѹ����                  ');
disp('====================================================');
disp('ĸ�����        ��ѹ��ֵ           ��ѹ���');
disp(voltage);
 
disp('====================================================');
disp('                    ֧·�й�����                     ');
disp('====================================================');
disp('    ��           ��            ����'); 
disp(branch);