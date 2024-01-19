function Y = solveY(branch,bus, num_node, num_branch);
Y=zeros(num_node,num_node);

for k=1:num_branch                    
    I=branch(k,1);  
    J=branch(k,2);   
    K=branch(k,9);                                               %��ѹ��֧·���
    
    Zs=branch(k,3)+j*branch(k,4);                                %֧·�迹
    if (Zs~=0)
    Ys=1/Zs;
    end
     Yp=j*branch(k,5);                                            %֧·����
    
	if (K==0)                                                     %��֧·Ϊ��֧ͨ·
	   Y(I,I)=Y(I,I)+Ys+Yp/2;
	   Y(J,J)=Y(J,J)+Ys+Yp/2;
	   Y(I,J)=Y(I,J)-Ys;
       Y(J,I)=Y(I,J); 
    else                                                                 %��֧·Ϊ��ѹ��֧·,���ñ�ѹ�����ε�Ч��·	
       Y(I,I)=Y(I,I)+Ys/K/K;       
       Y(J,J)=Y(J,J)+Ys+Yp;                                          
       Y(I,J)=Y(I,J)-Ys/K;           
       Y(J,I)=Y(I,J);
    end
    
end

for k=1:num_node                                                     %����ڵ�Եص���
      Y(k,k)=Y(k,k)+bus(k,5)+j*bus(k,6);
end

end
