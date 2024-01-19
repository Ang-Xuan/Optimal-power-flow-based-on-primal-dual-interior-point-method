function Y = solveY(branch,bus, num_node, num_branch);
Y=zeros(num_node,num_node);

for k=1:num_branch                    
    I=branch(k,1);  
    J=branch(k,2);   
    K=branch(k,9);                                               %变压器支路变比
    
    Zs=branch(k,3)+j*branch(k,4);                                %支路阻抗
    if (Zs~=0)
    Ys=1/Zs;
    end
     Yp=j*branch(k,5);                                            %支路导纳
    
	if (K==0)                                                     %若支路为普通支路
	   Y(I,I)=Y(I,I)+Ys+Yp/2;
	   Y(J,J)=Y(J,J)+Ys+Yp/2;
	   Y(I,J)=Y(I,J)-Ys;
       Y(J,I)=Y(I,J); 
    else                                                                 %若支路为变压器支路,采用变压器Π形等效电路	
       Y(I,I)=Y(I,I)+Ys/K/K;       
       Y(J,J)=Y(J,J)+Ys+Yp;                                          
       Y(I,J)=Y(I,J)-Ys/K;           
       Y(J,I)=Y(I,J);
    end
    
end

for k=1:num_node                                                     %补充节点对地导纳
      Y(k,k)=Y(k,k)+bus(k,5)+j*bus(k,6);
end

end
