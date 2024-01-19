function delta_X = delta(H_,dg_dx,dh_dx,L_Z,U_W,b,z,l,u,w,x,y)
b_z = b(1:length(z));
b_l = b(length(z)+1:length([z l]));
b_w = b(length([z l])+1:length([z l w]));
b_u = b(length([z l w])+1:length([z l w u]));
b_x = b(length([z l w u])+1:length([z l w u x]));
b_y = b(length([z l w u x])+1:length([z l w u x y]));

delta_xy = [H_, dh_dx; dh_dx', zeros(length(y),length(y))]\[b_x; b_y];
delta_x = delta_xy(1:length(x));
delta_y = delta_xy(length(x)+1:length([x y]));
delta_l = b_l+dg_dx'*delta_x;
delta_u = b_u-dg_dx'*delta_x;
delta_z = b_z-L_Z*delta_l;
delta_w = b_w-U_W*delta_u;
%% …Ë÷√≤Ω≥§
alpha_p = 1; 
alpha_d = 1;
for k = 1:length(l)
   if (delta_l(k)<0 && -l(k)/delta_l(k)<alpha_p)
      alpha_p = -l(k)/delta_l(k); 
   end
   if (delta_u(k)<0 && -u(k)/delta_u(k)<alpha_p)
      alpha_p = -u(k)/delta_u(k); 
   end
   if (delta_z(k)<0 && -z(k)/delta_z(k)<alpha_d)
      alpha_d = -z(k)/delta_z(k); 
   end
   if (delta_w(k)>0 && -w(k)/delta_w(k)<alpha_d)
      alpha_d = -w(k)/delta_w(k); 
   end
end
alpha_p = 0.9995*alpha_p;
alpha_d = 0.9995*alpha_d;
delta_X = [alpha_d*delta_z; alpha_p*delta_l; alpha_d*delta_w; alpha_p*delta_u; alpha_p*delta_x; alpha_d*delta_y];
end