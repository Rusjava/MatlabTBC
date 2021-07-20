function u=itteration2D_mod_stable(u, c, P, beta0, cntn, MMAX, edge1, edge2, beta)
 
%Top and bottom bondary conditions
SS=-edge1.'*beta;
SS1=-edge2.'*beta;
   
% Boundary condition at the lower boundary
Q(1)=-(u(3)-c*u(2)+u(1)-beta0*SS)/2;
     
% Preparation solving
for cntm=1:MMAX
    Q(cntm+1)=-(Q(cntm)+u(cntm+2)-c*u(cntm+1)+u(cntm))*P(cntm+1);
end
  
% Boundary condition at the upper boundary
u(MMAX+2)=(beta0*SS1+Q(MMAX)-(P(MMAX)+beta0)*Q(MMAX+1))/(1-(beta0+P(MMAX))*P(MMAX+1));
   
%Solving the system
for cntm=MMAX+1:-1:1
      u(cntm)=Q(cntm)-P(cntm)*u(cntm+1);
end