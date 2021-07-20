function [SY, SX, SZ]=Spmatrix3D_stable(MMAX,beta0,ci);

DIM=MMAX^3+6*MMAX^2;
DIM1=MMAX^2+4*MMAX;
DIM2=7*MMAX^3+18*MMAX^2;
DIM3=7*MMAX^2+12*MMAX;
DIM4=7*MMAX^3+15*MMAX^2;

SY=zeros(DIM2,1);
SX=zeros(DIM2,1);
SZ=zeros(DIM2,1);

%The forward and back of 3D
for ii=1:MMAX
  %forward
  SY(MMAX*(ii-1)+(1:MMAX))=MMAX*(ii-1)+(1:MMAX);
  SX(MMAX*(ii-1)+(1:MMAX))=MMAX*(ii-1)+(1:MMAX);
  SZ(MMAX*(ii-1)+(1:MMAX))=-1;  
  
  SY(MMAX^2+MMAX*(ii-1)+(1:MMAX))=MMAX*(ii-1)+(1:MMAX);
  SX(MMAX^2+MMAX*(ii-1)+(1:MMAX))=MMAX^2+MMAX+(MMAX+2)*(ii-1)+1+(1:MMAX);
  SZ(MMAX^2+MMAX*(ii-1)+(1:MMAX))=-beta0; 
  
  SY(2*MMAX^2+MMAX*(ii-1)+(1:MMAX))=MMAX*(ii-1)+(1:MMAX);
  SX(2*MMAX^2+MMAX*(ii-1)+(1:MMAX))=MMAX^2+DIM1+MMAX+(MMAX+2)*(ii-1)+1+(1:MMAX);
  SZ(2*MMAX^2+MMAX*(ii-1)+(1:MMAX))=1; 
  
  %back
  SY(DIM4+MMAX*(ii-1)+(1:MMAX))=MMAX^2+DIM1*MMAX+MMAX*(ii-1)+(1:MMAX);
  SX(DIM4+MMAX*(ii-1)+(1:MMAX))=MMAX^2+DIM1*(MMAX-2)+MMAX+(MMAX+2)*(ii-1)+1+(1:MMAX);
  SZ(DIM4+MMAX*(ii-1)+(1:MMAX))=-1;  
  
  SY(DIM4+MMAX^2+MMAX*(ii-1)+(1:MMAX))=MMAX^2+DIM1*MMAX+MMAX*(ii-1)+(1:MMAX);
  SX(DIM4+MMAX^2+MMAX*(ii-1)+(1:MMAX))=MMAX^2+DIM1*(MMAX-1)+MMAX+(MMAX+2)*(ii-1)+1+(1:MMAX);
  SZ(DIM4+MMAX^2+MMAX*(ii-1)+(1:MMAX))=beta0; 
  
  SY(DIM4+2*MMAX^2+MMAX*(ii-1)+(1:MMAX))=MMAX^2+DIM1*MMAX+MMAX*(ii-1)+(1:MMAX);
  SX(DIM4+2*MMAX^2+MMAX*(ii-1)+(1:MMAX))=MMAX^2+DIM1*MMAX+MMAX*(ii-1)+(1:MMAX);
  SZ(DIM4+2*MMAX^2+MMAX*(ii-1)+(1:MMAX))=1; 
end

%The middle of 3D
for jj=1:MMAX
  DIMj=4*MMAX^2+DIM3*(jj-1);
  DIMj1=MMAX^2+DIM1*(jj-1);
  DIMj2=MMAX^2+DIM1*jj;
  DIMj3=2*MMAX^2+DIM3*jj;
  
  for ii=1:MMAX
   %Forward
   SY(3*MMAX^2+DIM3*(jj-1)+MMAX*(ii-1)+(1:MMAX))=DIMj1+MMAX+(MMAX+2)*(ii-1)+1+(1:MMAX);
   SX(3*MMAX^2+DIM3*(jj-1)+MMAX*(ii-1)+(1:MMAX))=DIMj1-DIM1+MMAX+(MMAX+2)*(ii-1)+1+(1:MMAX);
   SZ(3*MMAX^2+DIM3*(jj-1)+MMAX*(ii-1)+(1:MMAX))=-1;
   %back
   SY(DIMj3+MMAX*(ii-1)+(1:MMAX))=DIMj1+MMAX+(MMAX+2)*(ii-1)+1+(1:MMAX);
   SX(DIMj3+MMAX*(ii-1)+(1:MMAX))=DIMj2+MMAX+(MMAX+2)*(ii-1)+1+(1:MMAX);
   SZ(DIMj3+MMAX*(ii-1)+(1:MMAX))=-1;
  end
  
 %The bottom 
  SY(DIMj+(1:MMAX))=DIMj1+(1:MMAX);
  SX(DIMj+(1:MMAX))=DIMj1+(1:MMAX);
  SZ(DIMj+(1:MMAX))=-1;

  SY(DIMj+MMAX+(1:MMAX))=DIMj1+(1:MMAX);
  SX(DIMj+MMAX+(1:MMAX))=DIMj1+(1:MMAX)+MMAX+1;
  SZ(DIMj+MMAX+(1:MMAX))=-beta0;

  SY(DIMj+2*MMAX+(1:MMAX))=DIMj1+(1:MMAX);
  SX(DIMj+2*MMAX+(1:MMAX))=DIMj1+(1:MMAX)+2*MMAX+3;
  SZ(DIMj+2*MMAX+(1:MMAX))=1;

%The middle
  for ii=1:MMAX
   SY(DIMj+3*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX))=DIMj1+MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(DIMj+3*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX))=DIMj1+MMAX+(MMAX+2)*(ii-2)+(1:MMAX)+1; 
   SZ(DIMj+3*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX))=-1;
   
   SY(DIMj+4*MMAX+(5*MMAX+6)*(ii-1)+1)=DIMj1+MMAX+(MMAX+2)*(ii-1)+1;
   SX(DIMj+4*MMAX+(5*MMAX+6)*(ii-1)+1)=DIMj1+MMAX+(MMAX+2)*(ii-1)+1;
   SZ(DIMj+4*MMAX+(5*MMAX+6)*(ii-1)+1)=-1;
   SY(DIMj+4*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+1)=DIMj1+MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(DIMj+4*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+1)=DIMj1+MMAX+(MMAX+2)*(ii-1)+(1:MMAX);
   SZ(DIMj+4*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+1)=-1;
   SY(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+2)=DIMj1+MMAX+(MMAX+2)*ii;
   SX(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+2)=DIMj1+MMAX+(MMAX+2)*ii-2;
   SZ(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+2)=-1;
   
   SY(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+3)=DIMj1+MMAX+(MMAX+2)*(ii-1)+1;
   SX(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+3)=DIMj1+MMAX+(MMAX+2)*(ii-1)+2;
   SZ(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+3)=-beta0;
   SY(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+3)=DIMj1+MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+3)=DIMj1+MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1; 
   SZ(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+3)=ci;
   SY(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+4)=DIMj1+MMAX+(MMAX+2)*ii;
   SX(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+4)=DIMj1+MMAX+(MMAX+2)*ii-1;
   SZ(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+4)=beta0;
   
   SY(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+5)=DIMj1+MMAX+(MMAX+2)*(ii-1)+1;
   SX(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+5)=DIMj1+MMAX+(MMAX+2)*(ii-1)+3;
   SZ(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+5)=1;
   SY(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+5)=DIMj1+MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+5)=DIMj1+MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+2;
   SZ(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+5)=-1;
   SY(DIMj+7*MMAX+(5*MMAX+6)*(ii-1)+6)=DIMj1+MMAX+(MMAX+2)*ii;
   SX(DIMj+7*MMAX+(5*MMAX+6)*(ii-1)+6)=DIMj1+MMAX+(MMAX+2)*ii;
   SZ(DIMj+7*MMAX+(5*MMAX+6)*(ii-1)+6)=1;
   
   SY(DIMj+7*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+6)=DIMj1+MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(DIMj+7*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+6)=DIMj1+MMAX+(MMAX+2)*ii+(1:MMAX)+1;
   SZ(DIMj+7*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+6)=-1;   
  end

% The second and MMAX+1 rows
  SX(DIMj+3*MMAX+(1:MMAX))=DIMj1+(1:MMAX); 
  SX(DIMj+7*MMAX+(5*MMAX+6)*(MMAX-1)+(1:MMAX)+6)=DIMj1+MMAX+(MMAX+2)*MMAX+(1:MMAX); 

%The top
  SY(DIMj3-3*MMAX+(1:MMAX))=DIMj2-MMAX+(1:MMAX);
  SX(DIMj3-3*MMAX+(1:MMAX))=DIMj2-3*MMAX+(1:MMAX)-3;
  SZ(DIMj3-3*MMAX+(1:MMAX))=-1;

  SY(DIMj3-2*MMAX+(1:MMAX))=DIMj2-MMAX+(1:MMAX);
  SX(DIMj3-2*MMAX+(1:MMAX))=DIMj2-2*MMAX+(1:MMAX)-1;
  SZ(DIMj3-2*MMAX+(1:MMAX))=beta0;

  SY(DIMj3-MMAX+(1:MMAX))=DIMj2-MMAX+(1:MMAX);
  SX(DIMj3-MMAX+(1:MMAX))=DIMj2-MMAX+(1:MMAX);
  SZ(DIMj3-MMAX+(1:MMAX))=1;
end
% The second and MMAX+1 planes
for ii=1:MMAX 
   SX(3*MMAX^2+MMAX*(ii-1)+(1:MMAX))=MMAX*(ii-1)+(1:MMAX);
   SX(2*MMAX^2+DIM3*MMAX+MMAX*(ii-1)+(1:MMAX))=MMAX^2+DIM1*MMAX+MMAX*(ii-1)+(1:MMAX);
end