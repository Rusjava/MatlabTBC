function [SY, SX, SZ]=Spmatrix2D_stable(MMAX,beta0,ci);

DIM=5*MMAX^2+12*MMAX;
DIM1=MMAX^2+3*MMAX;
SY=zeros(DIM,1);
SX=zeros(DIM,1);
SZ=zeros(DIM,1);

%The bottom
SY(1:MMAX)=(1:MMAX);
SX(1:MMAX)=(1:MMAX);
SZ(1:MMAX)=-1;

SY(MMAX+(1:MMAX))=(1:MMAX);
SX(MMAX+(1:MMAX))=(1:MMAX)+MMAX+1;
SZ(MMAX+(1:MMAX))=-beta0;

SY(2*MMAX+(1:MMAX))=(1:MMAX);
SX(2*MMAX+(1:MMAX))=(1:MMAX)+2*MMAX+3;
SZ(2*MMAX+(1:MMAX))=1;

%The middle
for ii=1:MMAX
   SY(3*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX))=MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(3*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX))=MMAX+(MMAX+2)*(ii-2)+(1:MMAX)+1; 
   SZ(3*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX))=-1;
   
   SY(4*MMAX+(5*MMAX+6)*(ii-1)+1)=MMAX+(MMAX+2)*(ii-1)+1;
   SX(4*MMAX+(5*MMAX+6)*(ii-1)+1)=MMAX+(MMAX+2)*(ii-1)+1;
   SZ(4*MMAX+(5*MMAX+6)*(ii-1)+1)=-1;
   SY(4*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+1)=MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(4*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+1)=MMAX+(MMAX+2)*(ii-1)+(1:MMAX);
   SZ(4*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+1)=-1;
   SY(5*MMAX+(5*MMAX+6)*(ii-1)+2)=MMAX+(MMAX+2)*ii;
   SX(5*MMAX+(5*MMAX+6)*(ii-1)+2)=MMAX+(MMAX+2)*ii-2;
   SZ(5*MMAX+(5*MMAX+6)*(ii-1)+2)=-1;
   
   SY(5*MMAX+(5*MMAX+6)*(ii-1)+3)=MMAX+(MMAX+2)*(ii-1)+1;
   SX(5*MMAX+(5*MMAX+6)*(ii-1)+3)=MMAX+(MMAX+2)*(ii-1)+2;
   SZ(5*MMAX+(5*MMAX+6)*(ii-1)+3)=-beta0;
   SY(5*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+3)=MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(5*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+3)=MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1; 
   SZ(5*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+3)=ci;
   SY(6*MMAX+(5*MMAX+6)*(ii-1)+4)=MMAX+(MMAX+2)*ii;
   SX(6*MMAX+(5*MMAX+6)*(ii-1)+4)=MMAX+(MMAX+2)*ii-1;
   SZ(6*MMAX+(5*MMAX+6)*(ii-1)+4)=beta0;
   
   SY(6*MMAX+(5*MMAX+6)*(ii-1)+5)=MMAX+(MMAX+2)*(ii-1)+1;
   SX(6*MMAX+(5*MMAX+6)*(ii-1)+5)=MMAX+(MMAX+2)*(ii-1)+3;
   SZ(6*MMAX+(5*MMAX+6)*(ii-1)+5)=1;
   SY(6*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+5)=MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(6*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+5)=MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+2;
   SZ(6*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+5)=-1;
   SY(7*MMAX+(5*MMAX+6)*(ii-1)+6)=MMAX+(MMAX+2)*ii;
   SX(7*MMAX+(5*MMAX+6)*(ii-1)+6)=MMAX+(MMAX+2)*ii;
   SZ(7*MMAX+(5*MMAX+6)*(ii-1)+6)=1;
   
   SY(7*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+6)=MMAX+(MMAX+2)*(ii-1)+(1:MMAX)+1;
   SX(7*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+6)=MMAX+(MMAX+2)*ii+(1:MMAX)+1;
   SZ(7*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+6)=-1;   
end

% The second and MMAX+1 rows
SX(3*MMAX+(1:MMAX))=(1:MMAX); 
SX(7*MMAX+(5*MMAX+6)*(MMAX-1)+(1:MMAX)+6)=MMAX+(MMAX+2)*MMAX+(1:MMAX); 

%The top
SY((1:MMAX)+DIM-3*MMAX)=(1:MMAX)+DIM1;
SX((1:MMAX)+DIM-3*MMAX)=(1:MMAX)+DIM1-2*MMAX-3;
SZ((1:MMAX)+DIM-3*MMAX)=-1;

SY((1:MMAX)+DIM-2*MMAX)=(1:MMAX)+DIM1;
SX((1:MMAX)+DIM-2*MMAX)=(1:MMAX)+DIM1-MMAX-1;
SZ((1:MMAX)+DIM-2*MMAX)=beta0;

SY((1:MMAX)+DIM-MMAX)=(1:MMAX)+DIM1;
SX((1:MMAX)+DIM-MMAX)=(1:MMAX)+DIM1;
SZ((1:MMAX)+DIM-MMAX)=1;