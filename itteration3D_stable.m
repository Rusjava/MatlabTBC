function u=itteration3D_stable(u, SX, SY, SZ, cci, beta0, MMAX, uleft, uright, ubottom, utop, beta);

   SH=zeros((MMAX+2)^2-4,1);
   SV=-ubottom*beta;
   SV1=-utop*beta;
   SS=-uleft*beta;
   SS1=-uright*beta;
   
   %Definition of sparse matrix
   SH(1:MMAX)=-beta0*SV(2:MMAX+1);
   for ii=1:MMAX
     SH(MMAX+(MMAX+2)*(ii-1)+(1:MMAX+2))=u(1:MMAX+2,ii+2)-cci*u(1:MMAX+2,ii+1)+u(1:MMAX+2,ii);
     SH(MMAX+(MMAX+2)*(ii-1)+1)=-beta0*SS(ii+1);
     SH(2*MMAX+(MMAX+2)*(ii-1)+2)=beta0*SS1(ii+1);
   end
   SH(MMAX+(MMAX+2)*MMAX+(1:MMAX))=beta0*SV1(2:MMAX+1);
  
   % Creating the sparse matrix and solving the system
   SR=sparse(SY,SX,SZ)\SH;
   
   %Determining u
   u(2:MMAX+1,1)=SR(1:MMAX);
   for ii=1:MMAX
     u(1:MMAX+2, ii+1)=SR(MMAX+(MMAX+2)*(ii-1)+(1:MMAX+2));
   end
   u(2:MMAX+1, MMAX+2)=SR(MMAX+(MMAX+2)*MMAX+(1:MMAX));