clear all%------------------------------------------3D Schrodinger equation with 2D propagation in facets with stable discrete TBC 
close all

alp0=0;
ww=0.1;
RMIN=80; %------------------------well radius
RMIN2=70;
RMAX=100;%------------Maximin r
ZMAX=2000;%----------------simulation length
NIT=300;%---------------------number of itterations
precision=1e-6;
alp1=-ww/(RMIN-RMIN2);

h=1.66;%----------------------------- <<lam*f/R
tau=20;%----------------------------- step

    sprsn=2;%----------------------------ARRAY THINNING (long range)
    sprsm=1;%----------------------------ARRAY THINNING

MMAX=round(2*RMAX/h);
muMAX=floor((MMAX+2)/sprsm);
MMIN=round(RMIN/h);
MMIN2=round(RMIN2/h);
NMAX=round(ZMAX/tau);
centr=(MMAX+3)/2;
rad2=MMIN^2;
rad22=MMIN2^2;

r=(0:MMAX+1)*h;
z(1:NMAX+1)=tau*(0:NMAX);

%Sizes of matrices and other auxiliary parameters
alp=zeros(MMAX+2,MMAX+2,MMAX+2);

%------------------------------Initial conditions
u0=ones(MMAX+2,MMAX+2,MMAX+2);

beam=0; %-----------------Illumination
switch beam
    case 0 %---------------------------------------GAUSSIAN
      brad=MMIN2^2;
      WAIST=brad/25;
      for cntm=1:MMAX+2
        for cntk=1:MMAX+2
          u(cntm,cntk,1:MMAX+2)=(sign(brad-(cntk-centr).^2-(cntm-centr).^2-((1:MMAX+2)-centr).^2)+1)/2.*exp(-(((1:MMAX+2)-centr).^2+(cntm-centr)^2+(cntk-centr)^2)/WAIST);
        end     
      end  
    case 1 %------------------------------------------------flat top
       for cntm=1:MMAX+2
        for cntk=1:MMAX+2
          u(cntm,cntk,1:MMAX+2)=(sign(rad22-(cntk-centr).^2-(cntm-centr).^2-((1:MMAX+2)-centr).^2)+1)/2;
        end     
      end  
end

%----------------------------------------------MARCHING
%Boundary conditions
utop(1:MMAX+2,1:MMAX+2,1)=u(:,:,MMAX+1);
ubottom(1:MMAX+2,1:MMAX+2,1)=u(:,:,2);
uleft(1:MMAX+2,1:MMAX+2,1)=u(2,:,:); 
uright(1:MMAX+2,1:MMAX+2,1)=u(MMAX+1,:,:);
uforward(1:MMAX+2,1:MMAX+2,1)=u(:,2,:); 
uback(1:MMAX+2,1:MMAX+2,1)=u(:,MMAX+1,:);

PLANE12(1:MMAX+2,1,1)=u(:,2,MMAX+1);
PLANE56(1:MMAX+2,1,1)=u(:,MMAX+1,MMAX+1);
PLANE15(1:MMAX+2,1,1)=u(2,:,MMAX+1);
PLANE26(1:MMAX+2,1,1)=u(MMAX+1,:,MMAX+1);
PLANE34(1:MMAX+2,1,1)=u(:,2,2);
PLANE78(1:MMAX+2,1,1)=u(:,MMAX+1,2);
PLANE37(1:MMAX+2,1,1)=u(2,:,2);
PLANE48(1:MMAX+2,1,1)=u(MMAX+1,:,2);
PLANE13(1:MMAX+2,1,1)=u(2,2,:);
PLANE24(1:MMAX+2,1,1)=u(MMAX+1,2,:);
PLANE57(1:MMAX+2,1,1)=u(2,MMAX+1,:);
PLANE68(1:MMAX+2,1,1)=u(MMAX+1,MMAX+1,:);

EDGE1(1,1,1)=u(2,2,MMAX+1); EDGE2(1,1,1)=u(MMAX+1,2,MMAX+1); EDGE3(1,1,1)=u(2,2,2); EDGE4(1,1,1)=u(MMAX+1,2,2);
EDGE5(1,1,1)=u(2,MMAX+1,MMAX+1); EDGE6(1,1,1)=u(MMAX+1,MMAX+1,MMAX+1); EDGE7(1,1,1)=u(2,MMAX+1,2); EDGE8(1,1,1)=u(MMAX+1,MMAX+1,2);

% Matrices of results
zplot(1)=0;
rplot(1:muMAX)=r(sprsm*(1:muMAX));
uplot(1:muMAX,1:muMAX,1:muMAX,1)=u(sprsm*(1:muMAX),sprsm*(1:muMAX),sprsm*(1:muMAX));  

%Common auxiliary parameters
I=eye(MMAX+2,MMAX+2);
DG=diag(ones(MMAX+1,1),1)+diag(ones(MMAX+1,1),-1);
nuu=0;
Dalp=h^2*(alp1-alp0)/2;
c0=2*i*h^2/tau;

ci=6-c0;
ci_2D=4-c0;
cci=I*(6+c0)-DG;
cci_2D=I*(4+c0)-DG;
c(1)=2-c0; c(2)=2+c0;
beta0=-i*2*sqrt(c0-c0^2/4);
phi=-1/2-(-1).^(0:NMAX)+((-1).^(0:NMAX))/2.*((1+c0/4)/(1-c0/4)).^(1:NMAX+1);
beta(1)=phi(1);
P(1)=-(c(1)-beta0)/2; 
for cntm=1:MMAX
    P(cntm+1)=-1/(c(1)+P(cntm));
end
%--------------------------------------------------------------------------
% Preparation of sparse matrixes
[SY, SX, SZ]=Spmatrix3D_stable(MMAX,0,0);

DIM=MMAX^3+6*MMAX^2;
DIM1=MMAX^2+4*MMAX;
DIM3=7*MMAX^2+12*MMAX;
DIM4=7*MMAX^3+15*MMAX^2;

SH=zeros(DIM,1);
%--------------------------------------------------------------------------

% Preparation of 2D sparse matrixes
[SY_2D, SX_2D, SZ_2D]=Spmatrix2D_stable(MMAX,beta0,ci_2D);

%--------------------------------------------------------------------------
potential=4;%-------------------------------------potential
tic
for cntn=1:NMAX
    alp=zeros(MMAX+2,MMAX+2,MMAX+2);
      switch potential
        case 1 % Sphere
          for cntm=1:MMAX+2
            for cntk=1:MMAX+2
              alp(cntm,cntk,1:MMAX+2)=h^2*alp0+(sign(rad2-(cntk-centr).^2-(cntm-centr).^2-((1:MMAX+2)-centr).^2)+1)*Dalp;
            end
          end  
        case 2 % Spherical barrier
          for cntm=1:MMAX+2
            for cntk=1:MMAX+2
              alp(cntm,cntk,1:MMAX+2)=h^2*alp0+(sign(rad2-(cntk-centr).^2-(cntm-centr).^2-((1:MMAX+2)-centr).^2)+1).*(1-sign(rad22-(cntk-centr).^2-(cntm-centr).^2-((1:MMAX+2)-centr).^2))*Dalp/2;
            end
          end  

        case 3 % Semi-sphere
          for cntm=1:MMAX+2
            for cntk=1:MMAX+2
              alp(cntm,cntk,1:MMAX+2)=h^2*alp0+(sign(rad2-(cntk-centr).^2-(cntm-centr).^2-((1:MMAX+2)-centr).^2)+1)*(sign(cntk-centr)+1)*Dalp/2;
            end
          end 
        
        case 4 % Semi-spherical barrier
          for cntm=1:MMAX+2
            for cntk=1:MMAX+2
              alp(cntm,cntk,1:MMAX+2)=h^2*alp0+(sign(rad2-(cntk-centr).^2-(cntm-centr).^2-((1:MMAX+2)-centr).^2)+1).*(1-sign(rad22-(cntk-centr).^2-(cntm-centr).^2-((1:MMAX+2)-centr).^2))*(sign(cntk-centr)+1)*Dalp/4;
            end
          end 
      end

   % Auxiliry parameters
   beta1=flipud(beta.');
    
  %Propagation on edges
   parfor s=1:cntn
      for t=1:cntn
        PLANE12(:,t,s)=itteration2D_mod_stable(PLANE12(:,t,s), c(2), P, beta0, cntn, MMAX, EDGE1(:,t,s), EDGE2(:,t,s), beta1);
        PLANE56(:,t,s)=itteration2D_mod_stable(PLANE56(:,t,s), c(2), P, beta0, cntn, MMAX, EDGE5(:,t,s), EDGE6(:,t,s), beta1);
        PLANE34(:,t,s)=itteration2D_mod_stable(PLANE34(:,t,s), c(2), P, beta0, cntn, MMAX, EDGE3(:,t,s), EDGE4(:,t,s), beta1);
        PLANE78(:,t,s)=itteration2D_mod_stable(PLANE78(:,t,s), c(2), P, beta0, cntn, MMAX, EDGE7(:,t,s), EDGE8(:,t,s), beta1);
        
        PLANE15(:,t,s)=itteration2D_mod_stable(PLANE15(:,t,s), c(2), P, beta0, cntn, MMAX, EDGE1(t,:,s).', EDGE5(t,:,s).', beta1);
        PLANE26(:,t,s)=itteration2D_mod_stable(PLANE26(:,t,s), c(2), P, beta0, cntn, MMAX, EDGE2(t,:,s).', EDGE6(t,:,s).', beta1);
        PLANE37(:,t,s)=itteration2D_mod_stable(PLANE37(:,t,s), c(2), P, beta0, cntn, MMAX, EDGE3(t,:,s).', EDGE7(t,:,s).', beta1);
        PLANE48(:,t,s)=itteration2D_mod_stable(PLANE48(:,t,s), c(2), P, beta0, cntn, MMAX, EDGE4(t,:,s).', EDGE8(t,:,s).', beta1);
        
        PLANE13(:,t,s)=itteration2D_mod_stable(PLANE13(:,t,s), c(2), P, beta0, cntn, MMAX, squeeze(EDGE3(t,s,:)), squeeze(EDGE1(t,s,:)), beta1);
        PLANE24(:,t,s)=itteration2D_mod_stable(PLANE24(:,t,s), c(2), P, beta0, cntn, MMAX, squeeze(EDGE4(t,s,:)), squeeze(EDGE2(t,s,:)), beta1);
        PLANE57(:,t,s)=itteration2D_mod_stable(PLANE57(:,t,s), c(2), P, beta0, cntn, MMAX, squeeze(EDGE7(t,s,:)), squeeze(EDGE5(t,s,:)), beta1);
        PLANE68(:,t,s)=itteration2D_mod_stable(PLANE68(:,t,s), c(2), P, beta0, cntn, MMAX, squeeze(EDGE8(t,s,:)), squeeze(EDGE6(t,s,:)), beta1);
      end
   end
  
   % Propagation on facets
   parfor s=1:cntn 
       ubottom(:,:,s)=itteration3D_stable(ubottom(:,:,s), SX_2D, SY_2D, SZ_2D, cci_2D, beta0, MMAX, PLANE37(:,:,s), PLANE48(:,:,s), PLANE34(:,:,s), PLANE78(:,:,s), beta1);
       utop(:,:,s)=itteration3D_stable(utop(:,:,s), SX_2D, SY_2D, SZ_2D, cci_2D, beta0, MMAX, PLANE15(:,:,s), PLANE26(:,:,s), PLANE12(:,:,s), PLANE56(:,:,s), beta1);
       uleft(:,:,s)=itteration3D_stable(uleft(:,:,s), SX_2D, SY_2D, SZ_2D, cci_2D, beta0, MMAX, squeeze(PLANE13(:,s,:)), squeeze(PLANE57(:,s,:)), squeeze(PLANE37(:,s,:)), squeeze(PLANE15(:,s,:)), beta1);
       uright(:,:,s)=itteration3D_stable(uright(:,:,s), SX_2D, SY_2D, SZ_2D, cci_2D, beta0, MMAX, squeeze(PLANE24(:,s,:)), squeeze(PLANE68(:,s,:)), squeeze(PLANE48(:,s,:)), squeeze(PLANE26(:,s,:)), beta1);
       uforward(:,:,s)=itteration3D_stable(uforward(:,:,s), SX_2D, SY_2D, SZ_2D, cci_2D, beta0, MMAX, PLANE13(:,:,s), PLANE24(:,:,s), squeeze(PLANE34(:,s,:)), squeeze(PLANE12(:,s,:)), beta1);
       uback(:,:,s)=itteration3D_stable(uback(:,:,s), SX_2D, SY_2D, SZ_2D, cci_2D, beta0, MMAX, PLANE57(:,:,s), PLANE68(:,:,s), squeeze(PLANE78(:,s,:)), squeeze(PLANE56(:,s,:)), beta1); 
   end
   
   SV(1:MMAX+2,1:MMAX+2)=0;
   SV1(1:MMAX+2,1:MMAX+2)=0;
   SS(1:MMAX+2,1:MMAX+2)=0;
   SS1(1:MMAX+2,1:MMAX+2)=0;
   SVV(1:MMAX+2,1:MMAX+2)=0;
   SVV1(1:MMAX+2,1:MMAX+2)=0;
   for jj=1:cntn
      SV(1:MMAX+2,1:MMAX+2)=SV(1:MMAX+2,1:MMAX+2)-ubottom(1:MMAX+2,1:MMAX+2,jj)* beta1(jj);
      SV1(1:MMAX+2,1:MMAX+2)=SV1(1:MMAX+2,1:MMAX+2)-utop(1:MMAX+2,1:MMAX+2,jj)* beta1(jj);
      SS(1:MMAX+2,1:MMAX+2)=SS(1:MMAX+2,1:MMAX+2)-uleft(1:MMAX+2,1:MMAX+2,jj)* beta1(jj);
      SS1(1:MMAX+2,1:MMAX+2)=SS1(1:MMAX+2,1:MMAX+2)-uright(1:MMAX+2,1:MMAX+2,jj)* beta1(jj);
      SVV(1:MMAX+2,1:MMAX+2)=SVV(1:MMAX+2,1:MMAX+2)-uforward(1:MMAX+2,1:MMAX+2,jj)* beta1(jj);
      SVV1(1:MMAX+2,1:MMAX+2)=SVV1(1:MMAX+2,1:MMAX+2)-uback(1:MMAX+2,1:MMAX+2,jj)* beta1(jj);
   end
   
   % Auxiliry parameter  
   beta(cntn+1)=(phi(1:cntn)*flipud(beta.')+phi(cntn+1))/(cntn+1);
     
   %Definition of a sparse matrix
   for ii=1:MMAX
    %Forward surface
     SZ(MMAX^2+MMAX*(ii-1)+(1:MMAX))=-beta0; 
     SH(MMAX*(ii-1)+(1:MMAX))=-beta0*SVV(2:MMAX+1,ii+1);
    %Back surface
     SZ(DIM4+MMAX^2+MMAX*(ii-1)+(1:MMAX))=beta0; 
     SH(MMAX^2+DIM1*MMAX+MMAX*(ii-1)+(1:MMAX))=beta0*SVV1(2:MMAX+1,ii+1);
   end
   % Middle
   for jj=1:MMAX 
    DIMj=4*MMAX^2+DIM3*(jj-1);
    DIMj1=MMAX^2+DIM1*(jj-1);
    DIMj2=DIMj1+DIM1;
    DIMj3=2*MMAX^2+DIM3*jj; 
   
    SZ(DIMj+MMAX+(1:MMAX))=-beta0;
    SH(DIMj1+(1:MMAX))=-beta0*SV(2:MMAX+1,jj+1);
    for ii=1:MMAX
      cconj=cci-diag(alp(:,jj+1,ii+1));
      d=u(:,jj+2,ii+1)+u(:,jj+1,ii+2)-cconj*u(:,jj+1,ii+1)+u(:,jj+1,ii)+u(:,jj,ii+1);
      d(1)=-beta0*SS(jj+1,ii+1);
      d(MMAX+2)=beta0*SS1(jj+1,ii+1);
   
      SZ(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+3)=-beta0;
      SZ(DIMj+5*MMAX+(5*MMAX+6)*(ii-1)+(1:MMAX)+3)=ci-alp((1:MMAX)+1,jj+1,ii+1);
      SZ(DIMj+6*MMAX+(5*MMAX+6)*(ii-1)+4)=beta0;     
      SH(DIMj1+MMAX+(MMAX+2)*(ii-1)+(1:MMAX+2))=d(:);
    end  
   SZ(DIMj3-2*MMAX+(1:MMAX))=beta0;
   SH(DIMj2-MMAX+(1:MMAX))=beta0*SV1(2:MMAX+1,jj+1);
   end
  
   % Creating the sparse matrix and solving the system
   SR=bicg(sparse(SY,SX,SZ),SH,precision,NIT); 
   
   %Determining u
   for ii=1:MMAX
     u(2:MMAX+1,1,ii+1)=SR((1:MMAX)+MMAX*(ii-1));
     u(2:MMAX+1,MMAX+2,ii+1)=SR(DIM-MMAX^2+(1:MMAX)+MMAX*(ii-1));
   end
   for jj=1:MMAX
     u(2:MMAX+1,jj+1,1)=SR(MMAX^2+DIM1*(jj-1)+(1:MMAX));
     for ii=1:MMAX
       u(1:MMAX+2,jj+1,ii+1)=SR(MMAX^2+DIM1*(jj-1)+MMAX+(MMAX+2)*(ii-1)+(1:MMAX+2));
     end
     u(2:MMAX+1,jj+1,MMAX+2)=SR(MMAX^2+DIM1*jj-MMAX+(1:MMAX));
   end
   
   % Values in the corners
   u(1,1,1:MMAX+2)=(u(2,1,1:MMAX+2)+u(1,2,1:MMAX+2))/2; 
   u(1,MMAX+2,1:MMAX+2)=(u(2,MMAX+2,1:MMAX+2)+u(1,MMAX+1,1:MMAX+2))/2; 
   u(MMAX+2,1,1:MMAX+2)=(u(MMAX+1,1,1:MMAX+2)+u(MMAX+2,2,1:MMAX+2))/2; 
   u(MMAX+2,MMAX+2,1:MMAX+2)=(u(MMAX+1,MMAX+2,1:MMAX+2)+u(MMAX+2,MMAX+1,1:MMAX+2))/2;
   u(1:MMAX+2,1,1)=(u(1:MMAX+2,2,1)+u(1:MMAX+2,1,2))/2; 
   u(1:MMAX+2,1,MMAX+2)=(u(1:MMAX+2,2,MMAX+2)+u(1:MMAX+2,1,MMAX+1))/2; 
   u(1:MMAX+2,MMAX+2,1)=(u(1:MMAX+2,MMAX+1,1)+u(1:MMAX+2,MMAX+2,2))/2; 
   u(1:MMAX+2,MMAX+2,MMAX+2)=(u(1:MMAX+2,MMAX+1,MMAX+2)+u(1:MMAX+2,MMAX+2,MMAX+1))/2;
   u(1,1:MMAX+2,1)=(u(2,1:MMAX+2,1)+u(1,1:MMAX+2,2))/2; 
   u(1,1:MMAX+2,MMAX+2)=(u(2,1:MMAX+2,MMAX+2)+u(1,1:MMAX+2,MMAX+1))/2; 
   u(MMAX+2,1:MMAX+2,1)=(u(MMAX+1,1:MMAX+2,1)+u(MMAX+2,1:MMAX+2,2))/2; 
   u(MMAX+2,1:MMAX+2,MMAX+2)=(u(MMAX+1,1:MMAX+2,MMAX+2)+u(MMAX+2,1:MMAX+2,MMAX+1))/2;
   
   %Storing the boundary values
   ubottom(:,:,cntn+1)=u(:,:,2);
   utop(:,:,cntn+1)=u(:,:,MMAX+1);
   uleft(:,:,cntn+1)=u(2,:,:);
   uright(:,:,cntn+1)=u(MMAX+1,:,:);
   uforward(:,:,cntn+1)=u(:,2,:);
   uback(:,:,cntn+1)=u(:,MMAX+1,:);
   
   %Facets boundary conditions
   PLANE12(1:MMAX+2,cntn+1,1:cntn+1)=utop(1:MMAX+2,2,1:cntn+1); PLANE12(1:MMAX+2,1:cntn+1,cntn+1)=uforward(1:MMAX+2,MMAX+1,1:cntn+1);
   PLANE56(1:MMAX+2,cntn+1,1:cntn+1)=utop(1:MMAX+2,MMAX+1,1:cntn+1); PLANE56(1:MMAX+2,1:cntn+1,cntn+1)=uback(1:MMAX+2,MMAX+1,1:cntn+1);
   PLANE34(1:MMAX+2,cntn+1,1:cntn+1)=ubottom(1:MMAX+2,2,1:cntn+1); PLANE34(1:MMAX+2,1:cntn+1,cntn+1)=uforward(1:MMAX+2,2,1:cntn+1);
   PLANE78(1:MMAX+2,cntn+1,1:cntn+1)=ubottom(1:MMAX+2,MMAX+1,1:cntn+1); PLANE78(1:MMAX+2,1:cntn+1,cntn+1)=uback(1:MMAX+2,2,1:cntn+1);
   
   PLANE15(1:MMAX+2,cntn+1,1:cntn+1)=utop(2,1:MMAX+2,1:cntn+1); PLANE15(1:MMAX+2,1:cntn+1,cntn+1)=uleft(1:MMAX+2,MMAX+1,1:cntn+1);
   PLANE26(1:MMAX+2,cntn+1,1:cntn+1)=utop(MMAX+1,1:MMAX+2,1:cntn+1); PLANE26(1:MMAX+2,1:cntn+1,cntn+1)=uright(1:MMAX+2,MMAX+1,1:cntn+1);
   PLANE37(1:MMAX+2,cntn+1,1:cntn+1)=ubottom(2,1:MMAX+2,1:cntn+1); PLANE37(1:MMAX+2,1:cntn+1,cntn+1)=uleft(1:MMAX+2,2,1:cntn+1);
   PLANE48(1:MMAX+2,cntn+1,1:cntn+1)=ubottom(MMAX+1,1:MMAX+2,1:cntn+1); PLANE48(1:MMAX+2,1:cntn+1,cntn+1)=uright(1:MMAX+2,2,1:cntn+1);
   
   PLANE13(1:MMAX+2,cntn+1,1:cntn+1)=uforward(2,1:MMAX+2,1:cntn+1); PLANE13(1:MMAX+2,1:cntn+1,cntn+1)=uleft(2,1:MMAX+2,1:cntn+1);
   PLANE24(1:MMAX+2,cntn+1,1:cntn+1)=uforward(MMAX+1,1:MMAX+2,1:cntn+1); PLANE24(1:MMAX+2,1:cntn+1,cntn+1)=uright(2,1:MMAX+2,1:cntn+1);
   PLANE57(1:MMAX+2,cntn+1,1:cntn+1)=uback(2,1:MMAX+2,1:cntn+1); PLANE57(1:MMAX+2,1:cntn+1,cntn+1)=uleft(MMAX+1,1:MMAX+2,1:cntn+1);
   PLANE68(1:MMAX+2,cntn+1,1:cntn+1)=uback(MMAX+1,1:MMAX+2,1:cntn+1); PLANE68(1:MMAX+2,1:cntn+1,cntn+1)=uright(MMAX+1,1:MMAX+2,1:cntn+1);
   
   %Edge boubdary conditions
   EDGE1(cntn+1,1:cntn+1,1:cntn+1)=PLANE12(2,:,:); EDGE1(1:cntn+1,cntn+1,1:cntn+1)=PLANE15(2,:,:); EDGE1(1:cntn+1,1:cntn+1,cntn+1)=PLANE13(MMAX+1,:,:);
   EDGE2(cntn+1,1:cntn+1,1:cntn+1)=PLANE12(MMAX+1,:,:); EDGE2(1:cntn+1,cntn+1,1:cntn+1)=PLANE26(2,:,:); EDGE2(1:cntn+1,1:cntn+1,cntn+1)=PLANE24(MMAX+1,:,:);
   EDGE5(cntn+1,1:cntn+1,1:cntn+1)=PLANE56(2,:,:); EDGE5(1:cntn+1,cntn+1,1:cntn+1)=PLANE15(MMAX+1,:,:); EDGE5(1:cntn+1,1:cntn+1,cntn+1)=PLANE57(MMAX+1,:,:);
   EDGE6(cntn+1,1:cntn+1,1:cntn+1)=PLANE56(MMAX+1,:,:); EDGE6(1:cntn+1,cntn+1,1:cntn+1)=PLANE26(MMAX+1,:,:); EDGE6(1:cntn+1,1:cntn+1,cntn+1)=PLANE68(MMAX+1,:,:);
   
   EDGE3(cntn+1,1:cntn+1,1:cntn+1)=PLANE34(2,:,:); EDGE3(1:cntn+1,cntn+1,1:cntn+1)=PLANE37(2,:,:); EDGE3(1:cntn+1,1:cntn+1,cntn+1)=PLANE13(2,:,:);
   EDGE4(cntn+1,1:cntn+1,1:cntn+1)=PLANE34(MMAX+1,:,:); EDGE4(1:cntn+1,cntn+1,1:cntn+1)=PLANE48(2,:,:); EDGE4(1:cntn+1,1:cntn+1,cntn+1)=PLANE24(2,:,:);
   EDGE7(cntn+1,1:cntn+1,1:cntn+1)=PLANE78(2,:,:); EDGE7(1:cntn+1,cntn+1,1:cntn+1)=PLANE37(MMAX+1,:,:); EDGE7(1:cntn+1,1:cntn+1,cntn+1)=PLANE57(2,:,:);
   EDGE8(cntn+1,1:cntn+1,1:cntn+1)=PLANE78(MMAX+1,:,:); EDGE8(1:cntn+1,cntn+1,1:cntn+1)=PLANE48(MMAX+1,:,:); EDGE8(1:cntn+1,1:cntn+1,cntn+1)=PLANE68(2,:,:);  
   
   % Sparse matrices
   if cntn/sprsn-floor(cntn/sprsn)==0
      nuu=nuu+1;
      zplot(nuu)=z(cntn);
      uplot(1:muMAX,1:muMAX,1:muMAX,nuu)=u(sprsm*(1:muMAX),sprsm*(1:muMAX),sprsm*(1:muMAX));  
    end
   progress=round(cntn/NMAX*100)
end
toc
rplot=rplot-RMAX;

[chush,szplot]=size(zplot);
STRING = [sprintf('|u|^2:R2 =%f T =%f',RMAX,szplot*tau)];
step=[1, floor((szplot-1)/3)+1, szplot-floor(szplot/3), szplot]; 
lbl={'(a)','(b)','(c)','(d)'};

% Movie making

% for s=1:szplot
%   pcolor(rplot*1e-3,rplot*1e-3,abs(squeeze(uplot(:,:,floor((muMAX+1)/2),s)).^2))
%   colormap('jet')
%   shading interp
%   colorbar
%   caxis ([0 0.1])
%   xlabel('x, \mum')
%   ylabel('y, \mum')
%   title(STRING)
%   F(s)=getframe;
% end
% movie(F,1,1);

% avi = avifile('E:\TEXDOCS\Boundary condition\XR_flat_cyl.avi','fps',3,'compression','Indeo3');
% for s=1:szplot
%   pcolor(rplot*1e-3,rplot*1e-3,abs(uplot(floor((muMAX+1)/2),:,:,s).^2))
%   colormap('jet')
%   shading interp
%   colorbar
%   caxis ([0 2])
%   xlabel('x, \mum')
%   ylabel('y, \mum')
%   title(STRING)
%   frame=getframe;
%   avi = addframe(avi,frame);
% end
% avi = close(avi);

% Projections

% pp=zeros(muMAX,szplot);
% pp(:,:)=uplot(:,floor((muMAX+1)/2),floor((muMAX+1)/2),:);
% figure
% pcolor(zplot*1e-6,rplot*1e-3,abs(pp.^2))
% colormap('jet')
% shading interp
% colorbar
% caxis ([0 0.3])
% xlabel('z, mm')
% ylabel('y, \mum')
% title(STRING)

figure('PaperOrientation','portrait','PaperType', 'A4')
for s=1:4
  %tt=strcat(int2str(zplot(step*(s-1)+1)));
  subplot(2,2,s)
  set(gca,'FontSize',16)
  pcolor(rplot,rplot,log10(squeeze(abs(uplot(:,:,floor((muMAX+1)/2),step(s)).^2)))) 
  colormap jet
  shading interp
  axis equal
  caxis ([-6 -2])
  colorbar('FontSize',16)
  ylabel('z (arb. units)','FontSize',18)
  xlabel('x (arb. units)','FontSize',18)
  title(lbl(s),'FontSize',20)
end

% figure('PaperOrientation','portrait','PaperType','A4')
% isos=0.01;
% for s=1:4
%   tt=strcat(int2str(zplot(step*(s-1)+1)*1e-3),' a.u.');
%   subplot(2,2,s)
%   p1 = patch(isosurface(rplot*1e-3,rplot*1e-3,rplot*1e-3,abs(uplot(:,:,:,step(s))).^2,isos),'FaceColor','red',...
%     'EdgeColor','none');
%   p2 = patch(isocaps(rplot*1e-3,rplot*1e-3,rplot*1e-3,abs(uplot(:,:,:,step(s))).^2,isos),'FaceColor','interp',...
%     'EdgeColor','none');
%   view(3); axis tight;
%   axis equal
%   axis([-RMAX*1e-3 RMAX*1e-3 -RMAX*1e-3 RMAX*1e-3 -RMAX*1e-3 RMAX*1e-3])
%   colormap jet
%   xlabel('x, a.u.')
%   ylabel('y, a.u.')
%   zlabel('z, a.u.')
%   camlight left; lighting gouraud
%   isonormals(rplot*1e-3,rplot*1e-3,rplot*1e-3,abs(uplot(:,:,:,step(s))).^2,p1)
%   title(tt,'FontSize',16)
% end

% % 3D movie making
% isos=1.4;
% for st=1:10
%   h=figure;
%   p1 = patch(isosurface(rplot*1e-3,rplot*1e-3,rplot*1e-3,abs(uplot(:,:,:,st)).^2,isos),'FaceColor','red',...
%     'EdgeColor','none');
%   p2 = patch(isocaps(rplot*1e-3,rplot*1e-3,rplot*1e-3,abs(uplot(:,:,:,st)).^2,isos),'FaceColor','interp',...
%     'EdgeColor','none');
%   view(3); axis tight;
%   axis equal
%   axis([-RMAX*1e-3 RMAX*1e-3 -RMAX*1e-3 RMAX*1e-3 -RMAX*1e-3 RMAX*1e-3])
%   colormap hsv
%   xlabel('x, \mum')
%   ylabel('y, \mum')
%   zlabel('z, \mum')
%   camlight left; camlight; lighting gouraud
%   isonormals(rplot*1e-3,rplot*1e-3,rplot*1e-3,abs(uplot(:,:,:,st)).^2,p1)
%   F(step)=getframe;
%   close(h);
% end
% movie(F,1,1);

% A big diagonal cross-section
% figure('PaperOrientation','portrait','PaperType','A4')
% for s=1:4
%    for t=1:muMAX
%      uuplot(t,s)=2*log(abs(uplot(t,t,t,step(s)))); 
%    end
%   tt=strcat(int2str(zplot(step(s))),' a.u.');
%   subplot(2,2,s)
%   plot(sqrt(3)*rplot,uuplot(:,s))
%   xlabel('R, a.u.')
%   ylabel('I, a.u.')
%   title(tt)
%   axis ([-sqrt(3)*RMAX sqrt(3)*RMAX -15 -2])
% end

%A small diagonal cross-section
% figure('PaperOrientation','portrait','PaperType','A4')
% for s=1:4
%   for t=1:muMAX
%      uuplot1(t,s)=2*log(abs(uplot(t,floor((muMAX+1)/2),t,step*(s-1)+1))); 
%   end
%   tt=strcat(int2str(zplot(step(s))),' a.u.');
%   subplot(2,2,s)
%   plot(sqrt(2)*rplot,uuplot1(:,s))
%   xlabel('\rho, a.u.')
%   ylabel('I, a.u.')
%   title(tt)
%   axis ([-sqrt(2)*RMAX sqrt(2)*RMAX -15 -2])
% end

% IT(1:muMAX,1)=sqrt(3)*rplot;
% IT(1:muMAX,2)=uuplot(1:muMAX,1);
% IT(1:muMAX,3)=uuplot(1:muMAX,2);
% IT(1:muMAX,4)=uuplot(1:muMAX,3);
% IT(1:muMAX,5)=uuplot(1:muMAX,4);
% ITT(1:muMAX,1)=sqrt(2)*rplot;
% ITT(1:muMAX,2)=uuplot1(1:muMAX,1);
% ITT(1:muMAX,3)=uuplot1(1:muMAX,2);
% ITT(1:muMAX,4)=uuplot1(1:muMAX,3);
% ITT(1:muMAX,5)=uuplot1(1:muMAX,4);
% 
% save 4D_gauss_semisphere_3corner_2.5_20.dat IT -ascii
% save 4D_gauss_semisphere_2corner_2.5_20.dat ITT -ascii