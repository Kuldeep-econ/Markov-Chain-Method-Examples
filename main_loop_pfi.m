%% Policy Function Iteration (Matlab code)
%% You need parameter file and policy file to run this

parameters=parameter();
si=parameters.si ; %CRRA utility parameter sigma
rho=parameters.rho ; %discount rate
r=parameters.r ; %interest rate
w=parameters.w ;

% time related
dt=parameters.dt;


% Grid related
da=parameters.da;
amin=parameters.amin ;
amax=parameters.amax ;
a =(amin:da:amax)';
I=length(a);  

dz=parameters.dz;
zmin=parameters.zmin ;
zmax=parameters.zmax ;
z =(zmin:dz:zmax);
J=length(z);  

the=parameters.the ;
sig2=parameters.sig2; 



[zz,aa]=meshgrid(z,a);

Pa_plus=zeros(I,J);     
Pa_minus=zeros(I,J);  



dz2 = dz^2;
mu = (-the*log(zz) + sig2/2).*zz; %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*zz.^2; %VARIANCE (FROM ITO'S LEMMA)


Pz_fixed=s2;
Pz_plus= dt*(max(mu,0))/dz + dt*Pz_fixed./(2*dz2);
Pz_plus(:,J)=0;
Pz_minus=dt*(max(-mu,0))/dz + dt*Pz_fixed./(2*dz2);
Pz_minus(:,1)=0;
if max(Pz_plus,[],'all')>1|| max(Pz_minus,[],'all')>1
    disp('Probablity z > 1')
end

Pz_up_waste=sparse(I,1);   %needed in construction of A matrix, check out spdiags in matlab


v0 = (w*zz + r.*aa).^(1-si)/(1-si)/rho;%to get an initial reasonable guess, i use this value function
c0=policy(v0,parameters); %inital consumption guess
c=c0;


maxit=1e+2;
crit=1e-8;


%% main loop
%tic;
for i=1:maxit
    C=c;
    
    
U=C.^(1-si)/(1-si);
u_stacked=reshape(U,I*J,1);


ss=w*zz+ r.*aa- C  ;
     
prob= dt/da;
Pa_plus(1:I-1,:)=prob*max(ss(1:I-1,:) , 0);  %P(a+delta_a,)
Pa_minus(2:I,:)= prob*max(-ss(2:I,:), 0) ;    %P(a-delta_a,)
   
     if max(abs(ss)*dt/da,[],'all')>1
         disp('Probability >1')
     end
                
 P_middle= -(Pa_plus+Pa_minus+Pz_plus+Pz_minus) ; %middle diagonal
   
 
 Ab=spdiags(Pa_minus(2:I*J)',-1,I*J,I*J);
 Af=spdiags([0;Pa_plus(1:I*J-1)'],1,I*J,I*J);%0 will be ignored because of the way spdiags is defined
 Ac=spdiags(reshape(P_middle,I*J,1),0,I*J,I*J);
 Af_z=spdiags([Pz_up_waste;Pz_plus(1:I*J-I)'],I,I*J,I*J); 
 Ab_z=spdiags(Pz_minus(I+1:I*J)',-I,I*J,I*J);
 
 Amain=Ab+Af+Ac+Af_z+Ab_z;  %A matrix in the slides
 
    BB=(speye(I*J)-exp(-rho*dt)*(speye(I*J)+Amain))/dt;
    
    Vnew=BB\u_stacked;   
    v=reshape(Vnew,I,J);
    
    Cnew=policy(v,parameters);
    
    
    dist=max(abs(Cnew-C),[],'all');
    
    c=Cnew;
    
    if dist<crit
        break
    end
    
end
%toc;
