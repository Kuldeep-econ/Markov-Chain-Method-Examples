function C=policy(V,parameters)

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




[zz,aa]=meshgrid(z,a);


dvaf = zeros(I,J);
dvab = zeros(I,J);



C_f=zeros(I,J);
S_f=zeros(I,J);
U_f=zeros(I,J);
U_f(I,:)=-inf  ;% Adjustment so that we never choose F term at grid end

C_b=zeros(I,J);
S_b=zeros(I,J);
U_b=zeros(I,J);
U_b(1,:)=-inf  ; % Adjustment so that we never choose B term at grid start


C_s=w*zz+ r.*aa;
U_s=C_s.^(1-si)/(1-si) ;

%% Main Loop


dvaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
dvab(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
dvab(1,:)=10^10;  % Check out notes
dvaf(I,:)=-10^10; %  

  
 
   C_f(1:I-1,:)= (exp(-rho*dt).*dvaf(1:I-1,:)).^(-1/si);
   S_f(1:I-1,:)=w*zz(1:I-1,:)+ r*aa(1:I-1,:)    - C_f(1:I-1,:)   ;
    U_f(1:I-1,:)= C_f(1:I-1,:).^(1-si)/(1-si);
    
    C_b(2:I,:)= (exp(-rho*dt).*dvab(2:I,:)).^(-1/si);
    S_b(2:I,:)=w*zz(2:I,:) +  r*aa(2:I,:)- C_b(2:I,:) ;
     U_b(2:I,:)= C_b(2:I,:).^(1-si)/(1-si);
    
     
   Pa_p_f=(max( S_f, 0)) ;
   Pa_p_b=(max( -S_f, 0)) ;
   
   Pa_m_b=(max(-S_b, 0)) ;  
   Pa_m_f=(max(S_b, 0)) ;  
   
  F = U_f + exp(-rho*dt)*(Pa_p_f.*dvaf-Pa_p_b.*dvab); 
  B= U_b + exp(-rho*dt)*(Pa_m_f.*dvaf-Pa_m_b.*dvab);
  S=U_s; 
    
  
 
    MM=cat(3,F,B,S);  %combines F,B,S into 3 d array or 3 matrices on 3 pages
    [~,index]=max(MM,[],3); %gives max of each element from all pages, 3 means z-axis
   
   
  
    C=(index==1).*C_f+ (index==2).*C_b+(index==3).*C_s ;

end