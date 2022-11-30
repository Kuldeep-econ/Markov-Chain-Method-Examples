%% Example 1: Growth Model policy function
function [c]=update_cons(V)

ga = 2;     %Utility parameter
a = 0.3;    %Production parameter alpha   
A = 1;      %TFP 
rho = 0.05; %discount rate
delta=0.15; %depreciation rate
%si=0.01;    %drift term

dt=10^-3;   % Delta_T 


kss = (a*A/(rho+delta))^(1/(1-a));    %Steady state capital with no shocks

kmin = 0.2*kss;
kmax = 2*kss;
dk =0.05;
k =(kmin:dk:kmax)';    % Capital Grid
I=length(k);

dvkf = zeros(I,1); %dv/dk term
dvkb = zeros(I,1);


C_f=zeros(I,1);
U_f=zeros(I,1);
U_f(I,1)=-inf;  % Adjustment so that we never choose C_f at grid end

C_b=zeros(I,1);
U_b=zeros(I,1);
U_b(1,1)=-inf;  % Adjustment so that we never choose C_b at grid start


C_s=A*k.^a-delta*k;
U_s=C_s.^(1-ga)/(1-ga);


% Notice dk term, you can also bring it later on if you want

 dvkf(1:I-1,1) = (V(2:I,1)-V(1:I-1,1))/dk;  
 dvkb(2:I,1) = (V(2:I,1)-V(1:I-1,1))/dk;
 dvkb(1,:)=10^10;  %adjustments you need at boundary, check slides
 dvkf(I,:)=-10^10;
 
    C_f(1:I-1,1)= (exp(-rho*dt).*dvkf(1:I-1,1)).^(-1/ga);
    U_f(1:I-1,1)= C_f(1:I-1,1).^(1-ga)/(1-ga); 
    
    
    C_b(2:I,1)= (exp(-rho*dt).*dvkb(2:I,1)).^(-1/ga);
    U_b(2:I,1)= C_b(2:I,1).^(1-ga)/(1-ga); 
     
    % Getting ready for construction of F,B,S in my slides
   Y_const= A.*k.^a-delta*k; % Common term in F,B
   
   
   % For F term
   P_k_f=max(Y_const-C_f,0);
   P_k_b=max(-(Y_const-C_f),0);
   
   % For B term
   M_k_f=max(Y_const-C_b,0);
   M_k_b=max(-(Y_const-C_b),0);
   
   
   
   % Ignoring common terms in construction below 
 F_term=dt*exp(-rho*dt)*(P_k_f.*dvkf -P_k_b.*dvkb);
 F =dt* U_f +F_term ;% F_I is -inf because U_f(I)=-inf 
 
 B_term=dt*exp(-rho*dt)*(M_k_f.*dvkf -M_k_b.*dvkb);
 B= dt*U_b +B_term;% B_1 is -inf because U_b(1)=-inf 
 
 
 S=dt*U_s ;
    
    MM=cat(3,F,B,S);   %  %combines F,B,S into  3 matrices on 3 pages
    
   [~,pos]=max(MM,[],3); % Tells you position of the maximimzer
   
   c=C_f.*(pos==1)+ C_b.*(pos==2)+C_s.*(pos==3);
   
    
  

end
