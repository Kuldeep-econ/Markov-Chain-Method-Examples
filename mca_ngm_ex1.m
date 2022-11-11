%% Markov Chain Approximation Method
% Example 1 : Growth Model

%Parameters


ga = 2;     %Utility parameter
a = 0.3;    %Production parameter alpha   
A = 1;      %TFP 
rho = 0.05; %discount rate
delta=0.15; %depreciation rate
si=0.01;    %diffusion term

dt=10^-3;   % Delta_T 


kss = (a*A/(rho+delta))^(1/(1-a));    %Steady state capital with no shocks

kmin = 0.2*kss;
kmax = 2*kss;
dk =0.05;
k =(kmin:dk:kmax)';    % Capital Grid
I=length(k);

maxit= 100;
crit = 1e-8;   %Stopping Criteria


% Initial Guess
c0 = A.*k.^a/2;    
c=c0;

Pk_plus=zeros(I,1);     %P(k+delta_k)
Pk_minus=zeros(I,1);     %P(k-delta_k) 



%% Main Loop 
tic;
for j=1:maxit

C=c;           % Step 1 
 
     U=C.^(1-ga)/(1-ga);
     ss= A.*k.^a-delta*k - C  ;
     
   Pk_plus(1:I-1,1)=dt*(max(  ss(1:I-1,1) , 0))/dk + dt*si^2/(2*dk^2);  %P(k+delta_k)
   Pk_minus(2:I,1)=dt*(max(-ss(2:I,1), 0))/dk  + dt*si^2/(2*dk^2);      %P(k-delta_k)  
   
  % Check if Probabilities are between (0,1)
     if max(Pk_plus)>1||max(Pk_minus)>1
          disp('Probability >1')
          break
      end
%              
  
  % A construction
  P_middle=-(Pk_plus+Pk_minus);
  
    Ab=spdiags([Pk_minus(2:I,1);0],-1,I,I);   % Check out spdiags, -1 here targets the diagnol we want
    Af=spdiags([0;Pk_plus(1:I-1,1)],1,I,I);
    Ac=spdiags(P_middle,0,I,I);
    
    Amain=Af+Ab+Ac;    %A matrix in my slides
    
    BB=(speye(I)-exp(-rho*dt)*(speye(I)+Amain))/dt;
    V=BB\U;             % Step 2
    
    
    C=update_cons(V);  % Step 3, (check update_cons file)
    
    Cchange = C - c;  
    c = C;
   
    dist(j)= max(abs(Cchange),[],'all');
    if dist(j)<crit    % Step 4
        disp('Policy Converged, Iteration = ')
        disp(j)
       break
    end 
end
toc;

%% Graphs




set(gca,'FontSize',14)
subplot(2,2,1)
plot(k,V,'LineWidth',1.5)
xlabel('k')
ylabel('Value')

subplot(2,2,2)
plot(k,c,'LineWidth',1.5)
xlabel('k')
ylabel('Consumption')
axis tight


subplot(2,2,3)
plot(k,ss,'LineWidth',1.5)
xlabel('k')
ylabel('Savings')
axis tight

