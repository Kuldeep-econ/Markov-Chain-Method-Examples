%% Function will give all the parameters 
function output =parameter()

%Utility function parameters
parameters.si = 2; %CRRA utility parameter sigma
parameters.rho = 0.05; %discount rate

parameters.r =0.03;
parameters.w=1;



dt=1e-2/2;%;1e-2;
parameters.dt=dt;


%Income shock
%ORNSTEIN-UHLENBECK PROCESS dlog(z) = -the*log(z)dt + sig2*dW
%STATIONARY DISTRIBUTION IS log(z) ~ N(0,Var) WHERE Var = sig2/(2*the)
Var = 0.07;
zmean = exp(Var/2); %MEAN OF LOG-NORMAL DISTRIBUTION N(0,Var)
Corr = 0.9;
the= -log(Corr);
parameters.the =the;
parameters.sig2 = 2*the*Var;
zmin = zmean*0.5;
zmax = zmean*1.5;

%a,h Grid related
da=0.05;
dz=0.05;
parameters.da=da;
parameters.dz=dz;

parameters.zmin =zmin; 
parameters.zmax=zmax;

amin=0;
amax=20;
parameters.amin =amin;
parameters.amax =amax;





output=parameters;
end