function R=basicgrowthmodel(k_init_ratio)

%% Advanced Macroeconomics 4 2022, NM
%%Solve the basic neoclassical growth model with CRRA utility
%%and Cobb-Douglas production function

%%Save this file somewhere where your Matlab finds them (you may use the File/Set path command in Matlab). 
%%Then in the command window write e.g. R=basicgrowthmodel(0.5) and press Enter. The output structure R contains the optimal allocation. 

%INPUT:
%      k_init_ratio: the initial capital stock relative to the steady state capital stock
%OUTPUT:
%      R.k: capital stock path
%      R.c: consumption path
%      R.n: labor supply path
%      also plots the paths
%OTHER FUNCTIONS:
%      util.m ss_equs.m trans_equs.m

%create a structure Params that collects all parameter values
Params.beta=0.9;%discount factor
Params.c_share=0.5;%consumption share in the utility function
Params.delta=0.1;%capital depreciation rate
Params.alpha=0.333;%capital share
Params.sigma=2;%CRRA parameter for the utility function, 1 corresponds to log utility
Params.T=60;%number of transition periods, should be big enough so that the paths converge.

%guess for the steady state
%%add one for the discount rate

guess=[0.5;0.5;0.5;0.5];%guess for the ss capital stock, consumption and labor supply
opts=optimset();%need to define this structure to add additional input arguments
%to ss_equs and trans_equs used by fsolve see MATLAB help for optimset and fsolve. 

%solve for the steady state
[ss, errors]=fsolve(@ss_equs,guess,opts,Params);

%check that ss equations are solved accurately
assert(max(abs(errors))<1e-5,'ss not solved')

%%steady state allocation
ss_k=ss(1);
ss_c=ss(2);
ss_n=ss(3);

disp('steady state capital, consumption and labour:')
[ss_k ss_c ss_n]

%determine the initial capital stock
Params.k_init=k_init_ratio*ss_k;

%guess for the transition. (a rather stupid initial guess, if you have
%problems, try to construct a more reasonable guess that converges slowly to the ss.)
guess=[ss_k*ones(Params.T-1,1);ss_n*ones(Params.T,1);ss_c*ones(Params.T,1)];

%solve for the transition
[tr, errors]=fsolve(@trans_equs,guess,opts,Params);

assert(max(abs(errors))<1e-4,'transition not solved')

T=Params.T;
%reshape the result (a long vector) into capital, labor and consumption
%paths
R.k=[Params.k_init;tr(1:T-1,1)]; %inital capital stock + capital stocks in periods 2,3,4,...
R.n=tr(T:2*T-1,1); %labor supply in periods 1,2,3,...
R.c=tr(2*T:3*T-1,1); %consumption in periods 1,2,3...

%make some figures
figure(1), plot(R.k,'-*r')
ylabel('Capital stock')
xlabel('Period')

figure(2), plot(R.n,'-*r')
ylabel('Labour')
xlabel('Period')

figure(3), plot(R.c,'-*r')
ylabel('Consumption')
xlabel('Period')

end

function z=ss_equs(x,Params)

k=x(1);
c=x(2);
n=x(3);
%extend x with the discount factor
disc = x(4);

f=k^Params.alpha*n^(1-Params.alpha);
fdk=Params.alpha*k^(Params.alpha-1)*n^(1-Params.alpha);
fdn=(1-Params.alpha)*k^Params.alpha*n^(-Params.alpha);

%add equation for the calibration


[~, udc, udl]=util(c,n,Params);

%steady state equations
z=zeros(3,1);

disc = 

z(1)=Params.beta*(1+fdk-Params.delta)-1;
z(2)=udc*fdn-udl;
z(3)=c+Params.delta*k-f;


end

function [u, duc, dul]=util(c,n,Params)

if Params.sigma==1
    u=Params.c_share*log(c)+(1-Params.c_share)*log(1-n);
    duc=Params.c_share./c;
    dul=(1-Params.c_share)./(1-n);
else
    aux1=c.^(Params.c_share).*(1-n).^(1-Params.c_share);
    u=aux1.^(1-Params.sigma)/(1-Params.sigma);
    aux2=aux1.^(-Params.sigma);
    duc=aux2.*Params.c_share.*c.^(Params.c_share-1).*(1-n).^(1-Params.c_share);
    dul=aux2.*(1-Params.c_share).*c.^Params.c_share.*(1-n).^(-Params.c_share);
end
end

function z=trans_equs(x,Params)

T=Params.T;
alpha=Params.alpha;
k=zeros(T,1);

k(2:T)=x(1:T-1,1);%  guess for the capital stock in periods 2,3,4,...,T
n=x(T:2*T-1,1);% guess for the labor supply in periods 1,2,3,...,T
c=x(2*T:3*T-1,1);% guess for consumption in periods 1,2,3,...,T
k(1)=Params.k_init;

%output and its derivatives
f=k.^alpha.*n.^(1-alpha);
fdk=alpha*k.^(alpha-1).*n.^(1-alpha);
fdn=(1-alpha)*k.^alpha.*n.^(-alpha);

aggresource=zeros(T-1,1);
laborfoc=zeros(T,1);
euler=zeros(T,1);

[~, udc, udl]=util(c,n,Params);

% specify a system of 3*T -1 equations that need to be solved
for t=1:T-1
    laborfoc(t)=udc(t)*fdn(t)-udl(t);
    euler(t)=udc(t)-Params.beta*udc(t+1)*(1+fdk(t+1)-Params.delta);
    aggresource(t)=c(t)+k(t+1)-f(t)-(1-Params.delta)*k(t);
end

laborfoc(T)=udc(T)*fdn(T)-udl(T);
euler(T)=udc(T)-Params.beta*udc(T)*(1+fdk(T)-Params.delta);

%stack the errors
z=[aggresource;laborfoc;euler];

end

