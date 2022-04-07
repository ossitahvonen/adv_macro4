
Params.delta=0.1;
Params.alpha=0.4;
Params.sigma=2;
Params.beta=0.96;
Params.c_share=0.4;
opts=optimset();%needed to add addional arguments when using fsolve or fzero

%The initial tax structure
Params.tk1=0.3;Params.tc1=0.2;Params.tn1=0.3;

%The new tax structure, tn will be adjusted to satisfy the (intertemporal)
%budget constraint. tk will equal tk1 for period 1 and tk2 afterwards tn
%and tc will be set equal to tn2 and tc2 from period 1 onwards.
Params.tk2=0.0;Params.tc2=0.2;

Params.T=100;%number of transition periods, should be big enough

%Guess for the steady state k,c,n and g
%We assume balanced budget in the SS, so g should be equal to tax revenue
guess=[1.5;0.5;0.5;0.2];

%set s

s=10


ss=guess
Params.k_init = 1;
ss_k=ss(1);ss_c=ss(2);ss_n=ss(3);g=ss(4);

x=[ss_k*ones(Params.T-1,1);ss_c*ones(Params.T,1);ss_n*ones(Params.T,1);0.2];

T=Params.T;alpha=Params.alpha;
k=x(1:T-1,1);
c=x(T:2*T-1,1);
n=x(2*T:3*T-1,1);
tn=x(3*T);%fixed tn that solves the gov. budget constraint

k=[Params.k_init;k];

%%output and its derivatives
f=k.^alpha.*n.^(1-alpha);
fdk=alpha*k.^(alpha-1).*n.^(1-alpha);
fdn=(1-alpha)*k.^alpha.*n.^(-alpha);

r=fdk-Params.delta;
w=fdn;
tk=Params.tk2*ones(T,1);


tk(1:s)=Params.tk1;
tn=tn*ones(T,1);
tc=Params.tc2*ones(T,1);