function R=taxreform
%% Advanced Macroeconomics 4 2022/NM
%% Neoclassical growth model with government bonds and flat-rate taxes 
%% Solve for the dynamics starting from a steady state and following a tax reform (constant tax rates). Present value budget constraint for the government. 

%OUTPUT:
%       Params: includes calibrated beta, g and gamma (c_share).
%       R.k: capital stock path
%       R.c: consumption path
%       R.n: labor supply path
%       also plots the paths

%OTHER FUNCTIONS:
%      util.m ss_equs.m trans_equs.m welfare_equ.m 

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

%Solve for the initial steady state
[ss,fval]=fsolve(@ss_equs,guess,opts,Params);
assert(max(abs(fval))<1e-5,"steady state not solved accurately")

%%the steady state allocation
ss_k=ss(1);ss_c=ss(2);ss_n=ss(3);g=ss(4);

%From now on we treat g and k1 as parameters.
Params.g=g;Params.k_init=ss_k;

%Add the initial steady state and parameters into the output structure
R.Params=Params;
R.ss_k=ss_k;R.ss_c=ss_c;R.ss_n=ss_n;

% Create a guess for the transition and fixed tn that solves the gov. budget constraint
guess=[ss_k*ones(Params.T-1,1);ss_c*ones(Params.T,1);ss_n*ones(Params.T,1);0.2];
[tr, fval]=fsolve(@trans_equs,guess,opts,Params);
assert(max(abs(fval))<1e-5,"transition not solved accurately")

%%Collect k,c, and n into R 
%%also add the steady state levels, this helps to check that T is big enough
R.k=[Params.k_init;tr(1:Params.T-1)];
R.c=tr(Params.T:2*Params.T-1);
R.n=tr(2*Params.T:3*Params.T-1);
R.tn=tr(3*Params.T);

[~,R.b,R.D,R.hh_budget]=trans_equs(tr,Params);%put the solution into trans_equs again to compute gov. bonds, deficits and to check the aggregate resource constraint

%create pictures
hold off
subplot(3,1,1),plot([ss_k;R.k],'-*');%period 1 in the figure corresponds to the initial steady state (period 0)
title('capital')
subplot(3,1,2),plot([ss_c;R.c],'-*');
title('consumption')
subplot(3,1,3),plot([ss_n;R.n],'-*');
title('labor')
xlabel('Period')

%finally, solve for the welfare effects (equivalent consumption compensation);
x_guess=0;
T=Params.T;
%steady state effect
R.welfare_eff_ss=fzero(@welfare_equ,x_guess,opts,R.ss_c,R.ss_n,R.c(T),R.n(T),Params);
%including transition
R.welfare_eff_tr=fzero(@welfare_equ,x_guess,opts,R.ss_c*ones(T,1),R.ss_n*ones(T,1),R.c,R.n,Params);
end

function z=ss_equs(x,Params)
z=zeros(4,1);

k=x(1);c=x(2);n=x(3);g=x(4);
tk=Params.tk1;tn=Params.tn1;tc=Params.tc1;
beta=Params.beta;

f=k^Params.alpha*n^(1-Params.alpha);
fdk=Params.alpha*k^(Params.alpha-1)*n^(1-Params.alpha);
fdn=(1-Params.alpha)*k^Params.alpha*n^(-Params.alpha);
r=fdk-Params.delta;
w=fdn;
[~, udc, udl]=util(c,n,Params);

z(1)=beta*(1+(1-tk)*r)-1;
z(2)=udc*w*(1-tn)/(1+tc)-udl;
z(3)=c*(1+tc)+k-(1+(1-tk)*r)*k-(1-tn)*w*n;
z(4)=f-c-g-Params.delta*k;

%deficit=g-tc*c-tn*w*n-tk*r*k 

end

function [z,b,D,hb]=trans_equs(x,Params)
%here we have the foc's and the agggregate resource constraint for the
%whole transition. 

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

tk(1)=Params.tk1;
tn=tn*ones(T,1);
tc=Params.tc2*ones(T,1);

%make sure that utility and its derivates are well defined.
%this can cause problems in some extreme cases! need to check that the
%solution is ok. 
c=max(c,1e-3);k=max(k,1e-3);n=max(n,1e-3);n=min(n,1-1e-3);

%get the derivatives of u, note that udn<0 (the deriv of u w.r.t to labor not leisure)
[~, udc, udl]=util(c,n,Params);

%preallocate
euler=zeros(T,1);%drop the last euler, because k(1) is given.
aggresource=zeros(T-1,1);

%euler equations
for t=1:T-1
    euler(t)=udc(t)/(1+tc(t))-Params.beta*(1+(1-tk(t+1))*r(t+1))*udc(t+1)/(1+tc(t+1));
end
    euler(T)=udc(T)/(1+tc(T))-Params.beta*(1+(1-tk(T))*r(T))*udc(T)/(1+tc(T));

%aggregate resource equations   
for t=1:T-1
aggresource(t)=c(t)+k(t+1)+Params.g-f(t)-(1-Params.delta)*k(t);
end


%intertemporal prices: p(t)
p=ones(T,1);

%get deficit (solve only after using fsolve)
D=zeros(T,1);

for t=1:T-1
    p(t+1)=p(t)/(1+(1-tk(t+1))*r(t+1));
end

%%government budget constraint
TaxRev=0;%present value of tax revenue
G=0;%present value of gov. spending
for t=1:T-1
   Tax=(tc(t)*c(t)+tn(t)*w(t)*n(t)+tk(t)*r(t)*k(t));
   TaxRev=TaxRev+p(t)*Tax;
   G=G+Params.g*p(t);
   D(t)=Params.g-Tax;
end
D(T)=D(T-1);


%% add new steady state...
R_ss=1+(1-tk(T))*r(T);
TaxRev=TaxRev+Tax*p(T)/(1-1/R_ss);%%The present value of tax revenues starting from the new steady state.
G=G+Params.g*p(T)/(1-1/R_ss);%The present value of gov. spending starting from the new steady state
BudConstr=TaxRev-G;

%labor supply focs, directly in vector format
laborfoc=udc.*w.*(1-tn)-udl.*(1+tc);

%collect the equations
z=[euler;laborfoc;aggresource;BudConstr];

if nargout==4
b=zeros(T,1);
hb=0;
for t=1:T-1
    b(t+1)=(1+(1-tk(t)))*r(t)*(D(t)+b(t));
    hb=hb+p(t)*(-(1+tc(t))*c(t)-k(t+1)+(1+(1-tk(t))*r(t))*k(t)+(1-tn(t))*w(t)*n(t));
end
hb=hb+p(T)/(1-1/R_ss)*(-(1+tc(T))*c(T)-k(T)+(1+(1-tk(T))*r(T))*k(T)+(1-tn(T))*w(T)*n(T));
end
end

function [u, duc, dul]=util(c,n,Params)
c_share=Params.c_share;
if(Params.sigma==1)
u=c_share*log(c)+(1-c_share)*log(1-n);
duc=c_share./c;
dul=(1-c_share)./(1-n);
else
aux1=c.^(c_share).*(1-n).^(1-c_share);
u=aux1.^(1-Params.sigma)/(1-Params.sigma);
aux2=aux1.^(-Params.sigma);
duc=aux2.*c_share.*c.^(c_share-1).*(1-n).^(1-c_share);
dul=aux2.*(1-c_share).*c.^c_share.*(1-n).^(-c_share);
end
end

function z=welfare_equ(x,c1,n1,c2,n2,Params)
U1=0;U2=0;T=length(c1);
for t=1:T
U1=U1+Params.beta^(t-1)*util((1+x)*c1(t),n1(t),Params);
U2=U2+Params.beta^(t-1)*util(c2(t),n2(t),Params);
end
U1=U1+Params.beta^(T)*util((1+x)*c1(T),n1(T),Params)/(1-Params.beta);
U2=U2+Params.beta^(T)*util(c2(T),n2(T),Params)/(1-Params.beta);
z=U1-U2;
end


