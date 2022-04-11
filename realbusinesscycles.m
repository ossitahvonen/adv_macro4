%%Advanced Macroeconomics 4 2022, NM
%Programs to solve a basic stochastic growth model with value funcion iteration.
%%Save this file somewhere where your Matlab finds them (you may use the File/Set path command in Matlab). 
%%Then in the command window write e.g. R=realbusinesscycles() and press
%%Enter. The output structure R contains the solution to the Bellman
%%equation and simulated paths.

function out=realbusinesscycles()
Params.N=100; %number of gridpoints for the asset state
Params.maxIter=100; %max number of value function iterations
Params.nSimulatedPeriods=5000;
Params.interptype='lin'; %Interpolation method (see interp1): linear is fast but requires more gridpoints
Params.k_min=1;%%min asset holdings, check that these bounds (min and max) are never binding
Params.k_max=3;%%max asset holdings
Params.sigma=2.0; %%risk-aversion
Params.beta=0.9;%%discount factor
Params.alpha=0.333;%capital share
Params.delta=0.1;%depreciation rate
Params.s=zeros(2,1);%%NOTE: for simplicity, these programs work only with two shocks. If you add shocks, you need to modify function value and the simulation part below a bit.
Params.z(1)=-0.05;%the first (bad) productivity shock (log)
Params.z(2)=0.05;%the second (good) productivity shock (log)
Params.P=[0.7 0.3; 0.3 0.7];%%the transition probability matrix for the income shock: P(i,j) is the prob that next income is j given that current shocks is i.
Params.n_shocks=length(Params.z);
assert(Params.n_shocks==size(Params.P,1));%just make sure P is of the correct size
Params.k_grid=linspace(Params.k_min,Params.k_max,Params.N); %linearily spaced gridpoints for capital
ValueFunction=zeros(Params.N,Params.n_shocks);
ValueFunction(:,1)=util(output(Params.k_grid,1,Params),Params)/(1-Params.beta);%our very first guess for the value function
ValueFunction(:,2)=util(output(Params.k_grid,1,Params),Params)/(1-Params.beta);

C=zeros(Params.N,2);%Consumption policy
Knext=zeros(Params.N,2);%Savings policy
Params.opts=optimset();%we need this for fminbnd

for iter=1:Params.maxIter
    ValueFunction_=ValueFunction;%our "guess" for the value function (we put it on the right hand side of the Bellman eq.)
    for ind_k=1:Params.N%go through all grid points for capital
        for ind_s=1:Params.n_shocks %... and all shocks
            state=[Params.k_grid(ind_k) ind_s];
            [c,knext,val]=Optimise(state,ValueFunction_,Params);
            C(ind_k,ind_s)=c;%update the policy functions
            Knext(ind_k,ind_s)=knext;
            ValueFunction(ind_k,ind_s)=val;%...and the value function
        end
    end
    disp(norm(ValueFunction-ValueFunction_))
    if norm(ValueFunction-ValueFunction_)<0.01%convergence criteria
        break %stop value function iteration
    end
end

%collect the results to the output structure
out.Bellman.ValueFunction=ValueFunction;
out.Bellman.C=C;
out.Bellman.Knext=Knext;
out.Params=Params;
out.Simulated=Simulate(Knext,Params,(Params.k_min+Params.k_max)*0.5);%simulate and add results to the output structure

end

function S=Simulate(PolicyKnext,Params,k_init)
k=k_init; %first period capital
i_s=1; %first period shock
S.k=zeros(Params.nSimulatedPeriods,1);S.i_s=S.k;S.y=S.k;% preallocating vectors
for t=1:Params.nSimulatedPeriods
    S.k(t)=k;
    S.i_s(t)=i_s;
    S.y(t)=output(k,i_s,Params);
    knext=interp1(Params.k_grid,PolicyKnext(:,i_s),k);%next period capital interpolated from the policy function
    i_s=randsample([1;2],1,true,Params.P(i_s,:)); %one way of drawing a random shock based on Params.P, see MATLAB documentation for randsample
    k=knext;
end
end

function[c,knext,val]=Optimise(state,ValueFunction,Params)
k=state(1,1);
s=state(1,2);
y=output(k,s,Params);
knext_min=Params.k_min;
knext_max=min(Params.k_max,y+(1-Params.delta)*k-1e-5);%max knext, must afford positive consumption
[knext,val]=fminbnd(@value,knext_min,knext_max,Params.opts,state,ValueFunction,Params);%fminbnd minimizes the output argument of "value" by choosing k subject to k_min<k<k_max
c=y+(1-Params.delta)*k-knext;
val=-val;
end

function val=value(knext,state,ValueFunction,Params)
k=state(1,1);s=state(1,2);
c=output(k,s,Params)+(1-Params.delta)*k-knext;%determine consumption from the aggregate resource constraint
val=util(c,Params);%current utility
for j=1:Params.n_shocks %get next period value given next period shock and sum up to expected value
    val=val+Params.beta*Params.P(s,j)*interp1(Params.k_grid,ValueFunction(:,j),knext,Params.interptype);
end
val=-val; %need to multiply by minus one because fminbdn minimizes.
end

function y=output(k,s,Params)
y=exp(Params.z(s))*k.^Params.alpha;
end

function z=util(c,Params)
if Params.sigma==1
    z=log(c);
else
    z=c.^(1-Params.sigma)/(1-Params.sigma);
end
end
