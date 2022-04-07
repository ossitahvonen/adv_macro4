%%Advanced Macroeconomics 4 2022, Niku Määttänen
%Programs to solve a basic stochastic life cycle savings model
%with a flat rate pension system. Save survivalprobs.xls in the same folder with this file.
%Write e.g. R=lifecycle() in the command window. To plot e.g. the
%average savings profile, use the command plot(mean(R.A')). To plot e.g.
%the first individual savings profíle, use the command plot(R.A(:,1)).

function out=lifecycle()
    global Params
    N=100;% number of gridpoints for the asset state
    Params.interptype='lin'; %Interpolation method (see interp1): linear is fast but requires more gridpoints
    Params.a_min=0;%%min asset holdings
    Params.a_max=20;%%max asset holdings, should be large enough so it is never close to being binding.
    Params.n_indivs=1000;%%number of simulated individuals
    Params.sigma=2.0; %%risk-aversion
    Params.beta=0.95;%%discount factor
    Params.R=1.04;%interest rate (1+r)
    Params.pension_benefit=0.3;%flat rate pension benefit. (Mean full time wage will be normalized to one.)
    Params.s=zeros(2,1);
    Params.s(1)=0.5;%the first income shock
    Params.s(2)=1.5;%the second income shock
    Params.P=[0.8 0.2; 0.2 0.8];%%the transition probability matrix for the income shock: P(i,j) is the prob that next income is j given that current shocks is i.
    Params.n_s=length(Params.s);
    
    Params.MaxAge=75;%We think of age=1 as corresponding to real age 25, so Max age corresponds to real age 99.
       for j=1:Params.MaxAge
           real_age=j+25;
           Params.w(j)=max(0,1570+47*real_age+0.80*real_age^2-0.016*real_age^3); %age-wage profile: simple regression based on data from Statistics Finland (average monthly wage, fulltime workers)
       end
    Params.RetirementAge=42;%corresponds to real age 66. No wage earnings after that.
    Params.w=Params.w/mean(Params.w(1:Params.RetirementAge));%normalize average wage for age group 25-66 to 1

    S=xlsread('survivalprobs.xls');%age-specific survival probabilities starting from real age 25
    Params.S=S(:,2);
    
    %Solve for the pension contribution rate (assuming PAYG system) 
    Params.tau=Params.pension_benefit*(Params.MaxAge-Params.RetirementAge+1)/(Params.RetirementAge-1)*(1/mean(Params.w(1:Params.RetirementAge-1)));
    
    %%Create the asset grid. We place the grid points logarithmically with logspace:
    Params.a_grid=logspace(log10(Params.a_min+1),log10(Params.a_max+1),N)-1;
        
    %Initialize value and policy functions as N*n_s*MaxAge arrays
    ValVec=zeros(N,Params.n_s,Params.MaxAge);%Value function. 
    Cons=zeros(N,Params.n_s,Params.MaxAge);%Consumption policy
    Anext=zeros(N,Params.n_s,Params.MaxAge);%Savings policy
    opts=optimset();%Need this for fminbnd

    %Value function iteration
    for age=Params.MaxAge:-1:1
        disp(age)
        for i=1:N
            for i_s=1:Params.n_s
                state=[Params.a_grid(i) i_s]; %note: "state" includes savings and income shock only, not age.
                c_min=0.001;
                c_max=Params.R*state(1,1)+Income(state,age,Params);%maximum consumption assuming no borrowing
                if age==Params.MaxAge
                    c=c_max;%In this case, always want to consume as much as possible
                    val=-util(c_max,Params);%we are going to multiply this by -1 ...
                else
                    [c,val,flag]=fminbnd(@value,c_min,c_max,opts,state,age,ValVec(:,:,age+1),Params);
                    assert(flag==1,'possible problem with fminbnd')
                end
                ValVec(i,i_s,age)=-val;%we used fminbnd to minimize -val...
                Cons(i,i_s,age)=c;
                Anext(i,i_s,age)=nextAsset(c,state,age,Params);
            end
        end
    end
    
    %Finally, simulate n_indivs life cycles 
    %Notice that we don't kill households before MaxAge.
    
    n_indivs=Params.n_indivs;
    A=zeros(Params.MaxAge,Params.n_indivs);
    C=zeros(Params.MaxAge,Params.n_indivs);
    S=zeros(Params.MaxAge,Params.n_indivs);

    disp('Simulating...')
    
    for i=1:n_indivs
        %choose initial shock randomly with equal
        %probablities for each shock
        i_s=randsample([1:1:Params.n_s]',1,true,ones(1,Params.n_s)/Params.n_s); %see MATLAB documentation for randsample
        state=[0 i_s];%initial state
        for age=1:Params.MaxAge
            c=interp1(Params.a_grid,Cons(:,state(2),age),state(1),Params.interptype,'extrap');
            A(age,i)=state(1);
            C(age,i)=c;
            S(age,i)=state(2);
            anext=nextAsset(c,state,age,Params);
            i_s_next=randsample([1:1:Params.n_s]',1,true,Params.P(i_s,:)); %next period shock 
            state=[anext i_s_next];%next state
        end
    end
    
    %collect the results into the output argument (structure)
    out=struct('A',A,'C',C,'ValVec',ValVec,'Cons',Cons,'Anext',Anext,'S',S,'Params',Params); %a concise way of defining a structure
end

function z=value(c,state,age,ValVec,Params)
    i_s=state(2);%current income shock
    anext=nextAsset(c,state,age,Params);%next asset given current state and consumption (c)
    %compute current utility plus expected discounted value given consumption
    val=util(c,Params);
    for i_s_next=1:Params.n_s %sum up next period values for different next period shocks to get expected value
    val=val+Params.S(age)*Params.beta*Params.P(i_s,i_s_next)*interp1(Params.a_grid,ValVec(:,i_s_next),anext,Params.interptype);
    end
    z=-val;%need to multiply by minus one because fminbdn minimizes.
end

function anext=nextAsset(c,state,age,Params)
    a=state(1);
    anext=a*Params.R+Income(state,age,Params)-c;
end

function z=Income(state,age,Params)
    i_s=state(2);
    if age>=Params.RetirementAge
        z=Params.pension_benefit;
    else
        z=(1-Params.tau)*Params.w(age)*Params.s(i_s);
    end
end

function z=util(c,Params)
    if(c<1e-3)
        z=-1e5;
    elseif(Params.sigma==1)
        z=log(c);
    else
        z=c^(1-Params.sigma)/(1-Params.sigma);
    end
end
