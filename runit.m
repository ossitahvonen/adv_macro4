%basicgrowthmodel(1.5)
R=realbusinesscycles()
%plot(R.Simulated.k)
%plot(R.Simulated.i_s)
%plot(R.Simulated.y)

res = ["Variable", "Mean", "SD";
    "Capital",mean(R.Simulated.k),std(R.Simulated.k);
    "Investment",mean(R.Simulated.inv),std(R.Simulated.inv);
    "Consumption",mean(R.Simulated.c),std(R.Simulated.c)]