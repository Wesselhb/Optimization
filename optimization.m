x0 = 120*ones(1, 120);

fmincon(@cost, x0, [], [])