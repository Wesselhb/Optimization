% Student dependent variables
Da = [4 6 5];
Db = [1 8 8];
global E2 E3
E1 = Da(1) + Db(1);
E2 = Da(2) + Db(2);
E3 = Da(3) + Db(3);

% Constants
global TAU MU Cr RHOm ALPHA K A Vf RHOc T L LAMBDA Dr;
TAU = 10/3600; % [hr]
MU = 80; % [km^2/hr]
Cr = 2000; % [veh/hr]
RHOm = 120; % [veh/hr/lane]
ALPHA = 0.1; % []
K = 10; % [veh/lane/km]
A = 2; % []
Vf = 110; % [km/hr]
RHOc = 33.5 + E1/3; % [veh/km/lane]
T = 10/3600; % [hr]
L = 1; % [km]
LAMBDA = 3; % [lanes]
Dr = 1500;

options = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 3e4);
%% Q4
z=ones(2,120);
z(1,:)=60*z(1,:);
lb = ones(2,120);
lb(1,:) = 60*lb(1,:);
ub = ones(2,120);
ub(1,:) = 120*ub(1,:);

Xbar = fmincon(@costfunc2,z,[],[],[],[],lb,ub,[],options);


%% Q5
z=ones(2,120);
z(1,:)=120*z(1,:);
lb = 60*ones(2,120);
lb(2,:) = 0*lb(1,:);
ub = ones(2,120);
ub(1,:) = 120*ub(1,:);
Xbar = fmincon(@costfunc2,z,[],[],[],[],lb,ub,[], options);

%% Q8

lb = ones(1,240);
ub = 4*ones(1,240);
Hbar = ga(@costfunc3,240,[],[],[],[],lb,ub,[],1:240);

for i=1:120
    Xbar(1,i) = 40+20*Hbar(i);
    Xbar(2, i) = Hbar(2*i)/5;
end
    
