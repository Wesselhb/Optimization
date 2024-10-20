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
%% Q4
z=ones(2,120);
z(1,:)=80*z(1,:);
lb = 60*ones(2,120);
lb(2,:) = 1*lb(1,:);
ub = ones(2,120);
ub(1,:) = 120*ub(1,:);

VSL_min = fmincon(@costfunc2,z,[],[],[],[],lb,ub);


%% Q5
z=ones(2,120);
z(1,:)=80*z(1,:);
lb = 60*ones(2,120);
lb(2,:) = 0*lb(1,:);
ub = ones(2,120);
ub(1,:) = 120*ub(1,:);

Xbar = fmincon(@costfunc2,z,[],[],[],[],lb,ub);

%% Q8
z=ones(2,120);
z(1,:)=80*z(1,:);
lb = 60*ones(2,120);
lb(2,:) = 0*lb(1,:);
ub = ones(2,120);
ub(1,:) = 120*ub(1,:);

Hbar = ga(@costfunc2,240,[],[],[],[],lb,ub);
