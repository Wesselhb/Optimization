% Student dependent variables
Da = [4 6 5];
Db = [1 8 8];
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

% Simulation
rho = zeros(5,  120);
v = zeros(5, 120);
omega = zeros(1,120);
rho(:,1) = 30 * ones(5,1);
v(:,1) = 80 * ones(5,1);

for k=1:120
    if k < 60
        q0 = (7000 + 100*E2)*0.6;
    else
        q0 = 2000 + 100*E3;
    end

    [qi, rho(1,k+1)] = nextrho(rho(1,k), q0, 0, v(1, k));
    v(1, k+1) = nextv(v(1, k), v(1,k), rho(1+1, k), rho(1, k), 120);

    for i=2:3
        [qi, rho(i,k+1)] = nextrho(rho(i,k), qi, 0, v(i, k));
        v(i, k+1) = nextv(v(i, k), v(i-1,k), rho(i+1, k), rho(i, k), Xbar(1,k));
    end
    
    qr = min([Xbar(2,k)*Cr; Dr+omega(k)/T; Cr*(RHOm - rho(4,k)/(RHOm-RHOc))]);
    omega(k+1) = omega(k) + T*(Dr - qr);
    [qi, rho(4, k+1)] = nextrho(rho(4,k), qi, qr, v(4, k));
    v(4, k+1) = nextv(v(4, k), v(3,k), rho(5,k), rho(4, k), 120);
    [~, rho(5,k+1)] = nextrho(rho(5,k), qi, 0, v(5, k));
    v(5, k+1) = nextv(v(5, k), v(5-1,k), rho(5, k), rho(5,k), 120);
    TTS(1,k) = T*omega(k) + T*L*LAMBDA*(sum(rho(:,k)));
end

 

%plotting
figure()
title('Velocity')
xlabel('Time')
ylabel('Velocity')
hold on
for i=1:5
    plot(v(i,:))
end
legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5')


figure()
title('density')
xlabel('Time')
ylabel('density')
hold on
for i=1:5
    plot(rho(i,:))
end
legend('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5')
%% Plot Q4
figure()
title('TTS')
xlabel('Time')
ylabel('Total time spent')
hold on
plot(TTS)

%% Plot Q5
figure()
title('VSL')
xlabel('Time')
ylabel('Speed Limit')
hold on
plot(Xbar(1,:))

figure()
title('ramp metering rate applied at time step k')
xlabel('Time')
ylabel('Ramp metering rate')
hold on
plot(Xbar(2,:))

