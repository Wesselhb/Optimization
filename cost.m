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

z=1:120;
lb = 60*ones(1,120);
ub = 120*ones(1,120);

x = fmincon(@costfunc2,z,[],[],[],[],lb,ub);



function f=costfunc2(V_SL)

global LAMBDA T L RHOc Vf A TAU MU K ALPHA E2 E3 Cr Dr RHOm;

    % functions
        function [qi, rhonext] = nextrho(rhoi, qprev, qr, vi)
            qi = LAMBDA*rhoi*vi;
            rhonext = min([rhoi + (T/(LAMBDA*L))*(qprev - qi + qr) RHOc]);
        end
    
    
        function vnext = nextv(vi, vprev, rhonext, rhoi, VSL)
            Vi = min([(1+ALPHA)*VSL; Vf*exp((-1/A)*(rhoi/RHOc)^A)]);
            posterm = vi + (T/TAU)*(Vi - vi) + (T/L)*vi*(vprev - vi);
            negterm = ((MU * T * (rhonext - rhoi))/(TAU * L * (rhoi + K)));
            vnext = posterm - negterm;
        end

    
    % Simulation
    rho = zeros(5, 120);
    v = zeros(5, 120);
    omega = zeros(1,120);
    rho(:,1) = 30 * ones(5,1);
    v(:,1) = 80 * ones(5,1);
    
    for k=1:120
        if k < 60
            q0 = 7000 + 100*E2;
        else
            q0 = 2000 + 100*E3;
        end
    
        [qi, rho(1,k+1)] = nextrho(rho(1,k), q0, 0, v(1, k));
        v(1, k+1) = nextv(v(1, k), v(1,k), rho(1+1, k), rho(1, k), 120);
    
        for i=2:3
            [qi, rho(i,k+1)] = nextrho(rho(i,k), qi, 0, v(i, k));
            v(i, k+1) = nextv(v(i, k), v(i-1,k), rho(i+1, k), rho(i, k), V_SL(k));
        end
        rx = 1;
        qr = min([rx*Cr; Dr+omega(k)/T; Cr*(RHOm - rho(4,k)/(RHOm-RHOc))]);
        omega(k+1) = omega(k) + T*(Dr - qr);
        [qi, rho(4, k+1)] = nextrho(rho(4,k), qi, qr, v(4, k));
        v(4, k+1) = nextv(v(4, k), v(3,k), rho(5,k), rho(4, k), 120);
        [~, rho(5,k+1)] = nextrho(rho(5,k), qi, 0, v(5, k));
        v(5, k+1) = nextv(v(5, k), v(5-1,k), rho(5, k), rho(5,k), 120);
    end
     f = T*sum(omega) + T*L*LAMBDA*(sum(sum(rho)));
end
