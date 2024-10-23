
function f=costfunc2(X)

global LAMBDA T L RHOc Vf A TAU MU K ALPHA E2 E3 Cr Dr RHOm;

    % functions
  function [qi, rhonext] = nextrho(rhoi, qprev, qr, vi)
        global LAMBDA T L RHOc;
        qi = LAMBDA*rhoi*vi; %calculate traffic flow for the current section, 
        % q(i,k)=LAMDA(i) RHO(i,k) velocity(i,k)
        rhonext = min([rhoi + (T/(LAMBDA*L))*(qprev - qi + qr) RHOc]); %calculate the traffic density for the next section 
        % rho(i,k+1) = rho(i,k) + {T/LAMBDA(i)*L(i)}{q(i-1,k)-q(i,k)+qr(i,k)}
        %if rho(i,k+1)>then the critical density RHOc then rho(i,k+1)=RHOc
    end
    
    
    function vnext = nextv(vi, vprev, rhonext, rhoi, VSL)
        global ALPHA Vf A RHOc T TAU L MU K;
        Vi = min([(1+ALPHA)*VSL; Vf*exp((-1/A)*(rhoi/RHOc)^A)]); %calculate the velocity for the section
        %V(i,k) = min([(1+ALPHA)*VSL(i,k); Vf*exp((-1/A)*(rho(i,k)/RHOc)^A)])
        posterm = vi + (T/TAU)*(Vi - vi) + (T/L)*vi*(vprev - vi);
        negterm = ((MU * T * (rhonext - rhoi))/(TAU * L * (rhoi + K)));
        vnext = posterm - negterm; %calculate v(i,k+1)
    end

    
    % Simulation
    rho = zeros(5, 120); %create a zero matrix for traffic denisty(RHO) with size [i,k]
    v = zeros(5, 120); %create a zero matrix for velocity with size [i,k]
    omega = zeros(1,120);  %create a zero matrix for OMEGA with size [1,k]
    rho(:,1) = 30 * ones(5,1); %Set the traffic denisty(RHO) for the first column (time 0) equal to 30;
    v(:,1) = 80 * ones(5,1); %Set the velocity for the first column (time 0) equal to 80;
    %It w
    for k=1:120
        if k < 60 %set the right q0 for the timestep
            q0 = 0.6*(7000 + 100*E2);
        else
            q0 = 2000 + 100*E3;
        end
    
        [qi, rho(1,k+1)] = nextrho(rho(1,k), q0, 0, v(1, k)); %Calculate the traffic flow and density for the next time step for i=1
        v(1, k+1) = nextv(v(1,k), v(1,k), rho(1+1, k), rho(1, k), 120); %calculate the velocoty for the next time step for i=1
    
        for i=2:3
            [qi, rho(i,k+1)] = nextrho(rho(i,k), qi, 0, v(i, k)); %Calculate the traffic flow and density for the next time step for i=2, 1=3
            v(i, k+1) = nextv(v(i, k), v(i-1,k), rho(i+1, k), rho(i, k), X(1,k)); %calculate the velocoty for the next time step for i=2, 1=3
        end
        
        qr = min([X(2,k)*Cr; Dr+omega(k)/T; Cr*(RHOm - rho(4,k)/(RHOm-RHOc))]);
        omega(k+1) = omega(k) + T*(Dr - qr);
        [qi, rho(4, k+1)] = nextrho(rho(4,k), qi, qr, v(4, k));
        v(4, k+1) = nextv(v(4, k), v(3,k), rho(5,k), rho(4, k), 120);
        [~, rho(5,k+1)] = nextrho(rho(5,k), qi, 0, v(5, k));
        v(5, k+1) = nextv(v(5, k), v(5-1,k), rho(5, k), rho(5,k), 120);
    end
     f = T*sum(omega) + T*L*LAMBDA*(sum(sum(rho))); %TTS = y=T omega(r,k) + T*sum(L(i)lamda(i)rho(i))
end
