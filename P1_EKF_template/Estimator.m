function [posEst,linVelEst,oriEst,driftEst,...
          posVar,linVelVar,oriVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,driftEst,...
%    posVar,linVelVar,oriVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2019
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch
%

const = EstimatorConst();

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % initial state mean
    posEst = [0 0]; % 1x2 matrix
    linVelEst = [0 0]; % 1x2 matrix
    oriEst = 0; % 1x1 matrix
    driftEst = const.GyroDriftStartBound; % 1x1 matrix
    
    % initial state variance
    posVar = 1/4*const.StartRadiusBound^2*ones(1,2); % 1x2 matrix
    linVelVar = [0 0]; % 1x2 matrix
    oriVar = 1/12*(const.RotationStartBound*2)^2; % 1x1 matrix
    driftVar = 0; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar'; linVelVar'; oriVar; driftVar]); 
    % estimator state
    estState.xm = [posEst'; linVelEst'; oriEst; driftEst];
    % time of last update
    estState.tm = tm;
    
    return;
end

%% Estimator iteration.

% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

% --- FUNCTION DEFINITIONS --- %

% State dynamics function
    function dxp_dt = State_update(x, u)
        dxp_dt = [x(3); 
                  x(4);
                  cos(x(5))*(tanh(u(1))-const.dragCoefficient*(x(3)^2+x(4)^2));
                  sin(x(5))*(tanh(u(1))-const.dragCoefficient*(x(3)^2+x(4)^2));
                  const.rudderCoefficient*u(2);
                  0];                 
    end

% Riccati Function
    function dPp_dt = Riccati(P, Q, x, u)
        P = reshape(P, size(estState.Pm));
        Ak = [0 0 1 0 0 0;...
              0 0 0 1 0 0;...
              0 0 -2*cos(x(5))*const.dragCoefficient*x(3) -2*cos(x(5))*const.dragCoefficient*x(4) -sin(x(5))*(tanh(u(1))-const.dragCoefficient*(x(3)^2+x(4)^2)) 0;
              0 0 -2*sin(x(5))*const.dragCoefficient*x(3) -2*sin(x(5))*const.dragCoefficient*x(4) cos(x(5))*(tanh(u(1))-const.dragCoefficient*(x(3)^2+x(4)^2)) 0;
              0 0 0 0 0 0;
              0 0 0 0 0 0];
        Lk = [0 0 0;
              0 0 0;
              -cos(x(5))*const.dragCoefficient*(x(3)^2+x(4)^2) 0 0;
              -sin(x(5))*const.dragCoefficient*(x(3)^2+x(4)^2) 0 0;
              0 const.rudderCoefficient*u(2) 0;
              0 0 1];
          
        dPp_dt = reshape(Ak*P + P*Ak' + Lk*Q*Lk', [numel(P) 1]);
    end
        
% Measurement function
    function zk = measurement(x)
        zk = [sqrt( (x(1)-const.pos_radioA(1))^2 + (x(2)-const.pos_radioA(2))^2 );
              sqrt( (x(1)-const.pos_radioB(1))^2 + (x(2)-const.pos_radioB(2))^2 );
              sqrt( (x(1)-const.pos_radioC(1))^2 + (x(2)-const.pos_radioC(2))^2 );
              x(5)+x(6);
              x(5)];
    end

% --- END OF FUNCTION DEFINITIONS --- %


% Covariance of process noise
Qk = diag([const.DragNoise, const.RudderNoise, const.GyroDriftNoise]);

% Covariance of observation noise
Rk = diag([const.DistNoiseA; const.DistNoiseB; const.DistNoiseC; const.GyroNoise; const.CompassNoise]);

% prior update
Pp = reshape(estState.Pm, [numel(estState.Pm) 1]);
xp = estState.xm;

interval = [estState.tm - dt, estState.tm];

[~,x] = ode23(@(t, x) State_update(xp, actuate), interval, xp);
[~,P] = ode23(@(t, P) Riccati(Pp, Qk, xp, actuate), interval, Pp);

xk = x(end,:)';
Pk = reshape(P(end,:), size(estState.Pm));

% measurement update
normA = norm(xk(1:2)-const.pos_radioA);
normB = norm(xk(1:2)-const.pos_radioB);
normC = norm(xk(1:2)-const.pos_radioC);

Hk = [(xk(1)-const.pos_radioA(1))/normA (xk(2)-const.pos_radioA(2))/normA 0 0 0 0; 
      (xk(1)-const.pos_radioB(1))/normB (xk(2)-const.pos_radioB(2))/normB 0 0 0 0;
      (xk(1)-const.pos_radioC(1))/normC (xk(2)-const.pos_radioC(2))/normC 0 0 0 0;
      0 0 0 0 1 1;
      0 0 0 0 1 0];

Mk = eye(length(sense));

yk = sense' - measurement(xk);

% Account for missing measurments of C
if sum(isinf(sense))
    Hk(3,:) = [];
    Mk = eye(length(Mk)-1);
    Rk(:,3) = [];
    Rk(3,:) = [];
    yk(3) = [];
end

Kk = Pk*Hk'*inv(Hk*Pk*Hk' + Mk*Rk*Mk');

xm = xk + Kk*yk;
Pm = (eye(length(Kk*Hk))-Kk*Hk)*Pk;

% Get resulting estimates and variances
estState.xm = xm;
estState.Pm = Pm;

% Output quantities
posEst = xm(1:2);
linVelEst = xm(3:4);
oriEst = xm(5);
driftEst = xm(6);

Var = diag(Pm);

posVar = Var(1:2);
linVelVar = Var(3:4);
oriVar = Var(5);
driftVar = Var(6);

end