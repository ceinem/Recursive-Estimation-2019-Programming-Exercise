function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==1, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index, scalar
%                       corresponds to continous time t = k*Ts
%                       If tm==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%
%
% Class:
% Recursive Estimation
% Spring 2019
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
%

% --- FUNCTIONS ---

    function [x,y] = init_pos_pdf(N)
        % First Circle A
        thetaA = rand(1,floor(N/2))*2*pi;
        rA = estConst.d*sqrt(rand(1,floor(N/2)));
        xA = rA.*cos(thetaA) + estConst.pA(1);
        yA = rA.*sin(thetaA) + estConst.pA(2);
        % Second Circle B
        thetaB = rand(1, ceil(N/2))*2*pi;
        rB = estConst.d*sqrt(rand(1, ceil(N/2)));
        xB = rB.*cos(thetaB) + estConst.pB(1);
        yB = rB.*sin(thetaB) + estConst.pB(2);
        % Put together to get x,y
        x = [xA xB];
        y = [yA yB];
    end

    function xk = state_update(x, u, v)
        xk = [x(1,:) + (u(1) + v(1,:)).*cos(x(3,:));
              x(2,:) + (u(1) + v(1,:)).*sin(x(3,:));
              x(3,:) + u(2) + v(2,:)];
    end


%     function d = distanceToWall(x, room)
%         far = 10;
%         
%         line = [x(1), x(2); x(1)+far*cos(x(3)), x(2)+far*sin(x(3))];
%         %tic
%         if isinterior(room, x(1), x(2))
%             tic
%             [~,C] = intersect(room, line);
%             toc
%             tic
%             d = norm(x(1:2)-C(1,:)');
%             toc
%         else
%             d = far;
%         end
%         %toc
%     end

    function C_intersect = intersect_points(P)
        far = 20;
        C_intersect = zeros(2,N_particles)*far;
        interior = isinterior(room, P(1,:), P(2,:));
        for j = 1:N_particles
            if interior(j)
                line = [P(1,j), P(2,j); P(1,j)+far*cos(P(3,j)), P(2,j)+far*sin(P(3,j))];
                [~,C] = intersect(room, line);
                C_intersect(:,j) = [C(1,1); C(1,2)];
            end
        end
    end

    function beta = get_weights(measurement, calculation)
        delta = measurement * ones(size(calculation)) - calculation;
        eps = estConst.epsilon;
        probabilities = [0, 1/(5*eps), 0, 2/(5*eps), 0, 1/(5*eps), 0];
        x = [-3*eps, -2.5*eps, -2*eps, 0 ,2*eps, 2.5*eps, 3*eps];

        pdf = interp1(x, probabilities, delta, 'linear');
        pdf(isnan(pdf)) = 0;
        
        beta = pdf;
        
        if (sum(beta) > 0)
            beta = beta/sum(beta);
        else
            warning('Particle weights were all zero')
            beta = ones(size(beta))/length(beta);
            K = K*10;
        end    
    end

% Set number of particles:
N_particles = 2750;
room = polyshape(estConst.contour);
K = 1e-2;

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    [postParticles.x_r, postParticles.y_r] = init_pos_pdf(N_particles); % [1xN_particles matrix, 1xN_particles matrix]
    postParticles.phi = rand(1,N_particles)*(2 * estConst.phi_0) - estConst.phi_0; % 1xN_particles matrix
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement your estimator here!

% Prior Update:
v = rand(2,N_particles) .* [estConst.sigma_f; estConst.sigma_phi] - [estConst.sigma_f/2; estConst.sigma_phi/2] ;

priorParticles = state_update([prevPostParticles.x_r; prevPostParticles.y_r; prevPostParticles.phi], act, v);

% Posterior Update:
C_intersect = intersect_points(priorParticles);
dist = vecnorm(C_intersect - priorParticles(1:2,:));

Beta = get_weights(sens, dist);

% Resample:
Particles = zeros(size(priorParticles));
rand_nums = rand(1,N_particles);
Beta_cum = cumsum(Beta);
for i = 1:N_particles
    Particles(:,i) = priorParticles(:,find(rand_nums(i) <= Beta_cum,1));
end

postParticles.x_r = Particles(1,:);
postParticles.y_r = Particles(2,:);
postParticles.phi = Particles(3,:);

% Roughening
E_x = max(postParticles.x_r)-min(postParticles.x_r);
E_y = max(postParticles.y_r)-min(postParticles.y_r);
E_phi = max(postParticles.phi)-min(postParticles.phi);
d = 3;
 
sigma_x = K * E_x * N_particles^(-1/d);
sigma_y = K * E_y * N_particles^(-1/d);
sigma_phi = K * E_phi * N_particles^(-1/d);

rough_x = normrnd(0,sigma_x,1,N_particles);
rough_y = normrnd(0,sigma_y,1,N_particles);
rough_phi = normrnd(0,sigma_phi,1,N_particles);

postParticles.x_r = postParticles.x_r + rough_x;
postParticles.y_r = postParticles.y_r + rough_y;
postParticles.phi = postParticles.phi + rough_phi;

% Use only if no measurement update
% postParticles.x_r = priorParticles(1,:);
% postParticles.y_r = priorParticles(2,:);
% postParticles.phi = priorParticles(3,:);
end % end estimator