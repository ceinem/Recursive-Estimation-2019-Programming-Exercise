function trackErrorNorm=run(simConst,estConst,doplot,seed)
%
%
% Main function for the Particle Filter programming exercise.
% It is used to simulate the true model, call the estimator and show the 
% results.
%
%
% Class:
% Recursive Estimation
% Spring 2019
% Programming Exercise 2
%
% ----------------------------------------------------------------------- %
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch
% ----------------------------------------------------------------------- %
% clear command window, close figures
clc;
close all;

if nargin==0
   % Define the simulation constants that are used in simulation, but not 
   % accessible to the estimator.  
   simConst = SimulationConst();
   
   % Define the constants accessible to the estimator.
   estConst = EstimatorConst();
   
   % Generate plots by default.
   doplot = true;
   
   % use random seed
   seed = 0;
end

%% Setup
% Set the random number generator state.
% Uncomment to make results reproducable.
if seed ~= 0
    rand('seed',seed);
    randn('seed',seed);
end

%% Simulation
% The function 'Simulator' simulates the robot dynamics and generates
% measurements.
[km, state, input, sense] = Simulator( simConst );

% state contains rolling out of [x_r, y_r, phi] for k={1,N}
% input contains rolling out of [u_f, u_phi] for k={1,N-1}
% sense contains rolling out of [z] for k={1,N}, where first measurement is
% infinity (to indicate no measurement at initial step)

%% Run the Estimator
N = simConst.N; 

% Initialize the estimator.
[posterior] = Estimator([],[],[],estConst,km(1));

% Get number of particles:
N_particles = size(posterior.x_r,2);

% initialize storage arrays for x, y and phi:
x_r = zeros(N, N_particles);
x_r(1,:) = posterior.x_r;
y_r = zeros(N, N_particles);
y_r(1,:) = posterior.y_r;
phi = zeros(N, N_particles);
phi(1,:) = posterior.phi;

% Call the estimator for the remaining time steps and measure the total 
% execution time (including the updates of the storage arrays):
if doplot
    disp('Generating data...')
end
nextDispFracComplete = 0.1;
tstart = tic;
for k = 2:N
    if doplot
        if (k/N ) > nextDispFracComplete
            nextDispFracComplete = nextDispFracComplete + 0.1;
            disp([num2str(round(k/N*100)),'% complete'])
        end
    end
    [posterior] = Estimator(posterior,sense(k,:)',input(k-1,:)',estConst,km(k));
    % update storage arrays:
    x_r(k,:) = posterior.x_r;
    y_r(k,:) = posterior.y_r;
    phi(k,:) = posterior.phi;
    
end
tEstAvg = toc(tstart)/(N-1)
if doplot
    disp('Done. Making plots.')
end

%% Plotting
if doplot
    figure(1)
    % Plot countour
    spacing = 0.6;
    contour = simConst.contour;
    axis([min(contour(:,1))-spacing,max(contour(:,1))+spacing,...
          min(contour(:,2))-spacing,max(contour(:,2))+spacing])
    hold all
    for i = 1:size(contour,1)-1
        plot([contour(i,1),contour(i+1,1)],[contour(i,2),contour(i+1,2)],'r--')
    end
    plot([contour(end,1),contour(1,1)],[contour(end,2),contour(1,2)],'r--')

    % Plot particles
    robot_points = 0.05*[-1,-1;-1,1;0.5,1;1,0;0.5,-1];
    robot_area = fill(robot_points(:,1)',robot_points(:,2)','r');
    distance_line = fill([0,0],[0,0],'r');
    hold all
    particle_circles = plot(NaN,NaN,'ro');
    pbaspect([1 1 1])
    hold off
    for k = 1:size(state,1)
        robot_area.XData = cos(state(k,3))*robot_points(:,1)...
            -sin(state(k,3))*robot_points(:,2)+state(k,1);
        robot_area.YData = sin(state(k,3))*robot_points(:,1)...
            +cos(state(k,3))*robot_points(:,2)+state(k,2);
        particle_circles.XData = x_r(k,:);
        particle_circles.YData = y_r(k,:);
        if k > 1
            distance_line.XData = cos(state(k,3))*[0;sense(k,:)]+state(k,1);
            distance_line.YData = sin(state(k,3))*[0;sense(k,:)]+state(k,2);
        end
        drawnow
    end 
end

%% The results
% initialize tracking error
trackError = zeros(N,1);

% Calculate the total tracking error as an indicator of your estimator's performance.
for k = 1:N
    trackError(k) = sqrt(mean((x_r(k,:) - state(k,1)).^2 + (y_r(k,:) - state(k,2)).^2));
end
trackErrorNorm = sqrt(trackError'*trackError/N);

return;