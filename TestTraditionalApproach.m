%% Test the traditional approach
% This script intents to show the traditional approach used to solve this kind 
% of problems and compare with our proposal.

%% Setting up the problem
%
% NOTE: YALMIP AND SOLVERS PATHS MUST BE INSTALLED AND LOADED
close all
clear all

% PROBLEM PARAMETERS

% Electicity rates parameters
Ts = 30;
Rates  = [11.87 14.11 82.05 14.11 82.05 11.87];  % Electricity rates 
Delta_C = 60/Ts*[6 1 3 8 4 2];                  % length in minutes of the constant price intervals
perMinuteRate = [11.87 * ones( 60*6 , 1); 14.11 * ones( 60*1, 1); 82.05 * ones( 60*3, 1); 14.11 * ones( 60*8, 1); 82.05 * ones( 60*4, 1 ); 11.87 * ones( 60* 2, 1) ];

% Cost of electricity in each timestep.
C=[]
for i = 1:length(Rates)
	C =[ C, Rates(i).*ones( 1, Delta_C(i))]
end

nAuxiliarTanks = 2;
nMainTanks = 1;

% System dynamics
A  = zeros(3,3);                                                    % La dinamica es la de un integrador               
B1 = -0.5;                                                          % Flow of Pump 1 
B2 = -0.6;                                                          % Flow of Pump 2 
B  = [B1 B2; -B1 0; 0 -B2]*(1/60);                                  % Control Actions Matrix
Bw = 2*[ 1/6; -(1/12)*1; -(1/12)*1]*(1/60);                         % Disturbance matrix [Qe Qs1 Qs2]
W = Bw * Ts;


% Initial Conditions and limits

x0 = [2; 1; 1];                 	      % Initial State. V1 V2 V3
xmax = [4; 2; 2];                     % Maximum Volume
xmin = [0.2; 0.2; 0.2];               % Minimum Volume
P1   = 5;                             % Power of Pump 1 
P2   = 6;                             % Power of Pump 2 
Pmss = [P1 P2];

% Solver Options
ops = sdpsettings('solver','glpk');
ops2  = sdpsettings('solver','lpsolve');
ops3 = sdpsettings('solver','lpsolve','cachesolvers',1);


[u, Energy, EnergyCost] = TraditionalApproach( C, Delta_C, Ts , x0, A, B , W, xmax, xmin, Pmss, ops);

Volumes = [];
for i = 1:size(u,1)
	Volumes = [ Volumes; x0' + sum(u(1:i,:) * -B') - i*Bw'];
end


%PlotVolumeEvolutionAndPumpsSignals( nMainTanks, nAuxiliarTanks, Volumes, u, perMinuteRate)