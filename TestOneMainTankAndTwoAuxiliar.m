%% Two step optimization algorithm
%  A case study that shows the applicability of this new approach to the solution of optimization problems with discrete control variables modeled by binary integer programs (BIP).
%  The solution of BIP is computationally demanding when the number of the BIP variables increase. Two instances that increase the number BIP variables in practical applications
%  are the reduction of the discretization sampling time and the increase of the optimization time period. The proposed approach transforms a single BIP optimization into a linear
%  program (LP) and $N$ feasibility BIP's, with less number of variables. The reduction of the number of variables increases the algorithm speed in providing a solution. The 
%  approach permits to solve optimization problems with longer time intervals and with a higher number of control variables, while being computationally tractable.
%%

%% Setting up the problem
%
% NOTE: YALMIP AND SOLVERS PATHS MUST BE INSTALLED AND LOADED
close all
clear all

% PROBLEM PARAMETERS

C  = [11.87 14.11 82.05 14.11 82.05 11.87];  % Electricity rates 
N  = length(C);
Delta_C = 60*[6 1 3 8 4 2];                  % length in minutes of the constant price intervals
perMinuteRate = [11.87 * ones( 60*6 , 1); 14.11 * ones( 60*1, 1); 82.05 * ones( 60*3, 1); 14.11 * ones( 60*8, 1); 82.05 * ones( 60*4, 1 ); 11.87 * ones( 60* 2, 1) ];
nAuxiliarTanks = 2;
nMainTanks = 1;

% System Dynamics

A  = zeros(3,3);                                                    % La dinamica es la de un integrador               
B1 = -0.5;                                                          % Flow of Pump 1 
B2 = -0.6;                                                          % Flow of Pump 2 
B  = [B1 B2; -B1 0; 0 -B2]*(1/60);                                  % Control Actions Matrix
Bw = 2*[ 1/6; -(1/12)*1; -(1/12)*1]*(1/60);                         % Disturbance matrix [Qe Qs1 Qs2]
Disturbance = [Bw(1)*ones(1,N); Bw(2)*ones(1,N); Bw(3)*ones(1,N)];
W = Disturbance.*[Delta_C; Delta_C; Delta_C];                       % Integral


% Initial Conditions and limits

x0=[2; 1; 1];                 	      % Initial State. V1 V2 V3
xmax = [4; 2; 2];                     % Maximum Volume
xmin = [0.2; 0.2; 0.2];               % Minimum Volume
P1   = 5;                             % Power of Pump 1 
P2   = 6;                             % Power of Pump 2 
Pmss =[P1 P2];

% Solver Options
ops = sdpsettings('solver','glpk');
ops2  = sdpsettings('solver','lpsolve');
ops3 = sdpsettings('solver','lpsolve','cachesolvers',1);


%% Test of the First Stage (LP)
% DESCRIPTIVE TEXT

%First Stage Call (Linear Programming problem)
[U,Energy,Cost]=FirstStageLP(C,x0,A,B,W,Pmss,xmax,xmin, ops);


%% Test of second Stage (BIP)
% 

% Changes in initial state
xmax = [4; 2.5; 2.5]*10;              % Maximum Volume
xmin = [0.2; 0.2; 0.2];               % Minimum Volume

% Rounding necessary for the BIP stage. 
Um = ceil(U);

% An L is fixed. L is the parameter that will determine the discretization timesteps in each one of the intervals.

L = 2;

% Obtain the optimal pumping in each electricity rate interval

signals = [];
for i = 1:6
   	u = SecondStageBIP( Um(:, i) , Delta_C(i), L , x0, A, B , W(:,i), xmax, xmin, ops2);
  	signals = [ signals; u ];
end

%Calculate Volumes in each timestep

Volumes = [];
for i = 1:size(signals,1)
	Volumes = [ Volumes; x0' + sum(signals(1:i,:) * -B') - i*Bw'];
end

%% Plot results
% Plot the volume evolution of each tank and the pumping.

PlotVolumeEvolutionAndPumpsSignals( nMainTanks, nAuxiliarTanks, Volumes, signals, perMinuteRate)