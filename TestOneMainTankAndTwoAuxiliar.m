addpath(genpath('~/Downloads/solvers'))

%%NOTE: YALMIP AND SOLVERS PATHS MUST BE INSTALLED AND LOADED

% Test of function First_Stage_LP
close all
clear all
% Data Input

% PROBLEM PARAMETERS

C  = [11.87 14.11 82.05 14.11 82.05 11.87];  % Electricity rates 
%C=[11.87 14.11 15.05 14.11 15.05 11.87]; 
N  = length(C);
Delta_C = 60*[6 1 3 8 4 2];                  % length in minutes of the constant price intervals


% System Dynamics
A  = zeros(3,3);                                                    % La dinamica es la de un integrador               
B1 = -0.5;                                                          % Flow of Pump 1 
B2 = -0.6;                                                          % Flow of Pump 2 
B  = [B1 B2; -B1 0; 0 -B2]*(1/60);                                  % Control Actions Matrix
Bw = 2*[ 1/6; -(1/12)*1; -(1/12)*1]*(1/60);                         % Disturbance matrix [Qe Qs1 Qs2]
Disturbance = [Bw(1)*ones(1,N); Bw(2)*ones(1,N); Bw(3)*ones(1,N)];
W = Disturbance.*[Delta_C; Delta_C; Delta_C];                       % Integral
%Bw = eye(3);


% Initial Conditions and limits

x0=[2; 1; 1];                 	      % Initial State. V1 V2 V3
xmax = [4; 2; 2];                     % Volumen Maximo
xmin = [0.2; 0.2; 0.2];               % Volumen Minimo
P1   = 5;                             % Power of Pump 1 
P2   = 6;                             % Power of Pump 2 
Pmss =[P1 P2];

% ----Solver Options
ops  = sdpsettings('solver','glpk');
%ops  = sdpsettings('solver','lpsolve');
%ops = sdpsettings('solver','lpsolve','cachesolvers',1);

%First Stage Call (Linear Programming problem)
[U,Energy,Cost]=FirstStageLP(C,x0,A,B,W,Pmss,xmax,xmin,ops);



% Changes in initial state
xmax = [4; 2.5; 2.5]*10;                     % Volumen Maximo
xmin = [0.2; 0.2; 0.2];               % Volumen Minimo


Um = ceil(U);
L = 2;
ops = sdpsettings('solver','lpsolve','cachesolvers',1);
signals = [];
for i = 1:6
   	[u, Umisum ] = SecondStageBIP( Um(:, i) , Delta_C(i), L , x0, A, B , W(:,i), xmax, xmin);
   	disp('----')
  	signals = [ signals; u ];
end

%Calculate Volumes for each timestep

Volumes = [];
for i = 1:size(signals,1)
	Volumes = [ Volumes; x0' + sum(signals(1:i,:) * -B') - i*Bw'];
end

%Data required for plotting. Should be somewhere else

perMinuteRate = [11.87 * ones( 60*6 , 1); 14.11 * ones( 60*1, 1); 82.05 * ones( 60*3, 1); 14.11 * ones( 60*8, 1); 82.05 * ones( 60*4, 1 ); 11.87 * ones( 60* 2, 1) ];
nAuxiliarTanks = 2;
nMainTanks = 1;

%%---------------Plot results------------------------

PlotVolumeEvolutionAndPumpsSignals( nMainTanks, nAuxiliarTanks, Volumes, signals, perMinuteRate)