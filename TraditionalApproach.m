% Function that implements the traditional BIP approach used for this kind of problems
%
%
%
% function [U, Energy, Cost] = TraditionalApproach(C, x0, A, B, W, Pmss, xmax, xmin)
%
% Inputs:
%   C:          Electricity price per period. Vector 1xN.
%   x0:         Initial state value. Vector nx1.
%   A:          Dynamic matrix. Matrix nxn.
%   B:          Control Matrix. Matrix nxM.
%   W:          Integral of disturbance per time interval. Matrix DxN.
%   Pmss:       Steady state power consmption of each actuator. Vector 1xM.
%   xmax:       Maximum state constraint. Matrix nx1.
%   xmin:       Minimum state constraint. Matrix nx1.
%   Ts:			sampling time. Integer
%   ops:		Yalmip's solver options.
% Outputs:
%   U:          Integral control action. Matrix MxN.
%   Energy:     Amount of energy consumed. Real.
%   Cost:       Energy cost. Real.

function [ u, Energy, EnergyCost ] = TraditionalApproach( Cost, Delta_C, Ts , x0, A, B , W, xmax, xmin, Pmss,ops);

% 0. - VARIABLE DEFINITION
N  = size(W, 2);
[n, M] = size(B);

C=[]
for i= 1: length(Cost)
    C = [C; Cost(i) * ones(M, 1)];
end

% 1.- SAMPLING TIME DEFINITION

K  = sum(Delta_C);       % Number of sampling intervals.


% 2.- DISCRETIZATION OF SYSTEM DYNAMICS
Scd = ss( A, B, eye(n), zeros(n, M));
Sdd = c2d(Scd,Ts,'zoh');
Ad  = Sdd.a;
Bd  = Sdd.b;

% 3.- CONSTRAINTS CONSTRUCTION
% Ax<=b
% Abarrad

Bbarrad=[];
aux=[];
for i = 1:K
    for j = K:-1:i
        aux=[ aux; Ad^(K-j) * Bd ];
    end
    Bbarrad=[ Bbarrad aux ];
    aux = [ zeros( n*(i), M) ];
end

Iden = ones(K,1);

Ikk=[];
for i=1:K
    Ikk = [Ikk; Ad^(i-1) * x0];
end

wk =  repmat( W/Ts, 1, K);

Wbar = [];
for i = 1:K
    Wbar = [ Wbar; sum( wk( :, 1:i ), 2 ) ];             
end


Alp = [ -Bbarrad; Bbarrad ];
blp = [-kron(xmin, Iden) + Ikk + Wbar; kron(xmax,Iden) - Ikk - Wbar ];

% 4.- SOLVER CALL


u = binvar( M*K, 1);

Constraints = [ Alp * u <= blp  ];
Objective =  C' * ( repmat(Pmss,1, K)' .* u );
solvesdp(Constraints, Objective,ops);
 
uu = double(u);

u=[];
for k = 1:2:2*K
    u = [u; uu(k:k+1)'];       % Extracts control actions u to an MxN matrix
end


% Minute-wise pump state discretization

u = kron( u, ones(K,1));
EnergyCost = C' * (repmat(Pmss, 1, K)' .* uu )
Energy = 0;
%Energy     = Pmss * uu ;