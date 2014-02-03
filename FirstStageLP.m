% Function that solves the first stage that is the LP of the two steps
% procedure.
%
%
%
% function [U,Energy,Cost]=First_Stage_LP(C, x0, A, B, W, Pmss, xmax, xmin)
%
% Inputs:
%   C:          Electricity price per period. Vector 1xN
%   x0:         Initial state value. Vector nx1
%   A:          Dynamic matrix. Matrix nxn
%   B:          Control Matrix. Matrix nxM
%   W:          Integral of disturbance per time interval. Matrix DxN
%   Pmss:       Steady state power consmption of each actuator. Vector 1xM
%   xmax:       Maximum state constraint. Matrix nx1
%   xmin:       Minimum state constraint. Matrix nx1
% Outputs:
%   U:          Integral control action. Matrix MxN
%   Energy:     Amount of energy consumed. Real
%   Cost:       Energy cost. Real


function [ U, Energy, EnergyCost ] = FirstStageLP( C, x0, A, B, W, Pmss, xmax, xmin)

N = length( C );  

% parameter declaration
U = sdpvar(2,N);
X = sdpvar(3,N);

V = [];
for i = 1:N
	V = [ V; X(:,i); U(:,i) ];
end

% Objective function
Objective = C * ( Pmss * U )';

% Alp matrix in Alp*x < blp
Abarra = [ A B ];
At     = kron( tril( ones( N, N), 0 ), Abarra );    % lower diagonal block matrix with A matrix in each one of its entries
Alp    = [ -At; At ];       


Iden = ones( N, 1 );

% W matrix 
Wbar = [];
for i = 1:N
    Wbar = [ Wbar; sum( W( :, 1:i ), 2 ) ];             
end

% blp matrix in Alp * x < blp
blp = [ kron( -Iden, (xmin-x0) ); kron( Iden,(xmax-x0)) ] + [ Wbar; -Wbar ];     % Matriz blp, Alp x <= blp

%Optimization call
Constraints = [Alp * V <= blp, U >= 0 ];
solvesdp( Constraints, Objective, ops);

%Output management

U = double(U);
Energy       = Pmss * ( sum( U' ) )';
EnergyCost   = C * ( Pmss * U )';
