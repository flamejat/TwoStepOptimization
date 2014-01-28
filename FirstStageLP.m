function [ U, Energy, EnergyCost ] = FirstStageLP( C, x0, A, B, W, Pmss, xmax, xmin ,ops)

N = length( C );  %Intervalos del coste

% parameter declaration
U = sdpvar(2,N);
X = sdpvar(3,N);

V = [];
for i = 1:N
	V = [ V; X(:,i); U(:,i) ];
end

%objective function
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
