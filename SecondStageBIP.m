% Function that solves the second stage that is the BIP of the two steps
% procedure.
%
%
%
% function [u]=Second_Stage_BIP(Ui,Delta_Ci,L,x0,A,B,Bw,W,xmax,xmin)
%
% Inputs:
%   Ui:         Optimal Energy to be deployed in the BIP. Matrix Mx1 with
%               at least one element greater than zero
%   Delta_Ci:   Time interval for the BIP. Integer
%   L:          Indirect samping time. Integer >=1. Ts=Ui/L
%   x0:         Initial state value. Vector nx1
%   A:          Dynamic matrix. Matrix nxn
%   B:          Control Matrix. Matrix nxM
%   Bw:         Disturbance Matrix. Matrix nxD
%   W:          Integral of disturbance at fast sampling. Matrix Dx1
%   xmax:       Maximum state constraint. Real
%   xmin:       Minimum state constraint. Real
% Outputs:
%   u:          Binary control action. Matrix MxK


function [u, Umisum, L, K] = SecondStageBIP( Ui, Delta_Ci, L, x0, A, B, Wi, xmax, xmin)


% variable definition
	N  = size(Wi, 2);
	[n, M] = size(B);


	% 1.- SAMPLING TIME DEFINITION
	Ui_dummy = Ui;
	if Ui_dummy == 0 
		Umisum = [ 0; 0 ];
		u = zeros(Delta_Ci,2);
	    return 
	end

	Ui_dummy(Ui_dummy == 0) = [];
	Ui_greatest = Ui_dummy(1);
	for k = 2:length(Ui_dummy)
	    Ui_greatest  = gcd(Ui_greatest,Ui_dummy(k)); % The greatest common divisor. Required for equating rational values.
	    %Ui_greatest = lcm(Ui_greatest,Ui_dummy(k)); % The greatest common divisor. Required for equating rational values.
	end
	%L  = ceil(Ui_greatest/2);
	Ts = Ui_greatest/L;
	 if Ts < 1
	 	L  = 1;
	 	Ts = Ui_greatest;
	 end
	 %L=1;
	 
	 K  = Delta_Ci*L/Ui_greatest;       % Number of sampling intervals.
	 K  = floor(K);                     % It must be an integer. 
	 %u  = K;

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
        aux=[aux; Ad^(K-j) * Bd];
    end
    Bbarrad=[Bbarrad aux];
    aux = [zeros( n*(i), M)];
end

Iden = ones(K,1);

Ikk=[];
for i=1:K
    Ikk=[Ikk; Ad^(i-1) * x0];
end

wk =  repmat( Wi/K, 1, K);

Wbar = [];
for i = 1:K
    Wbar = [ Wbar; sum( wk( :, 1:i ), 2 ) ];             
end


Alp = [ -Bbarrad; Bbarrad ];
blp = [-kron(xmin, Iden) + Ikk + Wbar; kron(xmax,Iden) - Ikk - Wbar ];

% solver call

SumMatrix = [ repmat([1, 0], 1, K);repmat([0, 1], 1, K)];

u = binvar( M*K, 1);
Constraints = [  SumMatrix * u == (Ui * K/Delta_Ci), Alp * u <= blp  ];

solvesdp(Constraints);
 
uu = double(u);

u=[];
for k = 1:2:2*K
    u = [u; uu(k:k+1)'];       % Extraer las acciones de control U a matriz MxN
end

%Discretización al minuto del estado de la bomba
u = kron( u, ones(Delta_Ci/K,1));

Umisum = SumMatrix* uu * Delta_Ci/K;

