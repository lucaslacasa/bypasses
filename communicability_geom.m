function [X An R]=communicability_geom(A, beta)

%Communicability_angle      
%	      Generates the matrices An and X of communicability angles 
%         and distances of a network.
%         The (i,j) entry of the matrix An corresponds to the                        
%         communicability angle between the position vectors of the nodes 
%         i and j in a hyperspherical embedding of the graph.
%
%         The entry (i,j) of the matrix X corresponds to the
%         communicability distance between the nodes i and j in the graph.
%
%   Input      A: adjacency matrix
%           beta: inverse temperature. Defaults to 1 
%           
%
%   Output  An: n by n symmetric hollow matrix of communicability angles.
%            X: n by n symmetric hollow matrix of communicability
%            distances.
%   
%   Reference:   Estrada, Ernesto, and Naomichi Hatano. 
%                "Communicability Angle and the Spatial Efficiency 
%                 of Networks." SIAM Rev (2016).
%
%   Example: [An, X] = communicability_angle(A,1);
 

if nargin <= 1
    beta = 1;
end;

% Precalculations
            
A=max(A,A')-diag(diag(A));  
n=length(A);
u=ones(n,1);

% Communicability

[V L]=eig(A);
G=(V*(expm(beta*L))*V');


%G=expm(beta*A);                   % Communicability matrix
sc=(diag(G));                       % Vector of self-communicabilities

% Comunicability angles matrix

An=acosd(G./((sc*u').*(u*sc')).^0.5);


% Communicability distance matrix

CD=(sc*u'+u*sc'-2*G);            %Squared Communicability distance matrix
X=CD.^0.5;                       %Communicability distance matrix
X_mean=mean(mean(X));

GG=expm(-beta*A);
s=diag(G);
a=ones(1,n)*GG*ones(n,1);
b=ones(1,n)*GG*s;
c=s'*GG*s;

R=(0.25*(c-(2-b)^2/a))^0.5;
