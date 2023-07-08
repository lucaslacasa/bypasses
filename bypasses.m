function [Energy Entropy S_max S_rel]=bypasses(A)

%         Bypasses calculates the global energy saving of using bypasses in      
%	  navigating a network/graph instead of using the shortest paths. 
%         That is, it first calculates the communicability distance between every
%         pair of vertices in the graph. Then, via geometrization os the graph
%         it obtain the shortest communicability path (SCP) between every pair of vertices.
%         At the same time the program calculates the shortest (topological) path (SP)
%         between every pair of vertices. Then, by using the ratio of the length of the
%         SCP to SP it obtaine the saving of using the SCP instead of the SP.
%
%         The function also calculates the walk entropy (as well as the maximum entropy and the relative one).
%
%
%         The program uses the Matlab function "communicability_geometry.m" which is provided
%         together with the function "bupasses.m".
%
%
%   Input      A: adjacency matrix
%            
%           
%
%   Output  Energy: the value of the total energy saved by using the SCP instead of SP in the network.
%            
%   
%   Reference:   Estrada, Ernesto, Gomez-Gardeñes, J, Lacasa, L. 
%                "Network bypasses sustain complexity"
%                 arXiv preprint arXiv:2207.06813.
%
%
%   Example: [E, S] = bypasses(A);


%Precalculations

A=max(A,A');
n=length(A);
beta=1;
[X, An, R]=communicability_geom(A, beta);
X=max(X,X');

%Generation of the weighted adjacency matrix of the network

B=X.*A;
B=real(B);
B = max(B,B');

% Creation of the graphs corresponding to the original network G, and that of the communicability-distance weighted 
% version of it

G=graph(A);
G1=graph(B);

%Calculation of the matrix of SP (D) and the matrix of SCP (D1) 

D=distances(G,'Method','unweighted');
D1=zeros(n,n);

for s=1:n
for t=1:n

D1(s,t)=length(shortestpath(G1,s,t))-1;

end;end;

%Calculation of the energy saving for the use of SCP instead of SP (Energy)

O=eye(n,n);
T=100*((D+O)./(D1+O)-ones(n,n));

T(isinf(T)) = 0;
Energy=sum(sum(T))/(n*(n-1));

%Calculation of the walk entropy of the network

P=abs(expm(A));
Q=P-diag(diag(P));


ZE=0.5*sum(sum(Q));
pE=(Q./ZE)+eye(n,n);

for i=1:n
for j=1:n

if pE(i,j)==0;
pE(i,j)=1;

end;end;end

Entropy=-0.5*sum(sum(pE.*log(pE)));
S_max=log(n*(n-1)/2);

S_rel=Entropy/S_max;


