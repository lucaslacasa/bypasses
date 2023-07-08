function [X P Q S]=communicability_shortest_path(A,beta,s,t)

%Communicability_shortest_path       
%	     Generates the communicability shortest path between two verties in a graph
%        and visualises it by comparing with the shortest (topological) one, in case 
%        they differ. 
%          
% 
%        The entry (i,j) of the matrix X corresponds to the 
%        communicability distance between the nodes i and j in the graph. 
%
%        Then, it obtained the communicability shortest path in a (directed)
%        network between the nodes s and t.
% 
%   Input         A: adjacency matrix 
%                 s: source node
%                 t: end node
% 
%   Output        X: n by n symmetric hollow matrix of communicability 
%                 distances. 
%                 P: nodes in the shortest communicability path starting at
%                 s and ending at t. 
%                 Q: nodes in the shortest path starting at s and ending at t.
%                 S: length of the shortest communicability path.  
%    
%   Reference:   Estrada, Ernesto, The communicability distance in graphs. 
%                Linear Algebra and its Applications, 436 2012, 4317-4328 
% 
%   Example: [X P Q S] = communicability_shortest_path(A,1,13,15);  



% Precalculations              
    A = max(A,A')-diag(diag(A));   
    
    
    n = length(A); 
    u = ones(n,1);
    
          
% Communicability distance matrix
    
[X An R]=communicability_geom(A, beta);

% Communicability shortest path

    B = X.*A;	% Weighted adjacency matrix based on communicability distance matrix 	
    B=max(B,B');
    H = graph(B);                     
    P = shortestpath(H,s,t);            % Shortest communicability path between s and t
    G=graph(A);
    Q = shortestpath(G,s,t);
% Length of the shortest communicability path

    L = length(P)-1; 
    S =0;                  
    
   for i=1:L
      S=S+X(P(i),P(i+1));
   end;
   
   % Visualization of the SCP and the SP between the origin-destination
   % selected
   
   p = plot(G,'Layout','force','Iterations',500,'MarkerSize',2,'LineWidth',2);
   p.EdgeColor='c';

     v=diag(expm(A));
     v_max=max(v);

    p.MarkerSize=10*v/v_max,
    p.NodeCData = v;
    p.NodeLabel = [];


highlight(p,P,'EdgeColor','g','LineWidth',8)
highlight(p,Q,'EdgeColor','r','LineWidth',8)

colormap jet;
colorbar 
set(gca,'visible','off');
set(gca,'LooseInset',get(gca,'TightInset'));


axis equal

