A = [0 1 1 0 0 0 0;
     1 0 1 0 0 0 0;
     1 1 0 0 0 0 0;
     0 0 0 0 0 0 1;
     0 0 0 0 0 1 1;
     0 0 0 0 1 0 1;
     0 0 0 1 1 1 0];
 
G = graph(A);
plot(G)

% Número de caminhos de tamanho d entre i e j
N = A^2;

% Busca em Largura (BFS)
BFS = bfsearch(G,1);

% Componentes Conexas
conncomp(G)

t = bfsearch(G,1, 'allevents', 'Restart', true);
%visualize_search(G,t)

%%
B = [0 0 1 0 0 0 0;
     0 0 1 0 0 0 0;
     1 1 0 1 1 0 0;
     0 0 1 0 1 1 1;
     0 0 1 1 0 0 1;
     0 0 0 1 0 0 0;
     0 0 0 1 1 0 0];

 G2 = graph(B);
 plot(G2)
 sum(sum(B^2 - diag(diag(B^2))))/2
 