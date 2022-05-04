function[K_T]=findKT(edges,cells,MicroTermIndexes,MacroTermIndexes,q_network,p_network,p_darcy)
G = graph(edges(:,2)',edges(:,3)');
d=distances(G);
dTerm = d(MicroTermIndexes,MacroTermIndexes);
K_T = zeros(size(cells,1),length(MacroTermIndexes));
for i = 1:size(dTerm,1)
    MinDist = min(dTerm(i,:));
    index = find(dTerm(i,:)==MinDist);
    ParentNode = MacroTermIndexes(index);
    EdgeNr = find(ParentNode==edges(:,3));
    q_I = q_network(EdgeNr);
    K_T(i,index) = q_I/(p_network(ParentNode)-p_darcy(i));
end
end