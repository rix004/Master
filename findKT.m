function[K_T,D]=findKT(edges,cells,cell_area,MicroTermIndexes,MacroTermIndexes,MicroTN,q_network,p_network,p_darcy,mu)
G = graph(edges(:,2)',edges(:,3)',edges(:,1)');
d=distances(G);
dTerm = d(MicroTermIndexes,MacroTermIndexes);
K_T = zeros(size(cells,1),length(MacroTermIndexes));
for MicroIndex = 1:size(dTerm,1)
    MinDist = min(dTerm(MicroIndex,:));
    MacroTindex = find(dTerm(MicroIndex,:)==MinDist);
    ParentNode = MacroTermIndexes(MacroTindex);
    q_I = q_network(MicroTN(MicroIndex,3));
    K_T(MicroIndex,MacroTindex) = q_I*mu/((p_network(ParentNode)-p_darcy(MicroIndex))*cell_area(MicroIndex));
    D(MicroIndex,MacroTindex) = MinDist;
end
end