function[K_T,D]=findKT1(edges,cell_area,MicroTermIndexes,MacroTermIndexes,MicroTN,q_network,p_network,mu)
G = graph(edges(:,2)',edges(:,3)',edges(:,1)');
d=distances(G);
dTerm = d(MicroTermIndexes,MacroTermIndexes);
K_T = zeros(length(MicroTermIndexes),length(MacroTermIndexes));
D = zeros(length(MicroTermIndexes),length(MacroTermIndexes));
for MicroIndex = 1:size(dTerm,1)
    MinDist = min(dTerm(MicroIndex,:));
    MacroTindex = find(dTerm(MicroIndex,:)==MinDist);
    ParentNode = MacroTermIndexes(MacroTindex);
    q_I = q_network(MicroTN(MicroIndex,3));
    K_T(MicroIndex,MacroTindex) = q_I/((p_network(ParentNode)-p_network(MicroTermIndexes(MicroIndex)))*cell_area(MicroIndex));
    D(MicroIndex,MacroTindex) = 1;
end
D = D';
end