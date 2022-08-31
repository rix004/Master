function[q,p_all]=SolveNetwork(Tree,kN,TermNodes,bc)
edges = Tree.edges;
nodes = Tree.nodes;
Ne = size(edges,1);
Nn = size(nodes,1);
Bc_nodes = TermNodes;  % Noder som du gir randverdier i
Bc_nodes(Tree.RootNodeIdx)=1;
Bc_vals = bc*Bc_nodes;  % Randverdier
Bc_vals(Tree.RootNodeIdx)=1;
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]); % Forutsetter at kolonne 2 og 3 i «edges» angir node-nummer
connections_trimmed = connections(:,Bc_nodes==0);
A = [spdiags(kN.^-1,0,Ne,Ne), connections_trimmed; connections_trimmed', sparse(size(connections_trimmed,2),size(connections_trimmed,2))];
rhs = [-connections*Bc_vals; zeros(sum(1-Bc_nodes),1)];
sol = A\rhs;
q = sol(1:Ne);
p = sol((Ne+1):end);
p_all = Bc_vals;
p_all(Bc_nodes==0) = p;
end