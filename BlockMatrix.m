function[A,rhs]=BlockMatrix(cells,edges,vertices,Tn,Te,f,k,g,Bc_nodes,Bc_vals,Bc,Bc_val,Ne,Ntn,Nin)
[LHS,~,cell_center,cell_edges,cell_area] = TPFA(cells,vertices,f,k,Bc,Bc_val);
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]);
connections_trimmed = connections(:,Bc_nodes==0);

I = spdiags(ones(Ntn,1),0,Ntn,Ntn);
q_term = zeros(Ntn,1);
q_term(Te)=1;
Q_term = spdiags(q_term,0,Ntn,Ne);

B = zeros(Ntn,Ne);
for i = 1:size(B,1)
    B(i,Tn(i,3))=1/(2*pi*k(cell_center(i,1),cell_center(i,2)))*log(sqrt(cell_area(i))*0.2/edges(Tn(i,3),4));
end
B = sparse(B);

C = zeros(Ntn,Ne);
for i = 1:size(C,1)
    C(i,Tn(i,3))=1;
end
C = sparse(C);

A = [LHS, sparse(Ntn,Ne), sparse(Ntn,Nin), -I, sparse(Ntn,Ntn);
    sparse(Ne,Ntn), spdiags(g.^-1,0,Ne,Ne), connections_trimmed, sparse(Ne,Ntn), sparse(Ne,Ntn);
    sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Nin), sparse(Nin,Ntn), sparse(Nin,Ntn);
    sparse(Ntn,Ntn), C, sparse(Ntn,Nin) -I, sparse(Ntn,Ntn);
    I, -B, sparse(Ntn,Nin), sparse(Ntn,Ntn), -I];
rhs = [zeros(size(LHS,1),1); -connections*Bc_vals; zeros(Nin,1); zeros(Ntn,1); zeros(Ntn,1)];
end