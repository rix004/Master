function[p_darcy,q_T,q_network,p_network]=SolveSystemUpS(Tree,TNlogic,TNinfo,Dir_network,Neu_network,mu,k,K_N,K_T,Term2Cell,LHS,RHS,cell_area,BC)
nodes = Tree.nodes;
edges = Tree.edges;
Ne = size(edges,1);
Nn = size(nodes,1);
Nin = sum(1-TNlogic)-1;                        % Interior nodes
Npn = sum(1-TNlogic);                          % Parent nodes
Ntn = sum(TNlogic);
Ncells = length(cell_area);
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]);
I = spdiags(ones(Ncells,1),0,Ncells,Ncells);
Volume=I.*cell_area;

if strcmp(BC, 'Dirichlet')
    Bc_nodes = zeros(Nn,1);
    Bc_nodes(Tree.RootNodeIdx)=1;
    %Term_edges = zeros(Ne,1);
    %Term_edges(TNinfo(:,3))=1;
    connections_trimmed = connections(:,(TNlogic+Bc_nodes)==0);

    I2 = zeros(Ne,Ntn);
    for i = 1:size(I2,2)
        I2(TNinfo(i,3),i)=1;
    end
    
    Term2Cell = Term2Cell.*cell_area;

    A = [LHS, sparse(Ncells,Ne), sparse(Ncells,Nin), sparse(Ncells,Ntn),-Volume;                        % Darcy´s lov + massebevaring i Darcydomenet (2.1 og 2.4)
    sparse(Ne,Ncells), spdiags(K_N.^-1,0,Ne,Ne), connections_trimmed, I2, sparse(Ne,Ncells);            % Poiseuilles lov (2.5)
    sparse(Nin,Ncells), connections_trimmed', sparse(Nin,Nin), sparse(Nin,Ntn), sparse(Nin,Ncells);     % Massebevaring i indre noder (2.2)
    sparse(Ntn,Ncells), I2', sparse(Ntn,Nin), sparse(Ntn,Ntn), -Term2Cell;                              % Massebevaring i terminalnoder (2.3)
    spdiags(nonzeros(K_T),0,Ncells,Ncells), sparse(Ncells,Ne), sparse(Ncells,Nin), -K_T, I];            % Flyt fra nettverk til domene (2.6)
    
    rhs = [RHS; -connections*Dir_network*Bc_nodes; zeros(Nin,1); zeros(Ntn,1); zeros(Ncells,1)];

    SOL = A\rhs;

    p_darcy = SOL(1:Ncells);
    q_network = SOL(Ncells+1:Ncells+Ne);
    p_int = SOL(Ncells+Ne+1:Ncells+Ne+Nin);
    p_tn = SOL(Ncells+Ne+Nin+1:Ncells+Ne+Nin+Ntn);
    q_T = SOL(Ncells+Ne+Nin+Ntn+1:end);

    int_indexes = find(Bc_nodes+TNlogic==0);
    p_network = Dir_network*Bc_nodes;
    p_network(int_indexes) = p_int;
    p_network(find(TNlogic==1)) = p_tn;

elseif strcmp(BC, 'Neumann')
    Bc_edges = [1;zeros(Ne-1,1)];
    connections_trimmed = connections(Bc_edges==0,(TNlogic+[1;zeros(Nn-1,1)])==0);
    connections_trimmed1 = connections(:,TNlogic==0);
    
    I_term = zeros(Ntn,Nn);
    for i = 1:size(I_term,1)
        I_term(i,MicroTermIndexes(i))=1;
    end

    P = zeros(Ntn,Ne-1);
    for i = 1:size(P,1)
        P(i,MicroTN(i,3)-1)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(MicroTN(i,3),4));
    end
    P = sparse(P);

    PM_vect = zeros(Ntn,1);
    for i = 1:size(PM_vect,1)
        PM_vect(i)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(MicroTN(i,3),4));
    end
    PM = spdiags(PM_vect,0,Ntn,Ntn);
    
    C= zeros(Ntn,Ne-1);
    for i = 1:size(C,1)
        C(i,MicroTN(i,3)-1)=1;
    end
    C = sparse(C);

    I_q = zeros(Ne,Ne-1);
    I_q(2:end,1:end)=spdiags(ones(Ne-1,1),0,Ne-1,Ne-1);

    A = [LHS, sparse(Ntn,Ne-1), sparse(Ntn,Npn), sparse(Ntn,Ntn), -I;
        sparse(Ne,Ntn), I_q, K_N.*connections, sparse(Ne,Ntn);
        sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Npn), sparse(Nin,Ntn), sparse(Nin,Ntn);
        sparse(Ntn,Ntn), C, sparse(Ntn,Npn), sparse(Ntn,Ntn), -I;
        I, sparse(Ntn,Ne-1), -I_term, PM];
    rhs = [RHS; -Neu_network*Bc_edges; [-Neu_network;zeros(Nin-1,1)]; zeros(Ntn,1); zeros(Ntn,1)];


    SOL = A\rhs;
    
    p_darcy = SOL(1:Ntn);
    q_network = [Neu_network;SOL(Ntn+1:Ntn+Ne-1)];
    p_pn = SOL(Ntn+Ne:Ntn+Ne+Nin);
    p_tn = SOL(Ntn+Ne+Nin+1:Ntn+Ne+Nin+Ntn);
    q_T = SOL(Ntn+Ne+Nin+Ntn+1:end);
    p_network = [p_pn;p_tn];
end