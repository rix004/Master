function[p_darcy,q_network,p_network]=SolveSystemEx(Tree,TNinfo,TNlogic,Dir_network,Neu_network,mu,k,K_N,LHS,RHS,cell_area,BC)
nodes = Tree.nodes;
edges = Tree.edges;
Ne = length(edges);
Nn = length(nodes);
Nin = sum(1-TNlogic)-1;                        % Interior nodes
Npn = sum(1-TNlogic);                          % Parent nodes
Ntn = size(TNinfo,1);
MicroTermIndexes = find(TNlogic==1);
connections = sparse([1:Ne, 1:Ne]', [edges(:,2); edges(:,3)], [-ones(Ne,1);ones(Ne,1)]);
I = spdiags(ones(Ntn,1),0,Ntn,Ntn);

I_term = zeros(Ntn,Nn);
for i = 1:size(I_term,1)
    I_term(i,MicroTermIndexes(i))=1;
end

if strcmp(BC, 'Dirichlet')
    Bc_nodes = zeros(Nn,1);
    Bc_nodes(Tree.RootNodeIdx)=1;
    connections_trimmed = connections(:,(TNlogic+Bc_nodes)==0);

    PC = zeros(Ntn,Ne);
    for i = 1:size(PC,1)
        PC(i,TNinfo(i,3))= mu/(2*k*pi)*log(0.2*sqrt(cell_area(i))/edges(TNinfo(i,3),4));
        
    end
    PC = sparse(PC);

    TE = zeros(Ntn,Ne);
    for i = 1:size(TE,1)
        TE(i,TNinfo(i,3))=1;
    end
    TE = sparse(TE);

    I_term_trimmed = I_term(:,Bc_nodes==0);
    TN = I_term(:,TNlogic==1);

    A = [LHS, -TE, sparse(Ntn,Nn-1);
    sparse(Ne,Ntn), spdiags(K_N.^-1,0,Ne,Ne), connections(:,Bc_nodes==0);
    sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Nn-1);
    TN, PC, -I_term_trimmed];
    
    rhs = [RHS; -connections*Dir_network*Bc_nodes; zeros(Nin,1); zeros(Ntn,1)];
    SOL = A\rhs;

    p_darcy = SOL(1:Ntn);
    q_network = SOL(Ntn+1:Ntn+Ne);
    p_network_temp = SOL(Ntn+Ne+1:end);
    p_int = SOL(Ntn+Ne+1:Ntn+Ne+Nin);
    p_tn = SOL(Ntn+Ne+Nin+1:Ntn+Ne+Nin+Ntn);

    p_network = zeros(Nn,1);
    count1 = 0;
    for i = 1:length(p_network)
        if Bc_nodes(i)==1
            p_network(i)=Dir_network;
        else
            count1 = count1 +1;
            p_network(i)=p_network_temp(count1);
        end
    end

elseif strcmp(BC, 'Neumann')
    Bc_edges = [1;zeros(Ne-1,1)];
    connections_trimmed = connections(Bc_edges==0,(TNlogic+[1;zeros(Nn-1,1)])==0);
    connections_trimmed1 = connections(:,TNlogic==0);

    PC = zeros(Ntn,Ne-1);
    for i = 1:size(PC,1)
        PC(i,TNinfo(i,3)-1)=mu/(2*k*pi)*log(sqrt(cell_area(i))*0.2/edges(TNinfo(i,3),4));
    end
    PC = sparse(PC);

    PM_vect = zeros(Ntn,1);
    for i = 1:size(PM_vect,1)
        PM_vect(i)=mu/(2*k*pi)*log(0.2*sqrt(cell_area(i))/edges(TNinfo(i,3),4));
    end
    PM = spdiags(PM_vect,0,Ntn,Ntn);
    
    TE= zeros(Ntn,Ne-1);
    for i = 1:size(TE,1)
        TE(i,TNinfo(i,3)-1)=1;
    end
    TE = sparse(TE);

    I_q = zeros(Ne,Ne-1);
    I_q(2:end,1:end)=spdiags(ones(Ne-1,1),0,Ne-1,Ne-1);

    A = [LHS, -TE, sparse(Ntn,Npn), sparse(Ntn,Ntn);
    sparse(Ne,Ntn), I_q, K_N.*connections;
    sparse(Nin,Ntn), connections_trimmed', sparse(Nin,Npn), sparse(Nin,Ntn);
    I, PC, -I_term];
        
    rhs = [RHS; -Neu_network*Bc_edges; [-Neu_network;zeros(Nin-1,1)]; zeros(Ntn,1)];

    SOL = A\rhs;
    
    p_darcy = SOL(1:Ntn);
    q_network = [Neu_network;SOL(Ntn+1:Ntn+Ne-1)];
    p_pn = SOL(Ntn+Ne:Ntn+Ne+Nin);
    p_tn = SOL(Ntn+Ne+Nin+1:Ntn+Ne+Nin+Ntn);
    p_network = [p_pn;p_tn];
end