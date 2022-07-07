function[r]=TotalResistance(MacroNode,MicroNode,edges,K_N)
e = find(edges(:,3)==MicroNode);
r = K_N(e);
n2 = edges(e,2);
while n2 > MacroNode
    e = find(edges(:,3)==n2);
    r = r + K_N(e);
    n2 = edges(e,2);
end
end