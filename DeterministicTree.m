% Deterministic tree data
DT.Levels = 4;
DT.StartPos = [0 0];
DT.StartAngle = 90;
DT.RotationAngle = 42.25;
DT.TrunkRadius = 0.1;    % mm
DT.RadiusRate = 0.7;
DT.TrunkLength = 0.7; %mm
DT.LengthRate = 0.65;

Tree= GetTree(DT);
DrawTree(Tree,10,[0.8500, 0.3250, 0.0980]);


nodes = Tree.nodes(:,1:2);
edges = Tree.edges;

nodes = FindLevels(nodes,edges);