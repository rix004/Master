clc;clear;

% Antall nivåer på treet, noder og kanter
num_levels = 4;
num_nodes = 2^num_levels-1;
num_edges = num_nodes -1;
nodes = [1:num_nodes];
edges = [1:num_edges];
k = 1;

% levels: Oversikt over hvilke noder som er på hvert nivå. 
levels = zeros(4,2^(num_levels-1));
for l = 1:num_levels
    levels(l,1:2^(l-1)) = nodes(2^(l-1):2^l-1);
end


% Matrise som angir forbindelser mellom noder og kanter. 1 angir en
% forbindelse på oversiden av kanten, -1 en forbindelse på undersiden, 0 angir ingen forbindelse
connections = zeros(num_nodes,num_edges); 
for i = 2:num_nodes
    connections(i, i-1) = -1;
end
k=1;
for i = 1:num_nodes
    if k <= num_edges
        connections(i,k)=1;
        connections(i, k+1)=1;
        k = k+2;
    elseif k > num_edges
        continue;
    end
end

% Matrise som angir posisjonen til nodene
position = zeros(num_nodes,2); x_pos = 0; y_pos=0;
% for j = 1:nodes
%     xPosition(j) = x_start + (-1)^(j-1) * (1/j);
%     for i = 1:levels
%         if j <= 2^i
%             yPosition(j)=y_pos-1;
%         else
%             
%     end
% end


% k = 0;
% m = 0;
% for i = 1:levels
%     for j = 1:nodes
%         if j > 2^(i-1)
%             continue
%         else
%             m = m+1;
%             xPosition(m) = x_pos + (1/m);
%             yPosition(m) = y_pos-i;
%         end
%     end
%     x_pos = x_pos - 1;
% end

% comb_position(:,1)=xPosition';
% comb_position(:,2)=yPosition';
position(1,1)=x_pos; position(1,2)=y_pos;
position(2,1)=x_pos-1; position(2,2)=y_pos-1;
position(3,1)=x_pos+1; position(3,2)=y_pos-1;
position(4,1)=x_pos-1-(1/2); position(4,2)=y_pos-2;
position(5,1)=x_pos-1+(1/2); position(5,2)=y_pos-2;
position(6,1)=x_pos+1-(1/2); position(6,2)=y_pos-2;
position(7,1)=x_pos+1+(1/2); position(7,2)=y_pos-2;
position(8,1)=x_pos-1-(1/2)-(1/4); position(8,2)=y_pos-3;
position(9,1)=x_pos-1-(1/2)+(1/4); position(9,2)=y_pos-3;
position(10,1)=x_pos-1+(1/2)-(1/4); position(10,2)=y_pos-3;
position(11,1)=x_pos-1+(1/2)+(1/4); position(11,2)=y_pos-3;
position(12,1)=x_pos+1-(1/2)-(1/4); position(12,2)=y_pos-3;
position(13,1)=x_pos+1-(1/2)+(1/4); position(13,2)=y_pos-3;
position(14,1)=x_pos+1+(1/2)-(1/4); position(14,2)=y_pos-3;
position(15,1)=x_pos+1+(1/2)+(1/4); position(15,2)=y_pos-3;

plot(position(:,1),position(:,2),'*')
