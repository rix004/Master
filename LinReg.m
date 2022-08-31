function[a,b]=LinReg(kT,r,iterations,MicroTermIndexes,MacroTermIndexes,TrueKT)
NotZeroVals = zeros(size(kT));
NotZeroVals(kT>0)=1;

%%% KUN NETTVERK %%%
rToPlot=[];
ktToPlot=[];
for iter = 1:iterations
    rToPlot = [rToPlot;r(NotZeroVals(:,iter)==1,iter)];
    ktToPlot = [ktToPlot;kT(NotZeroVals(:,iter)==1,iter)];
end


linreg = fitlm(rToPlot,log(ktToPlot));
a = linreg.Coefficients.Estimate(2);
b = linreg.Coefficients.Estimate(1);

figure()
plot(linreg);
set(gca,'YScale','log');
xlabel('R','FontSize',18)
ylabel('ln kT','Rotation',0,'FontSize',18);
title(' ')
set(legend,'visible','off')

%%% HELE SYSTEMET %%%%

% figure('Name','kT vs R')
% rToPlot=[];
% ktToPlot=[];
% aToPlot=[];
% for iter = 1:iterations
%     rToPlot = [rToPlot;r(T_cells1(:,iter)==1,iter)];
%     ktToPlot = [ktToPlot;kT1(T_cells1(:,iter)==1,iter)];
%     aToPlot = [aToPlot;A(T_cells1(:,iter)==1,iter)];
% end
% plot(rToPlot,log(ktToPlot),'.','MarkerSize',20);
% set(gca,'YScale','log');
% linreg = fitlm(rToPlot,log(ktToPlot),'Weights',aToPlot);
% plot(linreg);
% xlabel('|x^{N}_i - x^{D}|')
% ylabel('ln(K^{T}_i)')
 
% figure
% for i = 1:length(rToPlot)
%     plot(rToPlot(i),ktToPlot(i),'.','MarkerSize',aToPlot(i)*3500)
%     hold on
% end
% set(gca,'YScale','log');
% xlabel('|x^{N}_i - x^{D}|','FontSize',15)
% ylabel('K^{T}_i','FontSize',15)
% title('Radius rate = ',num2str(RandomTree.RadiusRate));
% mean_kt = zeros(length(MicroTermIndexes),1);
% mean_r = zeros(length(MicroTermIndexes),1);
% for i = 1:length(MicroTermIndexes)
%     mean_kt(i)=mean(kT(i,T_cells(i,:)==1));
%     mean_r(i)=mean(d(i,T_cells(i,:)==1));
% end

% std_dev = std(kT');
% errorbar(mean_r,mean_kt,std_dev,'x-','LineWidth',0.75,'Color',[0   0   0])
% hold on
% plot(mean_r(mean_r>0),mean_kt(mean_r>0),'-','LineWidth',3,'Color',[0 0 0])
% hold on
end