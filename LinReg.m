function[a,b]=LinReg(kT,r,Fig)
NotZeroVals = zeros(size(kT));
NotZeroVals(kT>0)=1;

%%% KUN NETTVERK %%%
rToPlot=[];
ktToPlot=[];
for iter = 1:size(kT,2)
    rToPlot = [rToPlot;r(NotZeroVals(:,iter)==1,iter)];
    ktToPlot = [ktToPlot;kT(NotZeroVals(:,iter)==1,iter)];
end
linreg = fitlm(rToPlot,ktToPlot);
a = linreg.Coefficients.Estimate(2);
b = linreg.Coefficients.Estimate(1);
if strcmp(Fig,'Figure')
    figure()
    plot(linreg);
    %set(gca,'YScale','log');
    xlabel('R','FontSize',18)
    ylabel('K^{T}','Interpreter','tex','Rotation',0,'FontSize',18);
    title(' ')
    set(legend,'visible','off')
end
