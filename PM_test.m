clc;
clear;
Num_tests = 2;
Num_doubles = 7;
for iter1 = 1:Num_tests
    for iter2 = 1:Num_doubles
        Np = 25*2^(iter2);
        [error,cell_size]=Calculations(Np);
        L2_error(iter1,iter2) = error;
        h(iter1,iter2) = cell_size;
    end
end

figure
% loglog([0.25 0.25/2 0.25/4 0.25/8],[0.1 0.05 0.025 0.025/2],'-','LineWidth',2.5,'Color','r')
% hold on

for iter3 = 1:size(L2_error,1)
    h(iter3,:)
    loglog(h(iter3,:),L2_error(iter3,:),'-','LineWidth',1)
    hold on
end

% mean_vect = zeros(1,size(L2_error,2));
% mean_h = zeros(1,size(L2_error,2));
% for iter4 = 1:size(L2_error,2)
%     mean_vect(iter4)=mean(L2_error(:,iter4));
%     mean_h(iter4)=mean(h(:,iter4));
% end
% 
% std_dev = std(L2_error);
% errorbar(mean_h,mean_vect,std_dev,'x-','LineWidth',3,'Color',[0 0 0])
% hold on
% % loglog(mean_h,mean_vect,'-','LineWidth',2.5,'Color',[0 0 0])
% % hold on
% 
% xlabel('(1/n_{cells})^{1/2}','FontSize',18)
% ylabel('Error','FontSize',17)
% lgd = legend('1st order reference');
% lgd.FontSize=(14);
% legend('Location','northwest')