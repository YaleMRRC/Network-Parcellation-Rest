function [Feature_imp_sum,mot_female,mot_male] = feature_importance_analysis(feature_imp,Mot_825,gender_825)

    % feature_imp = a 268x1 (or 188x1) vector including the importance
    % score of each feature (i.e. node) in the gradient boosting machine
    % sex classification. It is the output of the python code for sex
    % prediction.
    
    % Mot_825 = head motion values for all subejcts
    % gender_825 = sex labels for all subjects
    
    % Feature_imp_sum = the normalized feature importances across all
    % number of networks
    % mot_female and mote_male is the heaadmotion values for female and
    % male subjects respectively

    feature_imp(1,:)=[]
    A = feature_imp';
    A_max = max(sum(A,2));


    % close all
    for i=2:25
        C{i-1} = num2str(i);
    end
    figure('color',[1,1,1])


    xlim ([0,length(cortical)+1])
    legend(C,'FontWeight','bold','Location','northeastoutside')

    Feature_imp_sum=zeros(1,188);

    for K = 1:24
        Feature_imp_sum=Feature_imp_sum+feature_imp(K,:);
    end

    colormap(parula)
    [mm,nn] = sort(Feature_imp_sum,'descend')    
    temp = (A./A_max)*100;
    bar(temp(nn,:), 'stacked', 'LineWidth',2)
    
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 40)
    ax = gca;
    ax.YAxis.Exponent = 0
    set(gca,'xtick',[1,length(cortical)],'xticklabel',[1,length(cortical)])

    xlabel ('Node','FontSize',60,'FontWeight','bold');
    ylabel('Feature importance','FontSize',60,'FontWeight','bold');
    Feature_imp_sum = Feature_imp_sum/max(Feature_imp_sum)*100;
    max(Feature_imp_sum)
    xlim ([0,189])
    ylim([0,101])


    % Motion differences between male and female
       
%     load /Users/Mehraveh/Documents/MATLAB/Parcellation/October3/Mot_825.mat
%     load /Users/Mehraveh/Documents/MATLAB/Parcellation/October3/gender_825.mat
    mot = (Mot{1,1}+Mot{2,2}) / 2;
    mot_female=mot(gender(:,2)==1)'
    mot_male=mot(gender(:,2)==0)'

    figure,
    C = [mot_female; mot_male];
    grp = [zeros(1,length(mot_female)),ones(1,length(mot_male))];
    boxplot(C,grp)
    [h,p,ci,stats]=ttest2(mot_female,mot_male)
    set(gca,'xtick',1:2,'xticklabel',{'Female','Male'},'FontSize',25,'FontWeight','bold'); 
    format short e
    h = findobj(gca,'Type','line')
    set(h,'LineWidth',3)
    title({'Motion effect for female vs. male'; ['(p < ', num2str(p,'%.2e'),')']}, 'FontSize', 30)
end
