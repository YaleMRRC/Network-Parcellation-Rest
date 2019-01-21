function  [Inv_F_1st_sum, Inv_F_ratio_sum] = variability_analysis(index_global,K_max)
    
    % index_global is the individualized parcellation schemes for all subjects.
    % index_global is an output of the function "Parcellation_rest.m"
    % K_max = The number of networks, e.g. 25
    
    % Inv_F_1st_sum = 1/F1 
    % Inv_F_ratio_sum = F2/F1 

    
%     load /Users/Mehraveh/Documents/MATLAB/Parcellation/October3/BrainFigures_8_9_concat_825

    clear Maj

    l_cortical = size(index_global{2},2);

    for K=2:K_max
        [Maj(K,:), F_1st(K,:)] = mode(index_global{K});
    end
    for K=2:K_max
        for node = 1:l_cortical
            freq{K,node} = hist(index_global{K}(:,node),1:K);
        end
    end


    for K=2:K_max
        for node = 1:l_cortical
            [a, ind] = sort(freq{K,node},'descend');
            Maj_1{K}(node)=ind(1);
            Maj_2{K}(node)=ind(2);
            F_1{K}(node)=a(1);
            F_2{K}(node)=a(2);

        end
    end

    clear F_ratio F_1st
    for K=2:K_max

        F_1st(K,:) = F_1{K};
        F_2_temp = F_2{K};
        F_2_temp(F_2_temp == 0) = 1;
        F_ratio(K,:) = F_1{K}./F_2_temp;
        Inv_F_1st(K,:) = 1./F_1{K};
        Inv_F_ratio(K,:) = F_2{K}./F_1{K};
    end

    F_ratio(1,:)=[];
    F_1st(1,:)=[];
    Inv_F_1st(1,:)=[];
    Inv_F_ratio(1,:)=[];

    close all
    figure
    cc=jet(K_max);
    F_1st_sum=zeros(1,188);
    F_ratio_sum=zeros(1,188);
    Inv_F_1st_sum=zeros(1,188);
    Inv_F_ratio_sum=zeros(1,188);

    for K = 1:K_max-1
    F_1st_sum=F_1st_sum+F_1st(K,:);
    F_ratio_sum=F_ratio_sum+F_ratio(K,:);
    Inv_F_1st_sum=Inv_F_1st_sum+Inv_F_1st(K,:);
    Inv_F_ratio_sum=Inv_F_ratio_sum+Inv_F_ratio(K,:);
    end

    settt = {Inv_F_1st';Inv_F_ratio'};
    settt_sum = {Inv_F_1st_sum';Inv_F_ratio_sum'};
    settt_name = {'Inverse F1','F2:F1 ratio'};
    settt_y = {'Variability (1/F1)','Variability (F2/F1)'};

    for num=1:2
    A = settt{num}

    A_max= max(settt_sum{num});
    A_min= min(settt_sum{num})
    for i=2:K_max

        C{K_max-i+1} = num2str(i);
    end
    figure('Color',[1 1 1]);

    colormap(parula)

    [mm,nn] = sort(settt_sum{num},'descend')
    temp = (A./A_max*100);
    bar(temp(nn,:), 'stacked', 'LineWidth',2)

    set(gca, 'FontSize', 10)
    ax = gca;
    ax.YAxis.Exponent = 0


    xlabel ('Node','FontSize',60,'FontWeight','bold');
    ylabel(settt_y{num},'FontSize',60,'FontWeight','bold');
    xlim ([0,l_cortical+1])
    ylim([0,101])
    set(gca,'xtick',[1,188],'xticklabel',[1,l_cortical])
    axis('square'   )
    lgd = legend(C,'FontWeight','bold','Location','eastoutside','box','off','Orientation','vertical')

    end

    F_1st_sum = F_1st_sum/max(F_1st_sum)*100;
    max(F_1st_sum)
    F_ratio_sum = F_ratio_sum/max(F_ratio_sum)*100;
    max(F_ratio_sum)
    Inv_F_1st_sum = Inv_F_1st_sum/max(Inv_F_1st_sum)*100;
    max(Inv_F_1st_sum)
    Inv_F_ratio_sum = Inv_F_ratio_sum/max(Inv_F_ratio_sum)*100;
    max(Inv_F_ratio_sum)
end