function Hamming_dist_25_25_825_allreazliations = hamming_distance_analysis(M,HCP_subj,label_134_cort,K_max)

    % M is an output of the function "import_HCP_data.m" containing the
    % rest scan data with the following format:
    % M{1,1} = REST_LR, M{1,2} = REST_RL, M{2,1} = REST2_LR, M{2,2} = REST2_RL
    
    % HCP_subj is an output of the function "import_HCP_data.m" containing the
    % subject names with the following format:
    % HCP_subj{1,1} = REST_LR, HCP_subj{1,2} = REST_RL, HCP_subj{2,1} = REST2_LR, HCP_subj{2,2} = REST2_RL
     
    % label_134_cort = a 268x1 vector including node labels according to
    % their position: 1: Subcortical, 3: Cerebellum, 4: Cortical
    % K_max = The number of networks, e.g. 25
    
    % Hamming_dist_25_25_825_allreazliations = the hamming distance
    % analysis according to the pipeline described in the paper
    
    
    clc

    subjsLR = intersect(HCP_subj{1,1},HCP_subj{2,1});
    subjsRL = intersect(HCP_subj{1,2},HCP_subj{2,2});
    subjs_all = intersect(subjsLR,subjsRL);

    LR1 = arrayfun(@(x)find(HCP_subj{1,1}==x,1),subjs_all);
    LR2 = arrayfun(@(x)find(HCP_subj{2,1}==x,1),subjs_all);
    RL1 = arrayfun(@(x)find(HCP_subj{1,2}==x,1),subjs_all);
    RL2 = arrayfun(@(x)find(HCP_subj{2,2}==x,1),subjs_all);


    % load /home/mehraveh/documents/MATLAB/Parcellation/label_134_cort.mat  % 1:cortical 3:cerebellum, 4:cortical
    cortical_boolean = input('Do you want to consider only the cortical regions? (1: cortical, 0: whole brain)')
    if cortical_boolean == 1
        vect = [3];
    elseif cortical_boolean == 0
        vect = [1,3,4];
    end
    cortical = label_134_cort(ismember(label_134_cort(:,2),vect),1);

    l_full = length(LR1)

    for subj=1:l_full
        subj
        V = [M{1,1}{LR1(subj)}(:,cortical);M{2,1}{LR2(subj)}(:,cortical);M{1,2}{RL1(subj)}(:,cortical);M{2,2}{RL2(subj)}(:,cortical)];

        t = size(V,1);
        n = size(V,2);

        twoNorm = sqrt(sum(abs(V).^2,1));
        m = max(twoNorm);
        V = V/m;

        sqDistances_HCP{subj} = sqDistance(V);               
        D = sqDistances_HCP{subj};        
        e0 = zeros(t,1);
        e0(1)=e0(1)+3;
        d0 = sqDistance_Y(V,e0);

        if min(d0)<=max(max(D(1:n,1:n)))
            fprintf('No :( \n')
        end

        d0_HCP{subj} = d0;
    end    



    for numberofruns = 1:100
        numberofruns        
        l_full_shuffled = randperm(l_full);
        l_full_shuffled_100{numberofruns} = l_full_shuffled;
        f = zeros(n,l_full);
        S_opt_all = [];

        K = K_max;
        for involved_subj = 25:25:l_full  % The number of subjects in all tasks > 800
            temp_sqDsitance = sqDistances_HCP(l_full_shuffled(1:involved_subj));
            temp_d0 = d0_HCP(l_full_shuffled(1:involved_subj));
            S_opt = exemplar(temp_sqDsitance,temp_d0,t,n,K,involved_subj);
            S_opt_all{K,involved_subj}=S_opt;
        end
        temp = S_opt_all;

        for K=K_max:-1:2
            for involved_subj = 25:25:825 % The number of subjects in all tasks > 800
                S_opt_all{K,involved_subj} = temp{K_max,involved_subj};
                temp{K_max,involved_subj}=temp{K_max,involved_subj}(1:end-1);
            end
        end

        S_opt_all_100{numberofruns} = S_opt_all;

        inv=0;
        clear hamm_dist
        for K=2:K_max
            l=length(sqDistances_HCP);
            inv=0;
            for involved_subj = 25:25:825   % The number of subjects in all tasks > 800
                ind = S_opt_all{K,involved_subj};
                for subj=l_full_shuffled
                    D = sqDistances_HCP{subj};                
                    [D_sorted,index_global{K,involved_subj}(subj,:)] = min(D(ind,1:end),[],1);                          
                end
                inv=inv+1;
                [Maj,F] = mode(index_global{K,involved_subj}(l_full_shuffled(1:involved_subj),:));
                [Maj_full,F_full] = mode(index_global{K,involved_subj});
                hamm_dist(K,inv) = pdist2(Maj,Maj_full,'hamming')*length(cortical);
            end
        end
        Hamming_dist_25_25_825_allreazliations{numberofruns} = hamm_dist;
    end


    close all
    figure
    even = [1:12]*2;
    for i=even
        C{i/2} = num2str(i);
    end

    cc=jet(K_max);
    clear hamm
    for num=1:100
        for K=2:K_max
            hamm{K}(:,num) = Hamming_dist_25_25_825_allreazliations{num}(K,:);
        end
    end
    for num=1:100
        for K=even
            y_mean(K,:) = mean(hamm{K},2);
            err(K,:) = std(hamm{K},[],2);
        end
    end
    x = 25:25:825
    for K=even
        subplot(3,4,K/2)
        h = errorbar(x,y_mean(K,:),err(K,:),  '-s','MarkerSize',2.5,...
            'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',1, 'markersize',2.5); hold on
        xt = set(gca, 'Xtick',25:200:825,'xticklabel',{'25';'225';'425';'625';'825'});
        set(gca, 'FontSize', 17)
        set(gca, 'ylim', [0,35])
        set(gca, 'xlim', [0,850])   
        title(['K = ', num2str(K)],'FontSize',20,'FontWeight','bold');

    end

    axes( 'Position', [0.5, 0, 1, 1] ) ;
    set( gca, 'Color', 'None', 'XColor', 'none', 'YColor', 'none' ) ;
    text( 0, 0, 'Number of subjects', 'FontSize', 30', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
    axes( 'Position', [0.09, 0.5, 1, 1] ) ;
    set( gca, 'Color', 'None', 'XColor', 'none', 'YColor', 'none') ;
    h = text( 0, 0, 'Hamming distance', 'FontSize', 30', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
    set(h, 'rotation', 90)
end
