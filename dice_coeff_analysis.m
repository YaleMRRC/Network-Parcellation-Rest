function reproducibility = dice_coeff_analysis(M,HCP_subj,label_134_cort,K_max)

    % M is an output of the function "import_HCP_data.m" containing the
    % rest scan data with the following format:
    % M{1,1} = REST_LR, M{1,2} = REST_RL, M{2,1} = REST2_LR, M{2,2} = REST2_RL
    
    % HCP_subj is an output of the function "import_HCP_data.m" containing the
    % subject names with the following format:
    % HCP_subj{1,1} = REST_LR, HCP_subj{1,2} = REST_RL, HCP_subj{2,1} = REST2_LR, HCP_subj{2,2} = REST2_RL
     
    % label_134_cort = a 268x1 vector including node labels according to
    % their position: 1: Subcortical, 3: Cerebellum, 4: Cortical
    % K_max = The number of networks, e.g. 25
        
    % reproducibility = dice coefficients between two halves of data, 
    % calculated according to the pipeline described in the paper

    % Parcellation with two Halves of 400 - Rest8+Rest9+LR+RL concatenated all)
    subjsLR = intersect(HCP_subj{1,1},HCP_subj{2,1})
    subjsRL = intersect(HCP_subj{1,2},HCP_subj{2,2})
    subjs_all = intersect(subjsLR,subjsRL)

    LR1 = arrayfun(@(x)find(HCP_subj{1,1}==x,1),subjs_all)
    LR2 = arrayfun(@(x)find(HCP_subj{2,1}==x,1),subjs_all)
    RL1 = arrayfun(@(x)find(HCP_subj{1,2}==x,1),subjs_all)
    RL2 = arrayfun(@(x)find(HCP_subj{2,2}==x,1),subjs_all)

    S_opt = cell(2,1)
    load /home/mehraveh/documents/MATLAB/Parcellation/label_134_cort.mat  % 1:cortical 3:cerebellum, 4:cortical

    cortical_boolean = input('Do you want to consider only the cortical regions? (1: cortical, 0: whole brain)')
    if cortical_boolean == 1
        vect = [3];
    elseif cortical_boolean == 0
        vect = [1,3,4];
    end
    cortical = label_134_cort(ismember(label_134_cort(:,2),vect),1);

    % for this part we focus the analsyis on the first 800 subjcts
    l_full = 800; 
    l_half = 400;

    for subj=1:l_full
        V = [M{1,1}{LR1(subj)}(:,cortical);M{2,1}{LR2(subj)}(:,cortical);M{1,2}{RL1(subj)}(:,cortical);M{2,2}{RL2(subj)}(:,cortical)];

        t = size(V,1);
        n = size(V,2);

        twoNorm = sqrt(sum(abs(V).^2,1));
        m = max(twoNorm);
        V = V/m;

        sqDistances_HCP{subj} = sqDistance(V);

        D = sqDistances_HCP{subj};
        e0 = zeros(t,1);
        d0 = sqDistance_Y(V,e0);

        while min(d0)<=max(max(D(1:n,1:n)))
            e0(1)=e0(1)+0.1;
            d0 = sqDistance_Y(V,e0);
        end

        d0_HCP{subj} = d0;
    end

    for numberofruns = 1:100
        numberofruns        
        l_full_shuffled = randperm(l_full);
        l_full_shuffled_100{numberofruns} = l_full_shuffled;    
        S_opt_all = [];

        for half_subj = 1:2
            K=K_max;
            S_opt_all = [];
            if half_subj == 1
                temp_sqDsitance = sqDistances_HCP(l_full_shuffled(1:l_half));
                temp_d0 = d0_HCP(l_full_shuffled(1:l_half));
                S_opt = exemplar(temp_sqDsitance,temp_d0,t,n,K,l_half);            
                temp = S_opt;
                for K=K_max:-1:2
                    S_opt_all{K} = temp{K_max};
                    temp{K_max}=temp{K_max}(1:end-1);
                end
                index_global = [];             
                for K=2:K_max
                    l=length(sqDistances_HCP);    
                    for subj = l_full_shuffled
                        D = sqDistances_HCP{subj};       
                        ind = S_opt_all{K};
                        [D_sorted,index_global{K}(subj,:)] = min(D(ind,1:end),[],1);
                    end
                end
                clear Maj_half_1 F_half_1 Maj_half_21 F_half_21
                for K = 2:K_max
                    [Maj_half_1{K},F_half_1{K}] = mode(index_global{K}(l_full_shuffled(1:l_half),:));
                    [Maj_half_21{K},F_half_21{K}] = mode(index_global{K}(l_full_shuffled(l_half+1:l_full),:)); 
                    index_global_half_1{K} = index_global{K}(l_full_shuffled(1:l_half));
                    index_global_half_12{K} = index_global{K}(l_full_shuffled(l_half+1:l_full));
                end
                S_opt_all_half_1_allruns{numberofruns} = S_opt_all;
                S_opt_all_half_12_allruns{numberofruns} = S_opt_all;            
                Maj_half_1_allruns{numberofruns} = Maj_half_1;
                Maj_half_21_allruns{numberofruns} = Maj_half_21;
                index_global_half_1_allruns{numberofruns} = index_global_half_1;
                index_global_half_12_allruns{numberofruns} = index_global_half_12;



            elseif half_subj == 2            
                S_opt_all = [];            
                temp_sqDsitance = sqDistances_HCP(l_full_shuffled(l_half+1:l_full));
                temp_d0 = d0_HCP(l_full_shuffled(l_half+1:l_full));
                S_opt = exemplar(temp_sqDsitance,temp_d0,t,n,K,l_half);            
                temp = S_opt;
                for K=K_max:-1:2
                    S_opt_all{K} = temp{K_max};
                    temp{K_max}=temp{K_max}(1:end-1);
                end
                index_global = [];             
                for K=2:K_max
                    l=length(sqDistances_HCP);    
                    for subj = l_full_shuffled
                        D = sqDistances_HCP{subj};       
                        ind = S_opt_all{K};
                        [D_sorted,index_global{K}(subj,:)] = min(D(ind,1:end),[],1);
                    end
                end
                clear Maj_half_22 F_half_22
                for K = 2:K_max
                    [Maj_half_22{K},F_half_22{K}] = mode(index_global{K}(l_full_shuffled(l_half+1:l_full),:)); 
                end            
                S_opt_all_half_22_allruns{numberofruns} = S_opt_all;
                index_global_half_22_allruns{numberofruns} = index_global;
                Maj_half_22_allruns{numberofruns} = Maj_half_22;            
            end
        end
    end

    % Dice Coefficient to find the overlappings (Rest8+Rest9+LR+RL concatenated all)

    clc
    clear reproducibility n1 n2 d1 d2 ratio1 ratio2 dicecof1 dicecof2 index_Maj_1_22 index_Maj_22_1 onehotMaj_half_1 onehot_Maj_half_2_2
    goods = zeros(K_max,1);

    for numberofruns=1:100
        numberofruns
        for K = 2:K_max
            onehot_Maj_half_1{K} = oneHot(Maj_half_1_allruns{numberofruns}{K});
            onehot_Maj_half_22{K} = oneHot(Maj_half_22_allruns{numberofruns}{K});
            A = hist(Maj_half_1_allruns{numberofruns}{K},1:K);
            B = hist(Maj_half_22_allruns{numberofruns}{K},1:K);        
            for k=1:K
                n1{K}(k,:) = 2*sum(repmat(onehot_Maj_half_1{K}(:,k),1,K).*onehot_Maj_half_22{K},1);
                d1{K}(k,:) = A(k)+B;
                ratio1{K}(k,:) = n1{K}(k,:)./d1{K}(k,:);
                [dicecof1{K}(k,1), index_Maj_1_22{K}(k)] = max(ratio1{K}(k,:));
                n2{K}(k,:) = 2*sum(repmat(onehot_Maj_half_22{K}(:,k),1,K).*onehot_Maj_half_1{K},1);
                d2{K}(k,:) = B(k)+A;
                ratio2{K}(k,:) = n2{K}(k,:)./d2{K}(k,:);
                [dicecof2{K}(k,1), index_Maj_22_1{K}(k)] = max(ratio2{K}(k,:));
            end
            reproducibility(K,numberofruns) = mean((dicecof1{K} + dicecof2{K})/2);
            if max(hist(index_Maj_1_22{K},1:K))==1 || max(hist(index_Maj_22_1{K},1:K))==1
                goods(K) = 1;
            end
        end

    end
    reproducibility(1,:)=[]

    close all
    figure
    x = 2:K_max
    y_mean = mean(reproducibility,2)
    err = std(reproducibility,[],2)
    h = errorbar(x,y_mean,err,  '-s','MarkerSize',5,...
        'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',2, 'markersize',15); hold on
    xt = get(gca, 'XTick')
    set(gca, 'FontSize', 16)
    set(gca, 'xlim', [0,27])
    grid on
    xlabel ('Number of networks','FontSize',30,'FontWeight','bold');
    ylabel('Dice coefficient (400^{(11)} - 400^{(22)})','FontSize',30,'FontWeight','bold');
    title(['Dice Coefficient versus Number of Networks'],'FontSize',30,'FontWeight','bold');



    clc

    figure;hold on;
    y=reproducibility;
    mean_y = mean(y,2);
    std_y = std(y,[],2);

    H = shadedErrorBar(x, y', {@mean, @(x) 1*std(x) },{'-s', 'LineWidth', 2}, 0);


    legend([H.mainLine, H.patch], '\sigma', '\mu');
    xt = get(gca, 'XTick');

    set(gca, 'FontSize', 16)
    set(gca, 'xlim', [0,27])
    set(gca, 'ylim', [0.5,0.75])
    grid on
    xlabel ('Number of networks','FontSize',20,'FontWeight','bold');
    ylabel('Dice coefficient (400^1-400^{22})','FontSize',20,'FontWeight','bold');
    title(['Dice coefficient versus number of networks'],'FontSize',20,'FontWeight','bold');
end