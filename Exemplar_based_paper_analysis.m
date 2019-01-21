function [X_train,y_train,X_test,y_test] = cross_validated_parcellation_10fold(M,HCP_subj,S_opt_all,label_134_cort,K_max,gender_825)

    % M is an output of the function "import_HCP_data.m" containing the
    % rest scan data with the following format:
    % M{1,1} = REST_LR, M{1,2} = REST_RL, M{2,1} = REST2_LR, M{2,2} = REST2_RL
    
    % HCP_subj is an output of the function "import_HCP_data.m" containing the
    % subject names with the following format:
    % HCP_subj{1,1} = REST_LR, HCP_subj{1,2} = REST_RL, HCP_subj{2,1} = REST2_LR, HCP_subj{2,2} = REST2_RL
     
    % label_134_cort = a 268x1 vector including node labels according to
    % their position: 1: Subcortical, 3: Cerebellum, 4: Cortical
    % K_max = The number of networks, e.g. 25
    
    % S_opt_all = the global exempalrs calculated using the concatenation
    % of REST1 and REST2. This is an output of the function
    % "Parcellation_rest.m"
    
    % gender_825 = sex labels for all subjects


    subjsLR = intersect(HCP_subj{1,1},HCP_subj{2,1});
    subjsRL = intersect(HCP_subj{1,2},HCP_subj{2,2});
    subjs_all = intersect(subjsLR,subjsRL);

    LR1 = arrayfun(@(x)find(HCP_subj{1,1}==x,1),subjs_all);
    LR2 = arrayfun(@(x)find(HCP_subj{2,1}==x,1),subjs_all);
    RL1 = arrayfun(@(x)find(HCP_subj{1,2}==x,1),subjs_all);
    RL2 = arrayfun(@(x)find(HCP_subj{2,2}==x,1),subjs_all);


    cortical_boolean = input('Do you want to consider only the cortical regions? (1: cortical, 0: whole brain)')
    if cortical_boolean == 1
        vect = [3];
    elseif cortical_boolean == 0
        vect = [1,3,4];
    end
    cortical = label_134_cort(ismember(label_134_cort(:,2),vect),1);


    l_full = length(LR1)
    for subj=1:l_full
        V = [M{1,1}{LR1(subj)}(:,cortical);M{2,1}{LR2(subj)}(:,cortical);M{1,2}{RL1(subj)}(:,cortical);M{2,2}{RL2(subj)}(:,cortical)];

        t = size(V,1);
        n = size(V,2);
        mean_subtract = mean(V,2);
        V = V - repmat(mean_subtract,[1,n]);   
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

    new_order = randperm(l_full);
    sqDistanceNew = sqDistances_HCP(new_order);
    d0New = d0_HCP(new_order);
    genderNew = gender(new_order);

    indices = crossvalind('Kfold', l_full, 10)

    clear indice_train indice_test S_opt_all_train
    for fold = 1:10
        indice_train = (indices~=fold);
        indice_test = (indices==fold);
        l_train = sum(indice_train);
        S_opt = exemplar(sqDistanceNew(indice_train),d0New(indice_train),t,n,25,l_train);           

        temp = S_opt;
        for K=25:-1:2
            S_opt_all_train{fold,K} = temp;
            temp=temp(1:end-1);
        end
    end

    for fold = 1:10
        indice_train = (indices~=fold);
        indice_test = (indices==fold);
        l_train = sum(indice_train);

        for K=2:25  
            for subj = 1:l_full
                D = sqDistanceNew{subj};
                ind = S_opt_all_train{fold,K};
                [D_sorted,index_global_test{fold,K}(subj,:)] = min(D(ind,1:end),[],1);
            end

            X_train{fold,K} = index_global_test{fold,K}(indice_train,:);
            y_train{fold} = genderNew(indice_train,:);
            X_test{fold,K} = index_global_test{fold,K}(indice_test,:);
            y_test{fold} = genderNew(indice_test,:);  
        end
    end
end

