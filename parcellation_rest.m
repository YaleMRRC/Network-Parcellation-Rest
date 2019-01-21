function [S_opt_all, index_global, Maj, F] = parcellation_rest(M,HCP_subj,label_134_cort,K_max)
    

    % M is an output of the function "import_HCP_data.m" containing the
    % rest scan data with the following format:
    % M{1,1} = REST_LR, M{1,2} = REST_RL, M{2,1} = REST2_LR, M{2,2} = REST2_RL
    
    % HCP_subj is an output of the function "import_HCP_data.m" containing the
    % subject names with the following format:
    % HCP_subj{1,1} = REST_LR, HCP_subj{1,2} = REST_RL, HCP_subj{2,1} = REST2_LR, HCP_subj{2,2} = REST2_RL

    % label_134_cort = a 268x1 vector including node labels according to
    % their position: 1: Subcortical, 3: Cerebellum, 4: Cortical
    
    % K_max = The number of networks, e.g. 25
    
    
    % S_opt_all = exempalrs calculated using the entire population data
    % index_global = individualized parcellation scheme for each subject
    % Maj = group-level parcellation scheme, calculated using majority
    % voting.    
    % F = The frequency of the votes for eahc node-to-network assignment in
    % the group-level parcellation
    
    % Parcellation with all subjects - Rest8+Rest9+LR+RL concatenated all)

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

    S_opt = exemplar(sqDistances_HCP,d0_HCP,t,n,K,l_full);

    temp = S_opt;
    for K=K_max:-1:2
        S_opt_all{K} = temp;
        temp=temp(1:end-1);
    end

    for K=2:K_max
        for subj = 1:l_full
            D = sqDistances_HCP{subj};
            ind = S_opt_all{K};
            [D_sorted,index_global{K}(subj,:)] = min(D(ind,1:end),[],1);
        end
        [Maj{K},F{K}] = mode(index_global{K});  
    end            
end


