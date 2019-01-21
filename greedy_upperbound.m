function [GREEDY,UpperBound] = greedy_upperbound(M,HCP_subj,S_opt_all,label_134_cort,K_max)

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
    
    % GREEDY = the greedy estimation of optimum value
    % UpperBound = the upperbound to the optimum value calculated in a
    % data-drive manner
       

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

    K = K_max;   
    GREEDY = zeros(l_full,K);  
    UpperBound = zeros(l_full,K);
    delta = zeros(l_full,K,n);
    for subj = 1:l_full    
        subj
        D = sqDistances_HCP{subj};
        d0 = d0_HCP{subj};
        D(n+1,1:n) = d0;
        D(1:n,n+1) = d0;    
        for k = 1:K  

            S = S_opt_all{k};
            GREEDY(subj,k) = myfunc(D,S,n);      
            background = setdiff(1:188,S);
            for v = background
                delta(subj,k,v) = myfunc(D,[S,v],n) - myfunc(D,S,n);    
            end

            UpperBound(subj,k) =  GREEDY(subj,k) + sum(delta(subj,k));
        end
    end


    GREEDY(:,1) = [];
    UpperBound(:,1) = [];

    GREEDY_SUM = sum(GREEDY);
    GREEDY_SUM = (GREEDY_SUM-min(GREEDY_SUM))/(max(GREEDY_SUM)-min(GREEDY_SUM));
    UpperBound_SUM = sum(UpperBound);
    UpperBound_SUM = (UpperBound_SUM-min(UpperBound_SUM))/(max(UpperBound_SUM)-min(UpperBound_SUM));


    figure
    plot(sum(GREEDY)), hold on
    % plot((1-(1/exp(1)))*GREEDY_SUM), hold on
    plot(sum(UpperBound))

    figure
    plot(GREEDY_SUM), hold on
    % plot((1-(1/exp(1)))*GREEDY_SUM), hold on
    plot(UpperBound_SUM)
end
