function [DB,DB_Shen,DB_Wang_7,DB_Wang_17,Dunn,Dunn_Shen,Dunn_Wang_7,Dunn_Wang_17] = Comparison_DB_Dunn(M,HCP_subj,label_134_cort,K_max,index_global,Maj,index_global_Shen,Maj_Shen,index_global_Wang_7,Maj_Wang_7,index_global_Wang_17,Maj_Wang_17)

    % M is an output of the function "import_HCP_data.m" containing the
    % rest scan data with the following format:
    % M{1,1} = REST_LR, M{1,2} = REST_RL, M{2,1} = REST2_LR, M{2,2} = REST2_RL
    
    % HCP_subj is an output of the function "import_HCP_data.m" containing the
    % subject names with the following format:
    % HCP_subj{1,1} = REST_LR, HCP_subj{1,2} = REST_RL, HCP_subj{2,1} = REST2_LR, HCP_subj{2,2} = REST2_RL
     
    % label_134_cort = a 268x1 vector including node labels according to
    % their position: 1: Subcortical, 3: Cerebellum, 4: Cortical
    % K_max = The number of networks, e.g. 25
    
    % index_global = the individualzied parcellation scheme using
    % exemplar-based parcellation approach
    
    % Maj = the group-level parcellation scheme using exemplar-based
    % parcellation approach
    
    % index_global_Shen = the individualzied parcellation scheme from
    % Shen's parcellation algorithm
    
    % Maj_Shen = the group-level parcellation scheme from Shen parcellation
    % algorithm

    % index_global_7, index_global_17  = the individualzied parcellation
    % scheme using Wang's algorithm initializing from Yeo's 7 and 17 
    % networks respectively.
    
    % Maj_Wang_7, Maj_Wang_17 = the group-level parcellation scheme from Yeo's 7 and
    % 17 network

     
   
    % addpath(genpath('/home/mehraveh/documents/MATLAB/somtoolbox/'))
 
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

    l_full = length(LR1);
    l_full = 1:l_full;

    %%%% Calculating Davies-Bouldin index   
    clear DB_Wang_17 DB_Wang_7 DB_Shen DB
    for subj=1:l_full
       V = [M{1,1}{LR1(subj)}(:,cortical);M{2,1}{LR2(subj)}(:,cortical);M{1,2}{RL1(subj)}(:,cortical);M{2,2}{RL2(subj)}(:,cortical)];
       t = size(V,1);
       n = size(V,2);

       % Wang 17
       temp = index_global_Wang_17(subj,:);
       [t,r] = db_index(V', temp');
       DB_Wang_17(subj) = t;

       % Wang 7
       temp = index_global_Wang_7(subj,:);
       [t,r] = db_index(V', temp');
       DB_Wang_7(subj) = t;


       for K = 2:25

           % Exempalr-based
           temp = index_global{K}(subj,:);
           [t,r] = db_index(V', temp');
           DB(subj,K) = t;

           % Shen
           for alpha = 1:7
               temp = squeeze(index_global_Shen(subj,:,alpha,K));
            if sum(temp) == 0
                   continue
               else
                   [t,r] = db_index(V', temp');
                   DB_Shen(subj,K,alpha) = t;
               end
           end
       end

    end

    %%%% Calculating Dunn index

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
       e0(1) = 3;
       d0 = sqDistance_Y(V,e0); 

       d0_HCP{subj} = d0;
    end


    clear Dunn_Wang_17 Dunn_Wang_7 Dunn Dunn_Shen
    for subj=1:l_full
       subj
       D = sqDistances_HCP{subj};

        % Wang 17 
       if isempty(dunns(17,D,index_global_Wang_17(subj,:)))
           subj
           continue
       else
           Dunn_Wang_17(subj) =dunns(17,D,index_global_Wang_17(subj,:));

       end
       % Wang 7
       Dunn_Wang_7(subj) = dunns(7,D,index_global_Wang_7(subj,:));   


       for K = 2:25
           % Exemplar-based
           Dunn(K,subj) = dunns(K,D,index_global{K}(subj,:));
           % Shen
           for alpha = 1:7
               temp = squeeze(index_global_Shen(subj,:,alpha,K));
            if sum(temp) == 0
                   continue
               else
                   Dunn_Shen(K,subj,alpha) = dunns(K,D,index_global_Shen(subj,:,alpha,K));
               end
           end
       end
    end   
end
