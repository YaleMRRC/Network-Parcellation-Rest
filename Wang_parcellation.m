function [index_global_Yeo_17,Maj_Yeo_17,index_global_Yeo_7,Maj_Yeo_7] = wang_parcellation(M,HCP_subj,label_134_cort,Shen268_mapto_Yeo7network,Shen268_mapto_Yeo17network)

    % M is an output of the function "import_HCP_data.m" containing the
    % rest scan data with the following format:
    % M{1,1} = REST_LR, M{1,2} = REST_RL, M{2,1} = REST2_LR, M{2,2} = REST2_RL
    
    % HCP_subj is an output of the function "import_HCP_data.m" containing the
    % subject names with the following format:
    % HCP_subj{1,1} = REST_LR, HCP_subj{1,2} = REST_RL, HCP_subj{2,1} = REST2_LR, HCP_subj{2,2} = REST2_RL
     
    % label_134_cort = a 268x1 vector including node labels according to
    % their position: 1: Subcortical, 3: Cerebellum, 4: Cortical   
    
    % Shen268_mapto_Yeo7network, Shen268_mapto_Yeo17network = the 268x1 
    % vector including the network labels from Yeo's parcellation with 7
    % and 17 networks respectively
        
    % index_global_Yeo_17, index_global_Yeo_7 = individualized parcellation
    % from Wang's algorithm intialized from Yeo's 17 and 7 networks
    % respectively

    % Maj_Yeo_17, Maj_Yeo_7 = Group-level parcellation from Yeo's 7 and 17
    % networks
    
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
    l_full_shuffled = 1:l_full;

    
    
    %%%%%%%%% 17 networks
    Shen268_mapto_Yeo17network = Shen268_mapto_Yeo17network(cortical,:);
    clear B
    for k=1:17
        B_yeo{k} = find(Shen268_mapto_Yeo17network(:,k)==1);
    end
    for subj=1:l_full
        V = [M{1,1}{LR1(subj)}(:,cortical);M{2,1}{LR2(subj)}(:,cortical);M{1,2}{RL1(subj)}(:,cortical);M{2,2}{RL2(subj)}(:,cortical)];
        t = size(V,1);
        n = size(V,2);
        for k=1:17
            B{subj,k} = find(Shen268_mapto_Yeo17network(:,k)==1);
        end
        for k=1:17
            Ref{1}(subj,k,1:t) =  mean(V(:,B{subj,k}),2); % iteration
        end
    end

    h=zeros(0,0);   
    labels = zeros(825,length(cortical));
    for i = 1:10
        i
        for subj=1:l_full
            V = [M{1,1}{LR1(subj)}(:,cortical);M{2,1}{LR2(subj)}(:,cortical);M{1,2}{RL1(subj)}(:,cortical);M{2,2}{RL2(subj)}(:,cortical)];
            t = size(V,1);
            if t<4800
                V(4800,:) = zeros(1,length(cortical));
            end
            corr_with_ref(subj,:,:) = corr(V,squeeze(Ref{i}(subj,:,:))');
            [temp1 , temp2] = sort(squeeze(corr_with_ref(subj,:,:)),2,'descend');
            if temp1(:,2)==0
                confidence(subj,:) = 4;
            else
                confidence(subj,:) = temp1(:,1)./temp1(:,2);
            end
            h(subj) = pdist2(temp2(:,1)',labels(subj,:),'hamming')*length(cortical);
            labels(subj,:) = temp2(:,1);

            for k = 1:17
                B{subj,k} = find(labels(subj,:) == k);
            end
            for k=1:17
                B_confident{k} = B{subj,k}(find(confidence(subj,B{subj,k})>1.1));
                Core{i}(subj,k,:) =  mean(V(:,B_confident{k}),2); % iteration
            end
            Ref{i+1}(subj,:,:) = (2*Core{i}(subj,:,:) + Ref{i}(subj,:,:))/3;
        end
    end

    %
    index_global_Yeo_17 = zeros(825,length(cortical));
    for k = 1:17
        for subj = l_full_shuffled
            index_global_Yeo_17(subj,B{subj,k}) = k;
        end
    end

    Maj_Yeo_17 = zeros(1,length(cortical));
    for k = 1:17
        Maj_Yeo_17(1,B_yeo{k}) = k;
    end


    
    %%%%%%%% 7 Networks
    Shen268_mapto_Yeo7network = Shen268_mapto_Yeo7network(cortical,:);
    clear B

    for k=1:7
        B_yeo{k} = find(Shen268_mapto_Yeo7network(:,k)==1);
    end

    for subj=1:l_full    
        V = [M{1,1}{LR1(subj)}(:,cortical);M{2,1}{LR2(subj)}(:,cortical);M{1,2}{RL1(subj)}(:,cortical);M{2,2}{RL2(subj)}(:,cortical)];

        t = size(V,1);
        n = size(V,2);
        for k=1:7
            B{subj,k} = find(Shen268_mapto_Yeo7network(:,k)==1);
        end
        for k=1:7
            Ref{1}(subj,k,1:t) =  mean(V(:,B{subj,k}),2); % iteration
        end
    end

    h=zeros(0,0);
    labels = zeros(825,length(cortical));

    for i = 1:10
        i
        for subj=1:l_full
            V = [M{1,1}{LR1(subj)}(:,cortical);M{2,1}{LR2(subj)}(:,cortical);M{1,2}{RL1(subj)}(:,cortical);M{2,2}{RL2(subj)}(:,cortical)];
            t = size(V,1);
            if t<4800
                V(4800,:) = zeros(1,length(cortical));
            end
            corr_with_ref(subj,:,:) = corr(V,squeeze(Ref{i}(subj,:,:))');
            [temp1 , temp2] = 	sort(squeeze(corr_with_ref(subj,:,:)),2,'descend');
            if temp1(:,2)==0
                confidence(subj,:) = 4;
            else
                confidence(subj,:) = temp1(:,1)./temp1(:,2);
            end
            h(subj) = pdist2(temp2(:,1)',labels(subj,:),'hamming')*length(cortical);
            labels(subj,:) = temp2(:,1);

            for k = 1:7
                B{subj,k} = find(labels(subj,:) == k);
            end
            for k=1:7
                B_confident{k} = B{subj,k}(find(confidence(subj,B{subj,k})>1.1));
                Core{i}(subj,k,:) =  mean(V(:,B_confident{k}),2); % iteration
            end
            Ref{i+1}(subj,:,:) = (2*Core{i}(subj,:,:) + Ref{i}(subj,:,:))/3;
        end
    end

    index_global_Yeo_7 = zeros(825,length(cortical));
    for k = 1:7
        for subj = 1:l_full
            index_global_Yeo_7(subj,B{subj,k}) = k;
        end
    end

    Maj_Yeo_7 = zeros(1,length(cortical));
    for k = 1:7
        Maj_Yeo_7(1,B_yeo{k}) = k;
    end    
end