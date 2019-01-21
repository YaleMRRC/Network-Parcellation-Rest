function [disimilirities_global,reproducibility] = rest1_rest2_reproducibility(M,HCP_subj,S_opt_all,label_134_cort,K_max,local_global_flag)

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

    % local_global_flag = the flag indicating whether to compare REST1 and REST2
    % with the same exempalrs across the two runs (global, '1') or using different
    % exemplars across the two runs (local, '0') 
    % For global exemplars, we laready have done REST1 and REST2 using
    % concatenation.
    
    % similirities_global = parcellation scheme similarity between REST1
    % and REST2 calculated by hammign distance using global exemplars
    
    % reproducibility = dice coefficient between parcellation schemes of
    % REST1 and REST2 using global or local exemplars


    cortical_boolean = input('Do you want to consider only the cortical regions? (1: cortical, 0: whole brain)')
    if cortical_boolean == 1
        vect = [3];
    elseif cortical_boolean == 0
        vect = [1,3,4];
    end
    cortical = label_134_cort(ismember(label_134_cort(:,2),vect),1);

    
    subjsLR = intersect(HCP_subj{1,1},HCP_subj{2,1});
    subjsRL = intersect(HCP_subj{1,2},HCP_subj{2,2});
    subjs_all = intersect(subjsLR,subjsRL);

    LR1 = arrayfun(@(x)find(HCP_subj{1,1}==x,1),subjs_all);
    LR2 = arrayfun(@(x)find(HCP_subj{2,1}==x,1),subjs_all);
    RL1 = arrayfun(@(x)find(HCP_subj{1,2}==x,1),subjs_all);
    RL2 = arrayfun(@(x)find(HCP_subj{2,2}==x,1),subjs_all);


    HCP_subj_common{1,1} = LR1;
    HCP_subj_common{1,2} = RL1;
    HCP_subj_common{2,1} = LR2;
    HCP_subj_common{2,2} = RL2;

    if local_global_flag == 0 
    % Calculating the REST1 and REST2 separately - Local exemplars
        for rest = 1:2 
            rest
            LR = HCP_subj_common{rest,1};
            RL = HCP_subj_common{rest,2};
            l_full = length(LR);

            for subj=1:l_full        
                V = [M{rest,1}{LR(subj)}(:,cortical);M{rest,2}{RL(subj)}(:,cortical)];
                t = size(V,1);
                n = size(V,2);
                mean_subtract = mean(V,2);
                V = V - repmat(mean_subtract,[1,n]);   

                twoNorm = sqrt(sum(abs(V).^2,1));
                m = max(twoNorm);
                V = V/m;

                sqDistances_HCP{rest,subj} = sqDistance(V);

                D = sqDistances_HCP{rest,subj};
                e0 = zeros(t,1);
                e0(1)=e0(1)+3;
                d0 = sqDistance_Y(V,e0);

                if min(d0)<=max(max(D(1:n,1:n)))
                    fprintf('No :( \n')
                end

                d0_HCP{rest, subj} = d0;
            end
        end


        clear Maj F index_global 
        for rest = 1:2
            rest
            temp_sqDsitance = sqDistances_HCP(rest,:);
            temp_d0 = d0_HCP(rest,:);
            S_opt{rest} = exemplar(temp_sqDsitance,temp_d0,t,n,K_max,l_full);


            temp = S_opt{rest};
            for K=K_max:-1:2
                S_opt_all{K,rest} = temp;
                temp=temp(1:end-1);
            end

            for K=2:K_max  
                for subj = 1:l_full
                    D = sqDistances_HCP{rest,subj};
                    ind = S_opt_all{K,rest};
                    [D_sorted,index_global{K,rest}(subj,:)] = min(D(ind,1:end),[],1);
                end
                [Maj{K,rest},F{K,rest}] = mode(index_global{K,rest});
            end            
        end
    end    
    
    if local_global_flag == 1
    % Calculating the REST1 and REST2 separately - Global exemplars
        for rest = 1:2 
            rest
            LR = HCP_subj_common{rest,1};
            RL = HCP_subj_common{rest,2};
            l_full = length(LR);

            for subj=1:l_full        
                V = [M{rest,1}{LR(subj)}(:,cortical);M{rest,2}{RL(subj)}(:,cortical)];

                t = size(V,1);
                n = size(V,2);
                mean_subtract = mean(V,2);
                V = V - repmat(mean_subtract,[1,n]);   

                twoNorm = sqrt(sum(abs(V).^2,1));
                m = max(twoNorm);
                V = V/m;

                sqDistances_HCP{rest,subj} = sqDistance(V);

                D = sqDistances_HCP{rest,subj};
                e0 = zeros(t,1);
                e0(1)=e0(1)+3;
                d0 = sqDistance_Y(V,e0);

                if min(d0)<=max(max(D(1:n,1:n)))
                    fprintf('No :( \n')
                end

                d0_HCP{rest, subj} = d0;
            end
        end

%         load /home/mehraveh/documents/MATLAB/Parcellation/Oct3/BrainFigures_8_9_concat_825_cortical.mat
    S_opt_all_all_rest = S_opt_all;

        for rest = 1:2
            for K=2:K_max  
                for subj = 1:l_full
                    D = sqDistances_HCP{rest,subj};
                    ind = S_opt_all_all_rest{K};
                    [D_sorted,index_global_all_rest{K,rest}(subj,:)] = min(D(ind,1:end),[],1);
                end
                [Maj_all_rest{K,rest},F_all_rest{K,rest}] = mode(index_global_all_rest{K,rest});
            end            
        end
    end
    
    if local_global_flag == 1
    % Comparing REST1 and REST2 for each subject with global exemplars

%     load /home/mehraveh/documents/MATLAB/Parcellation/2017/BrainFigures_rest_1_2_cortical_global_exemplars_per_rest.mat
        clear similirities_global
        l_full = length(index_global_all_rest{2,1});

        for K=2:K_max
            K
            for subj=1:l_full     
                disimilirities_global(K,subj) = pdist2(index_global_all_rest{K,1}(subj,:),index_global_all_rest{K,2}(subj,:),'hamming');               
            end
        end
    end

    % Dice Coefficient to compare REST1 and REST2 for both global and local exempalrs

    load /home/mehraveh/documents/MATLAB/Parcellation/2017/BrainFigures_rest_1_2_cortical_local_exemplars_per_rest.mat
    load /home/mehraveh/documents/MATLAB/Parcellation/2017/BrainFigures_rest_1_2_cortical_global_exemplars_per_rest.mat
    l_full = length(index_global_all_rest{2,1})

    clc
    clear reproducibility n1 n2 d1 d2 ratio1 ratio2 dicecof1 dicecof2
    goods = zeros(K_max,1);


    for K = 2:K_max
        K
        for subj=1:l_full
            if local_global_flag == 0
%             local exemplars
                onehot_REST1 = oneHot(index_global{K,1}(subj,:));
                onehot_REST2 = oneHot(index_global{K,2}(subj,:));
                A = hist(index_global{K,1}(subj,:),1:K);
                B = hist(index_global{K,2}(subj,:),1:K);
            elseif local_global_flag == 1
%             global exemplars
                onehot_REST1 = oneHot(index_global_all_rest{K,1}(subj,:));
                onehot_REST2 = oneHot(index_global_all_rest{K,2}(subj,:));
                A = hist(index_global_all_rest{K,1}(subj,:),1:K);
                B = hist(index_global_all_rest{K,2}(subj,:),1:K);
            end
            for k=1:K
                n1{K}(k,:) = 2*sum(repmat(onehot_REST1(:,k),1,K).*onehot_REST2,1);
                d1{K}(k,:) = A(k)+B;
                ratio1{K}(k,:) = n1{K}(k,:)./d1{K}(k,:);
                [dicecof1{K}(k,1), index_REST1_REST2{K}(k)] = max(ratio1{K}(k,:));
                n2{K}(k,:) = 2*sum(repmat(onehot_REST2(:,k),1,K).*onehot_REST1,1);
                d2{K}(k,:) = B(k)+A;
                ratio2{K}(k,:) = n2{K}(k,:)./d2{K}(k,:);
                [dicecof2{K}(k,1), index_REST2_REST1{K}(k)] = max(ratio2{K}(k,:));
            end
            reproducibility(K,subj) = mean((dicecof1{K} + dicecof2{K})/2);        
        end
    end
    reproducibility(1,:)=[]
end
