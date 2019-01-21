function [M, HCP_subj] = import_HCP_data(datadir)
    
    % datadir = the path to the HCP data set scuh that REST1 
    % and REST2 with Left-Right and Right-Left acquisitions are stored in
    % separate folders: REST_LR, REST_RL, REST2_LR, REST2_RL.
    % e.g. datadir = '/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/'    
    
    % addpath(genpath('/home/mehraveh/documents/MATLAB/Parcellation/'))
    % filesRest1LR = dir(['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST_LR/matrices/*_GSR_roimean.txt']);
    % filesRest1RL = dir(['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST_RL/matrices/*_GSR_roimean.txt']);
    % filesRest2LR = dir(['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST2_LR/matrices/*_GSR_roimean.txt']);
    % filesRest2RL = dir(['/mnt/store1/mridata2/mri_group/HCP_data/HCP_900_DATA/REST2_RL/matrices/*_GSR_roimean.txt']);

    addpath(genpath(datadir))
    filesRest1LR = dir([datadir,'REST_LR/matrices/*_GSR_roimean.txt']);   
    filesRest1RL = dir([datadir,'REST_RL/matrices/*_GSR_roimean.txt']);
    filesRest2LR = dir([datadir,'REST2_LR/matrices/*_GSR_roimean.txt']);
    filesRest2RL = dir([datadir,'REST2_RL/matrices/*_GSR_roimean.txt']);

    files = {'filesRest1LR','filesRest1RL';'filesRest2LR','filesRest2RL'}


    for task = 1:2
        for lr = 1:2
        file_cur = eval(files{task,lr});    
            for subj = 1:length(file_cur)
                temp = importdata(file_cur(subj).name);
                M{task,lr}{subj} = temp.data(:,2:end);
            end
        end
    end
    
    save ([datadir,'M_Rests.mat','M');

    for task = 1:2
        for lr=1:2
            subj_data=zeros(0,1);
            file_cur = eval(files{task,lr});
            for i = 1:length(file_cur)
                if (file_cur(i).bytes) ~= 0
                    out=regexp(file_cur(i).name,'\d+','match');
                    subj_data=[subj_data;str2double(cat(1,out{1}))];
                end
            end
            HCP_subj{task,lr}=subj_data;
        end
    end
    save ([datadir,'HCP_subj.mat','HCP_subj');
end
