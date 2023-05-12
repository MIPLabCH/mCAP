%% Finds the moments of (de)activation in a group of fMRI subjects
%% Based on the overlapping weighted timecourses and union . we use mCAP_w_CAP_get_seed_traces function here.!
% Inputs
% tcvox: cell aray with a seed signal in each cell (time points x masked voxels)
% seed: masks for the seeds used (masked voxels x n_seed)
% T: threshold for the retention of active or inactive frames
% FDall: framewise displacement traces for all subjects (time points x n_subjects)
% Mot_thresh: threshold (in mm) to use for scrubbing
%
% Outputs
% Xonp: cell array (each cell dimension masked voxels x n_retained_frames)
% for active frames brain patterns
% p: 5xn_subject matrix with percentages of scrubbed frames, scrubbed and
% active frames, scrubbed and inactive frames, retained frames for active
% analysis and retained frames for inactive analysis
function [Xonp,p,Indices,idx_sep_seeds,Xonp_scrub,S] = mCAP_w_find_activity(tcvox,...
    seed,T,FDall,Mot_thresh,SelMode,SeedType,SignMat,seedmask)
    % Contains the indices (logical) of the selected frames for each seed
    % subcase
    idx_sep_seeds = nan(size(FDall,1),size(FDall,2),size(seed,2));

    % Computes the indices of the peak values that must be excluded
    flag = FDall>Mot_thresh;
    
    % We want to store the indices of scrubbed frames and return it later
    Indices.scrubbed = logical(flag);
    
    
    % Same for the retained frames
    switch SeedType
        case 'Union'
            Indices.kept.active = logical(zeros(size(FDall,1),size(FDall,2)));
            Indices.scrubbedandactive = logical(zeros(size(FDall,1),size(FDall,2)));
        otherwise
            Indices.kept.active = logical(ones(size(FDall,1),size(FDall,2)));
            Indices.scrubbedandactive = logical(ones(size(FDall,1),size(FDall,2)));
    end
    
    % Each cell of flag will contain a time x 1 vector of logicals (1 if
    % the frame is to be censored, 0 otherwise)
    flag = num2cell(flag,1);
    
    % 1 x n_subject vector containing the percentage of scrubbed frames (throughout the whole scan) 
    p_scrubbed = cell2mat(cellfun(@(x) sum(x)/length(x)*100, flag,'un',0));
    
    % If we wanted an average seed, we can either be in the
    % subject-specific seed case (we know there is only one seed type), or
    % in the case with several different seeds (up to three) to average
    % together
    if strcmp(SeedType,'Average')
        
        errordlg('Should never happen...');
        
        
    % If Union or Intersection was wanted instead, then computations for
    % different seeds need to be carried independently
    elseif strcmp(SeedType,'Union')
        for idx = 1:size(seed,2)
      		weight_cap=seed;
      		weightseed=weight_cap(seedmask==1,:);
       		maxv=max(weightseed')';
      		for iv = 1 : size(weightseed,1)
      	    
          		maskforweight(iv,:)=weightseed(iv,:)==maxv(iv,1);
            end
           
      		weightedmask=weightseed.*maskforweight;



            %%% Farnaz: modifed the function to use weighted mean on the
            %%% seedmap
            S = mCAP_w_CAP_get_seed_traces(tcvox,weightedmask(:,idx),SignMat(idx,:),logical(seedmask));
            if strcmp(SelMode,'Threshold')

                % Computes the indexes at which we have a seed activity of interest
                xindp = cellfun(@(x) x>T, S, 'un', 0);

            elseif strcmp(SelMode,'Percentage')

                errordlg('Should never happen...');
            end

            % flag now contains the traces with high activity AND high motion
            flag_active = cellfun(@(x,y) x & y,xindp, flag,'un',0);

            % Vector (1xn_subj) with the percentage of traces removed because of too high
            % motion and being selected as active
            p_scrubactive = cell2mat(cellfun(@(x) sum(x)/length(x)*100, flag_active,'un',0));
            
            % My indices of active/inactive frames now contain only the non
            % corrupted frames
            xindp = cellfun(@(x,y) x & ~y, xindp,flag_active,'un',0);
            
            % Updates the indices of the frames to keep
            switch SeedType
                case 'Intersection'
                    
                    Indices.kept.active = (Indices.kept.active) & cell2mat(xindp);
                    Indices.scrubbedandactive = (Indices.scrubbedandactive) & cell2mat(flag_active);
                    
                case 'Union'
                    Indices.scrubbedandactive = (Indices.scrubbedandactive) | cell2mat(flag_active);
                    Indices.kept.active = (Indices.kept.active) | cell2mat(xindp);
            end
            
            idx_sep_seeds(:,:,idx) = cell2mat(xindp);
            
        end
    end
    
    % Each cell contains the frames selected as active or as inactive (if
    % inactive, the sign is reversed, i.e. inactivation is a positive
    % signal). Size: masked voxels x n_retained_frames
    Xonp = cellfun(@(x,y) x(y,:)',tcvox,num2cell(Indices.kept.active,1),'un',0);
    Xonp_scrub = cellfun(@(x,y) x(y,:)',tcvox,num2cell(logical(Indices.scrubbedandactive),1),'un',0);
    
    % Percentage of active and inactive frames retained per subject
    p_active = cell2mat(cellfun(@(x) size(x,2)/size(FDall,1)*100, Xonp,'un',0));
    
    % Matrix containing all the interesting probabilities
    p = [p_scrubbed; p_scrubactive; p_active];
end
