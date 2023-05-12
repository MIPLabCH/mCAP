%% Creating a log file in output directory
logFile=fullfile(out_dir,'pipleline_log.txt');
fid=fopen(logFile,'w');

%% Setting initial variables
isconverged=0 ;
it=1;

%% running iterative mCAPs
while isconverged==0 && it <= max_iteration
    %%% selecting frames based on whole seed
    fprintf(fid,['Running iteration number ',num2str(it),'...\n']);
    
    SignMatrix=repmat(SM,size(seed_vol{it,1},2),1);
    Tr=T;
    if it == 1;
        Tr=0.5; %%% initial treshold is lowered to allow maximal frame selection
        [Xon,~,ind_all{it},iss_all{it},~]= CAP_find_activity(TC,seed_vol{it,1},Tr,FD,Tmot,...
            SelMode,SeedType,SignMatrix);
        
    else
        
        [Xon,~,ind_all{it},iss_all{it},~]= mCAP_w_find_activity(TC,seed_vol{it,1},Tr,FD,Tmot,...
            SelMode,SeedType,SignMatrix,seedmask);
    end
    %%% Temporal clustering
    fprintf(fid,['___Kmeans...\n']);
    
    [idx{it},Cm{it},SUMD{it},D{it}] = kmeans(cell2mat(Xon)',ik,'Options',options,'distance','correlation','replicates'...
        ,n_rep,'empty','drop','maxiter',n_iter);
    frames=cell2mat(Xon);
    
    
    %%% Normalizing cap voxels on their temporal variation
    for ii = 1 : ik
        mean_c(:,ii)=mean(frames(:,idx{it}==ii),2);
        std_c(:,ii)=std(frames(:,idx{it}==ii),[],2);
    end
    
    CAP_norm{it}=mean_c./std_c; clear mean_c std_c
    
    %%% creating seed mask based on the normalized cap
    
    seed_map=CAP_norm{it}(seedmask,:);

    
    newseed=repmat(seedmask,1,ik); %% the new seed has ik parcels
    newseed(seedmask==1,:)=seed_map;
    
    
    
    it=it+1;
    seed_vol{it,1}=newseed; clear newseed ;
    
    
    %%% assesing convergence comparimg new seed with previous seed (only done for second seed onwards)
    if it >2
        [dis_final] = mCAP_CAP_CompareSeeds(seed_vol{it,1}(seedmask==1,:),seed_vol{it-1,1}(seedmask==1,:));
        mean_dis(it,1)=mean(dis_final);
        if mean_dis(it) <= convergence_seed_value
            isconverged=1;
        end
    end
  
    
end

fprintf(fid,['******** WHILE LOOP Ended...\n']);

%% Gathering variables to save
mean_dis(1:2,:)=[]; %% we didnt asses in the first 2 iterations so these are by default zero

%%% gathering variables to save
mat_save.CAP_norm=CAP_norm; clear CAP_norm;
mat_save.Cm=Cm;clear Cm;
mat_save.D=D; clear D;
mat_save.idx=idx; clear idx;
mat_save.ind_all=ind_all; clear ind_all;
mat_save.iss_all=iss_all; clear iss_all;
mat_save.mean_dis=mean_dis; clear mean_dis;
mat_save.seed_vol=seed_vol; clear seed_vol;
mat_save.SUMD=SUMD; clear SUMD;
mat_save.lastXON=Xon; clear Xon frames;


if save_all == 1;
fprintf(fid,['__________________ saving variables...\n']);

%%% saving resuts fo ik
save_dir=fullfile(out_dir,['mCAPs_conv_',char(num2str(isconverged)),'_K_',char(num2str(ik))]);
mkdir(save_dir);
save(fullfile(save_dir,['mCAP_K',char(num2str(ik)),'.mat']),'mat_save','it','-v7.3'); clear mat_save it;
fprintf(fid,['Done without errors...\n']);
end
