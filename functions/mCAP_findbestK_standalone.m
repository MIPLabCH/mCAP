clear all
close all
%%%
serverdef='/media/miplab-nas2/';

%% 1. Path definition
SPM_path=fullfile(filesep,serverdef,'Code','Younes','spm12');
TbCAPs_path=fullfile('./functions/tbCAPs_functions');
mCAPs_path=fullfile('./functions');

%% 2. Data addressing
data_path=fullfile('./data/data_test.mat');
out_dir=fullfile('/TT_result_test');

%% 3. Specifying running options
%%% Saving all results
save_all=1;
%%% Doing parallel processing
% Paralell processing option for running the Kmeans iterations faster
parallel_processing=1; % use parfor for Kmeans
numb_pool=7; % number of cores to use in parfor
%% 3. Specifying the CAPs and Kmeans parameters
% Selected ranage of K
K_all=(2:10);
% Max iteration number
max_iteration=6;
% The distance define as convergence across seeds
convergence_seed_value=0.005;
% Threshold above which to select frames
T =1;
% Threshold of FD above which to scrub out the frame and also the t-1 and
% t+1 frames 
Tmot = 0.5;
% Number reps and iterations for Kmeans
N = 40;
n_rep=100;
n_iter=300;


%% 4. Adding paths
addpath(SPM_path);
addpath(TbCAPs_path);
addpath(mCAPs_path)
%%% Making output directory
if ~isdir(out_dir);  mkdir(out_dir); end


%% 5. Loading data

load(data_path);
data_path=fullfile('./data/data_test.mat'); % data should be saved as describe below 
%A. seedmask = A logical vector of number of voxels of whole brain included
% in the study where 1 indicates the voxel within the initial seed -
% bilateral thalami in case of test data. (Dimension : nVox x 1)
%B. TC = Compatible with the work of Bolton.T.W et al 2019. A cell sturcture of timecourses (Dimension: 1 x nParticipants)
%which each cell contains a matri of preprocessed BOLD time-course for
%whole brain (Dimension : lenght_of_timecourse x nVox)
%C. FD = Compatible with the work of Bolton.T.W et al 2019.A matrix of
%frame-wise displacment values (Dimension : lenght_of_timecourse x nParticipants)

%%% deviding data into test and train frames


nsub=size(TC,2);
%randomize
order_rand=randperm(nsub);
for i =1: nsub
TC_new{1,i}=TC{1,order_rand(i)};
end
TC=TC_new; clear TC_new;
% deviding in half 
n1=floor(nsub/2);
n2=nsub-n1;
TC_all=TC; clear TC; 
FD_all=FD(:,order_rand); clear FD; 
for iin= 1 : n1 
    TCa{iin} = TC_all{iin};
    FDa(:,iin)=FD_all(:,iin);
end
for iin = n1+1 : nsub 
    TCb{iin-n1} = TC_all{iin};
    FDb(:,iin-n1)=FD_all(:,iin);
end
clear FD_all; clear TC_all;

%% 6.  Unchange parameters of mCAP
% Selection mode ('Threshold' )
SelMode = 'Threshold';
% Type of used seed information: 'Union' is flagged as active if any of the
% seeds are active
SeedType = 'Union';
% retain frames if seed is active
SM = [1,0]; %% activation

% Percentage of positive-valued voxels to retain for clustering
Pp = 100; %%check this value!!
% Percentage of negative-valued voxels to retain for clustering
Pn =100;
% Percentage of frames to use in each fold of consensus clustering
Pcc = 80;
%% 7. Initialize parpool if indicated 

if parallel_processing==1
    
    delete(gcp('nocreate'));
    parpool(numb_pool);
    
    
    stream = RandStream('mlfg6331_64');  % Random number stream
    options = statset('UseParallel',1,'UseSubstreams',1,...
        'Streams',stream);
    
else
    options=[];
end


%% 8. Main
for kinvestigate = 1 : length(K_all)
seed_vol{1,1}=logical(seedmask);
seedmask=logical(seedmask);   
ik=K_all(kinvestigate);

%%% creating a log file for this K

logFile=fullfile(out_dir,['mCAP_TT_k',char(num2str(ik)),'_log.txt']);
fid=fopen(logFile,'w');




clear mean_dis
isconverged=0 ;
it=1;
while isconverged==0 && it <max_iteration
    
%%% selecting test frames based on whole seed
fprintf(fid,['Running iteration number ',num2str(it),'...\n']);

SignMatrix=repmat(SM,size(seed_vol{it,1},2),1);
Tr=T;
if it == 1 
    Tr= 0.5; 
    [Xonb,~,~,~]= CAP_find_activity(TCb,seed_vol{it,1},Tr,FDb,Tmot,...
   SelMode,SeedType,SignMatrix);
else
[Xonb,~,~,~]= mCAP_w_find_activity(TCb,seed_vol{it,1},Tr,FDb,Tmot,...
   SelMode,SeedType,SignMatrix,seedmask);
end

testframes=cell2mat(Xonb); clear Xonb;
    
    
%%% selecting train frames based on whole seed
fprintf(fid,['_____ CAPs...\n']);

SignMatrix=repmat(SM,size(seed_vol{it,1},2),1);
Tr=T;
if it == 1 
    Tr= 0.5; 
    [Xon,~,ind_all{it,1},iss_all{it,1}]= CAP_find_activity(TCa,seed_vol{it,1},Tr,FDa,Tmot,...
   SelMode,SeedType,SignMatrix);

else

[Xon,~,ind_all{it,1},iss_all{it,1}]= mCAP_w_find_activity(TCa,seed_vol{it,1},Tr,FDa,Tmot,...
   SelMode,SeedType,SignMatrix,seedmask);
end


%%% Temporal clustering 
fprintf(fid,['___Kmeans...\n']);

[idx{it},Cm{it},SUMD{it},D{it}] = kmeans(cell2mat(Xon)',ik,'Options',options,'distance','correlation','replicates'...
   ,n_rep,'empty','drop','maxiter',n_iter);
frames=cell2mat(Xon);

%%% obtaining caps improvement 

[emptk(it,1),maxSUMD(it,1),frame_assign{it},~]=mCAP_tframe_m(testframes,Cm{it});

%%% normalizing caps
 for ii = 1 : ik
     mean_c(:,ii)=mean(frames(:,idx{it}==ii),2);
     std_c(:,ii)=std(frames(:,idx{it}==ii),[],2);
 end

CAP_norm{it}=mean_c./std_c; clear mean_c std_c
 
%%% creating input for weighted timecourse

newseed=CAP_norm{it}.*seedmask;
it=it+1;
seed_vol{it,1}=newseed; clear newseed;

%%% assesing convergence 
if it >2
    mnew=seed_vol{it,1}; mold=seed_vol{it-1,1};
    [dis_final] = mCAP_CAP_CompareSeeds(mnew(seedmask==1,:),mold(seedmask==1,:));
    mean_dis(it,1)=mean(dis_final);
 
%%%comment if you want to run through a fixed number of iterations
     if mean_dis(it) <=convergence_seed_value
         isconverged=1;
     end
end


end


fprintf(fid,['WHILE LOOP FINISHED_ Last K-means...\n']);

%%% gathering variables to save
mat_save.frame_assign=frame_assign; clear frame_assign;
mat_save.CAP_norm=CAP_norm; clear CAP_norm;
mat_save.Cm=Cm;clear Cm;
mat_save.D=D; clear D;
mat_save.idx=idx; clear idx;
mat_save.ind_all=ind_all; clear ind_all;
mat_save.iss_all=iss_all; clear iss_all;
mat_save.mean_dis=mean_dis; %clear mean_dis;
%mat_save.mean_dis_cap_first=mean_dis_cap_first;
mat_save.seed_vol=seed_vol; clear seed_vol;
mat_save.SUMD=SUMD; clear SUMD;
mat_save.lastXON=Xon; clear Xon frames;
mat_save.last_testframes=testframes; clear testframes;
mat_save.emptk=emptk; 
mat_save.maxSUMD=maxSUMD; mean_max_sumD=mean(maxSUMD); clear maxSUMD;

fprintf(fid,['__________________ saving variables...\n']);


%%% saving resuts fo ik 
save_dir=fullfile(out_dir,['TTmCAP_conv',char(num2str(isconverged)),'_K',char(string(ik))]);
mkdir(save_dir);
save(fullfile(save_dir,['mCAPTT_K',char(num2str(ik)),'.mat']),'emptk','mean_max_sumD','mat_save','mean_dis','it','-v7.3'); clear mat_save it;
fprintf(fid,['Done without errors...\n']);


end




% Copyright 2023 Farnaz Delavari
%
% Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

