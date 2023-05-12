
%% 1. Adding paths
addpath(SPM_path);
addpath(TbCAPs_path);
addpath(mCAPs_path)
%%% Making output directory
if ~isdir(out_dir);  mkdir(out_dir); end


%% 2. Loading data

load(data_path);


%% 3.  Unchange parameters of mCAP
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
%% 4. Initialize parpool if indicated 

if parallel_processing==1
    
    delete(gcp('nocreate'));
    parpool(numb_pool);
    
    
    stream = RandStream('mlfg6331_64');  % Random number stream
    options = statset('UseParallel',1,'UseSubstreams',1,...
        'Streams',stream);
    
else
    options=[];
end

%% 5.  Defining the whole seed is first seed. 
seed_vol{1,1}=logical(seedmask);
seedmask=logical(seedmask);

