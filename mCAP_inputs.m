clear all
close all
%%%
serverdef='media/miplab-nas2/';

%% 1. Path definition
SPM_path=fullfile(filesep,serverdef,'Code','Younes','spm12');
TbCAPs_path=fullfile('./functions/tbCAPs_functions');
mCAPs_path=fullfile('./functions');

%% 2. Data addressing
data_path=fullfile('./data/data_test.mat'); % data should be saved as describe below 
%A. seedmask = A logical vector of number of voxels of whole brain included
% in the study where 1 indicates the voxel within the initial seed -
% bilateral thalami in case of test data. (Dimension : nVox x 1)
%B. TC = Compatible with the work of Bolton.T.W et al 2019. A cell sturcture of timecourses (Dimension: 1 x nParticipants)
%which each cell contains a matri of preprocessed BOLD time-course for
%whole brain (Dimension : lenght_of_timecourse x nVox)
%C. FD = Compatible with the work of Bolton.T.W et al 2019.A matrix of
%frame-wise displacment values (Dimension : lenght_of_timecourse x nParticipants)
out_dir=fullfile('./result_test');

%% 3. Specifying running options
%%% Saving all results
save_all=1;
%%% Doing parallel processing
% Paralell processing option for running the Kmeans iterations faster
parallel_processing=0; % use parfor for Kmeans
numb_pool=7; % number of cores to use in parfor
%% 3. Specifying the CAPs and Kmeans parameters
% Selected K
ik=8;
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

