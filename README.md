# Script used for the micro Co-activation pattern analysis
1- Data should be saved as describe below 
A. seedmask = A logical vector of number of voxels of whole brain included in the study where 1 indicates the voxel within the initial seed -  bilateral thalami in case of test data -(Dimension : nVox x 1)
B. TC = Compatible with the work of Bolton.T.W et al 2019. A cell sturcture of timecourses (Dimension: 1 x nParticipants) which each cell contains a matri of preprocessed BOLD time-course for whole brain (Dimension : lenght_of_timecourse x nVox)
C. FD = Compatible with the work of Bolton.T.W et al 2019.A matrix of

2- The scripts should be run after modification of mCAP_input to this order :
1.mCAP_input
2.mCAP_initialize
3.mCAP_main

3- Prerequisite
SPM12
