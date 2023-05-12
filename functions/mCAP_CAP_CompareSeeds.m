function [dis_final] = mCAP_CAP_CompareSeeds(seed1,seed2)
    seed1=double(seed1);
    seed2=double(seed2);
    % Cosine distance between seeds
    DIS = pdist2(seed1',seed2','cosine');

    n_voxels = sum(sum(seed1));
    
    % Matching pairs
    IDX = munkres(DIS);
    
    idxk = 1;

    for k=1:size(DIS,1)
        if sum(IDX(k,:)) == 1
            dis_final(idxk) = DIS(k,IDX(k,:));
            idxk = idxk + 1;
        else

        end
    end

    dis_final = dis_final;%/n_voxels;
    
    
    %% check all the points 
end
