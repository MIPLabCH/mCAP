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
