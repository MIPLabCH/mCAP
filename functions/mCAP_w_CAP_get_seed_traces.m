%% From the raw voxelwise data, computes a spatially WEIGHTED averaged seed signal
%% 1.09.2021 modified t include sum of abs values in weight demin
%% 4.10.2021 modified to include squared sum ob values in weight demin 
% tcvox is a cell array (each cell has time x masked voxels dimension)
% seedmask is the seed mask (masked voxel x 1)
% seeda is the signal within the seed normalized based on temporal std
% which we use as weight for the weighted mean here.
% S is a cell array (each cell is a time x 1 vector)
function [S] = mCAP_w_CAP_get_seed_traces(tcvox,seeda,SignMat,seedmask)
%seedm=zscore(seeda(seedmask),[],1);

powerofsum1=1;
powerofsum2=1;
weightvalue=seeda;

weightpower=weightvalue.^powerofsum1;
denim=sum((abs(weightvalue)).^powerofsum2);
%denim=sum((weightvalue).^powerofsum2);

%%% remove the negative values and replace them with zero
%seedm(seedm<0)=0;
%seedm=seeda(seedmask)-mean(seeda(seedmask));
    if SignMat(1)

        S = cellfun(@(x) zscore((sum((x(:,seedmask).*weightpower'),2))./denim), tcvox, 'un', 0);
        % Computation of the spatially averaged seed signal
        %S = cellfun(@(x) zscore(sum((x(:,seedmask).*abs(seedm)'),2)./sum((abs(seedm)))), tcvox, 'un', 0);
        %S2 = cellfun(@(x) sum(x(:,seed),2), tcvox, 'un', 0);
        %S = cellfun(@(x) zscore(sum((x(:,seedmask).*(seedm)'),2)./sum(seedm.^2)), tcvox, 'un', 0);
    elseif SignMat(2)
        errordlg('PROBLEM WITH SIGN MAT !!!');
        % Computation of the spatially averaged seed signal
%         S = cellfun(@(x) (-1)*zscore(mean(x(:,seedmask).*(seeda(seedmask))',2)), tcvox, 'un', 0);
    else
        errordlg('PROBLEM WITH SIGN MAT !!!');
    end
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
