function [emptk,maxsumd,frame_assign,sumda]=mCAP_tframe_m(frame,Cma)

k=size(Cma,1);
f= size(frame,2);

for im = 1 : f
   for ic = 1 : k
        dis_mat(ic,1)=pdist2(frame(:,im)',Cma(ic,:),'correlation'); %%% same way of calculating distance in Kmeans        
   end
   
    targetcap=find(dis_mat==min(dis_mat));
    frame_assign(im,1)=targetcap(1);
    dis_frame_assign(im,1)=min(dis_mat);

end


for ic = 1: k 
    sumda(ic,1)=sum(dis_frame_assign(frame_assign==ic));
end
emptk=k-size(unique(frame_assign),1);
maxsumd=max(sumda);
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
