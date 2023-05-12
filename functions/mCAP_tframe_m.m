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
