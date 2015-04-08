function [matrix,positions,connection,connectible,fc_abs,w,fc]=PartialOverlap_ray_parallel_gp2(n,flagtaken,flagsize,flagdistance,mplaces,range,group,overlap)
%function [matrix,positions,connection,connectible,fc_abs,w,fc]=PartialOverlap_ray_parallel_gp2(n,flagtaken,flagsize,flagdistance,mplaces,range,group,overlap)
% wire n nodes by establishing a connection
% neurons(nodes) are positioned randomly in the given space by
% neuron_positions.m 
% all nodes are positioned before connections are established
% neurons are divided into 2 groups and they have overlapping time windows
% for axon growth; neurons within the same group start growing their axons
% at the same time and one group starts axon growth and 
% all the possible synapses that all neuron can establish are detected using
% the intersection between the 'ray' or the growth direction and the
% neurons in the vicinity considering the radii of neurons and the
% proxlimity rule(see below)
% n : the number of neurons
% flagsize: if flagsize is set, the size each neuron is assigned
% randomly between 0 and the size
% flagtaken: if flagtaken is set, neurons have limited number of incoming
% connections, leading to competition among neurons
% flagdistance: if flagdistance is set, the distance between the centers of
% neurons takes into account the radii of the neurons that are to be
% connected
% limit: limit of the embedding space
% mplaces: the number of incoming connections for a neuron to accomodate a
% synapses
% range: the proximity rule to connect a neuron; if the distance between
% the growth cone and a neuron within vicinity is less than or equal to
% this value, it is assumed that the growth cone would connect the neuron
% within the range.
% output: matrix: n x n adjacency matrix
%         positions: 3D coordinates (range 0..limit)
%         w: Euclidean distances of all connections
%         fc: filling fraction (see Stepanyants2002)
%		  connectible: indices of connectible neurons
%	      connection: indices of connected neurons
%         fc_abs: the absolute number of connectible spots
% Author: Sol Lim 

RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',sum(100*clock)));


MAXPLACES=mplaces; % was 1


size=cellradius;
[finalpos,check,cellsize]=neuron_positions(n,flagsize,size,limit); %placing neurons
positions=finalpos;

fc_abs=0;
matrix = zeros(n);
free_places = MAXPLACES * ones(n,1); % free places for synapses
d=dist(positions'); cind=cell(n,1);

for i=1:n
       
allocate random directions       
        z=-1+2*rand(1);
        theta=2*pi*rand(1);
        r=sqrt(1-z^2);
        ray=[r*cos(theta), r*sin(theta), z];
        
        rv=positions-repmat(positions(i,:),n,1); 
        dis=d(:,i); 
        cosine=(rv*ray')./dis;       %    n*1
        pcosine=find(cosine>=0); 
        
        if flagdistance==0
            d_ray=dis(pcosine).*sqrt(1-cosine(pcosine).^2); %the distance btw ray and centers of spheres
        else
            tempcellsize=cellsize*ones(n,1); % cellsize increased by 0.02
            d_ray=dis(pcosine).*sqrt(1-cosine(pcosine).^2)-tempcellsize(pcosine);% distance given cellsize
        end
        temp=[pcosine d_ray];
        connect= d_ray<=range;%finding neurons within connectable range
        connect=temp(connect,:); % connectable neurons and their original position indicies
        
        if isempty(connect)==1
            connect_ind=0;
            connect=zeros(1,3);
            
        else
            
            d_along=dis(connect(:,1)).*cosine(connect(:,1));
            temp1=[connect d_along]; % original indices and their distances from the active neuron making connections
            temp2=sortrows(temp1,3);
            connect=temp2;
            connect_ind=temp2(:,1);
        end;
        
        if i<=n/group;
            cind(i,1)={[i*ones(length(connect_ind),1) connect(1:length(connect_ind),:)]};% [from neuron x, to neuron y, the closest distance, the distance from neuron x to y (btw centers)
        else
            connect(1:length(connect_ind),3)=connect(1:length(connect_ind),3)+overlap; % implementing overlap
            cind(i,1)={[i*ones(length(connect_ind),1) connect(1:length(connect_ind),:)]};
        end
end        

   cind_mat=cell2mat(cind); % [from neuron x, to neuron y, the closest distance, the distance from neuron x to y (btw centers)  
   
   cind_sort=sortrows(cind_mat,4); % sort cind_mat according to the distance between connectible neurons
   
   connectible=cell(length(cind_sort) ,1);     
           for j=1:length(cind_sort) 

                if cind_sort(j,2)~=0 
                        fc_abs=fc_abs+1;
						connectible{j}=[cind_sort(j,1),cind_sort(j,2)];
                        if (flagtaken == 0) ||(free_places(cind_sort(j,2)) ~= 0)
							matrix(cind_sort(j,1),cind_sort(j,2)) = 1;
							connection{j}=[cind_sort(j,1),cind_sort(j,2)];
							free_places(cind_sort(j,2)) =free_places(cind_sort(j,2))- 1;
                        end    

                 end; % if                 
            end;  



w = nonzeros(matrix.*d);
fc = length(w) / fc_abs;

return



