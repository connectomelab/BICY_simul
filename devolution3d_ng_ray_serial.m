function [matrix,positions,connection,clength,connectible,fc_abs,w,fc]=devolution3d_ng_ray_serial(n,limit,flagtaken,flagsize,cellradius,flagdistance,mplaces,range);
% function[matrix,positions,connection,clength,connectible,fc_abs,w,fc]=devolution3d_ng_ray_serial(n,limit,flagtaken,flagsize,cellradius,flagdistance,mplaces,range);
% wire n nodes by establishing a connection
% neurons(nodes) are positioned randomly in the given space by
% neuron_positions.m 
% all nodes are positioned before connections are established
% neurons take turns to grow its axons serially; when a neuron finishes
% growing as the growth cone goes beyond the embedding space, the next
% neuron starts growing its axons. The order of growth is random.
% all the possible synapses that a neuron can establish are detected using
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




RandStream.setGlobalStream ...
      (RandStream('mt19937ar','seed',sum(100*clock)));
  
MAXPLACES=mplaces; % the maximum number of incoming connections for a neuron to accomodate

size=cellradius;
[finalpos,check,cellsize]=neuron_positions(n,flagsize,size,limit);%placing neurons
positions=finalpos;

matrix = zeros(n);
free_places = MAXPLACES * ones(n,1); % free places for synapses
connection=cell(n,1);% and can not establish another oned=dist(positions');
clength=cell(n,1);
connectible=cell(n,1);
d=dist(positions');
fc_abs=0;

for i=1:n
        %allocate random directions   
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
            tempcellsize=cellsize;%*ones(n,1); % dimension of cellsize is herited from neuron_position here
            d_ray=dis(pcosine).*sqrt(1-cosine(pcosine).^2)-tempcellsize(pcosine);% distance given cellsize 
        end
        temp=[pcosine d_ray];
        connect= d_ray<=range;%finding neurons within connectable range
        connect=temp(connect,:); % connectable neurons and their original position indicies
            if isempty(connect)==1
                connect_ind=0;
               
            else 
           
               d_along=dis(connect(:,1)).*cosine(connect(:,1));
               temp1=[connect d_along]; % original indices and their distances from the active neuron making connections
               temp2=sortrows(temp1,3);
               connect_ind=temp2(:,1);
				connectible{i,1}=connect_ind; % the indices of connectible neurons with neuron i
            end;   
        
        
            for j=1:length(connect_ind) 
                if connect_ind==0
                   connect=[0 0 0];
                   break;
                else fc_abs =fc_abs+1;
                
                     if (flagtaken == 0) ||(free_places(connect_ind(j)) ~= 0) 
                     matrix(i,connect_ind(j)) = 1;
                     free_places(connect_ind(j)) =free_places(connect_ind(j))- 1;
                     end
                end    
                
            end 
        count= nnz(matrix(i,:));
        connection(i)={connect(1:count,1)}; %the index of connected neurons(arrival)
        clength(i)={connect(1:count,2)};% the distance btw the connected neurons, connection length
                 
end; % for i

w = nonzeros(matrix.*d);
fc = length(w) / fc_abs;
return


