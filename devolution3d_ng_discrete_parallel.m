function [matrix,positions,connection,connectible,fc_abs,w,fc]=devolution3d_ng_discrete_parallel(n,limit,flagtaken,flagsize,cellradius,flagdistance,mplaces,range)
%function [matrix,positions,connection,connectible,fc_abs,w,fc]=devolution3d_ng_discrete_parallel(n,limit,flagtaken,flagsize,flagdistance,mplaces,range)
% wire n nodes by establishing a connection
% neurons(nodes) are positioned randomly in the given space by
% neuron_positions.m 
% all nodes are positioned before connections are established
% all neurons grow their axons simultaneously
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
% Author: Sol Lim, this was modified version of Marcus Kaiser's devolution3d_ng.m  


RandStream.setGlobalStream ...
     (RandStream('mt19937ar','seed',sum(100*clock)));


MAXPLACES=mplaces; % the maximum number of incoming connections for a neuron to accomodate
size=cellradius;
[finalpos,check,cellsize]=neuron_positions(n,flagsize,size,limit);%placing neurons
positions=finalpos;

matrix = zeros(n);count=0;
free_places = MAXPLACES * ones(n,1); % free places for synapses

for i=1:n
        % allocate random directions       
        z=-1+2*rand(1);
        theta=2*pi*rand(1);
        r=sqrt(1-z^2);
        ray=[r*cos(theta), r*sin(theta), z];
        pos = positions(i,:);
    
    while withingrid(pos,limit)     
                
        pos = pos + ray;
      
        
           nn = nearestnode(positions,multipos(pos,n),i,cellsize,flagdistance,range);
        
        
          if (nn(1) ~= i) && (nn(1) ~= 0)
           
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                count=count+1;
            
                connectible(count,1)={[i nn]};

           end; % if
    end; % while
end; % for i

fc_abs=0;
candidates=cell2mat(connectible);

sorted_can=sortrows(candidates,3);
size_can=size(candidates);

for j=1:size_can(1)
    
        fc_abs=fc_abs+1;
         if (flagtaken==0) || (free_places(sorted_can(j,2)))~=0
            matrix(sorted_can(j,1),sorted_can(j,2))=1;
            connection{j,1}=sorted_can(j,:);
            free_places(sorted_can(j,2))=free_places(sorted_can(j,2))-1;
            
         end

end
            
 
w = nonzeros(matrix.*dist(positions'));
fc = length(w) / fc_abs;
 
% fprintf('Number of edges: %d\n', length(w));
% fprintf('Edge density: %f\n', length(w)/(n*(n-1)));
% fprintf('Percentage of filled slots: %f\n', n/(LIMIT*LIMIT*LIMIT));
% %fprintf('Exponential parameter mu: %f\n', expfit(w));
% 
%l=wirelengths2(celegans131matrix,celegans131positions);[y,x]=hist(l);figure(1);plot(x,y);figure(2);plot(x,length(l)*exppdf(x,expfit(l)));expfit(l)
% l=[];for i=1:20;[matrix,positions,w]=devolution(400,1);l=[l w];end;[y,x]=hist(nonzeros(l));cftool(x,y)
return




function nn=withingrid(pos,limit)
nn = 0;
if (pos(1)>0) & (pos(2)>0) & (pos(3)>0) & (pos(1)<limit) & (pos(2)<limit) & (pos(3)<limit)
   nn = 1;
end;
return;

%Euclidean distance
function nn = nearestnode(positions,pos, i,cellsize,flagdistance,range)
p2 = abs(positions-pos);
%a = min((p2'));
d=sqrt(sum(p2.^2,2));
[val,loc] = min((d));
if flagdistance==1
    val_givensize=val-cellsize;
else
    val_givensize=val;
end
if (loc ~= i) & (val_givensize <=range)
    nn =[ loc val_givensize];
else
    nn =[ 0 0];
end;
return;

function mpos = multipos(pos,n)
mpos=[];
for i=1:n
    mpos = [mpos; pos];
end;