function [matrix,positions,connection,clength,connectible,fc_abs,w,fc]=devolution3d_ng_ray_serialN(n,flagtaken,flagsize,flagdistance,mplaces,range)%,msynapses)
% function [connection,fc_abs,w,fc]=devolution3d_ng(n,flagtaken,flagsize,flagdistance);
% wire n nodes by establishing a connection
% nodes of size 1 are randomly positioned on a
% 34x34x34 grid
% if flag is set, the development is static,
% that means that all nodes are positioned before
% connections are established
% output: matrix: n x n
%         positions: 3D coordinates (range 0..limit)
%         w: metric lengths of all connections
%         fc: filling fraction (see Stepanyants2002)
%		  connectible: indices of connectible neurons
%	      connection: indices of connected neurons
%         clengths
%         fc_abs: the absolute number of connectible spots

% Author: Sol Lim, Marcus Kaiser  Date:  
global positions para cellsize
 cellsize=flagsize;
% initialisation of node positions
LIMIT=34; % 39304 positions, about 4 times as many as for 2D 
          % in order have a comparable chance for establishing 
          % a cell position without cell overlap
MINSIZE=0; % was 2
MAXPLACES=mplaces; % was 1
% MAXSYNAPSES=msynapses; % was 1
fc_abs=0;
% RandStream.setDefaultStream ...
%      (RandStream('mt19937ar','seed',sum(100*clock)));

limit=LIMIT;
% [finalpos,l_index,n_step,t,check,cellsize]=neuron_positions(n,flagsize,limit);
% positions=finalpos;

% if flagsize == 0
%   cellsize = ones(n,1); % unit cell size
% else
%   cellsize = 2*rand(n,1); % random cell size in[0,2], average 1
% end
% 
% positions = zeros(n,3);
% positions(1,:) = [LIMIT/2, LIMIT/2, LIMIT/2];
% for i=2:n
%     positions(i,:) = LIMIT*rand(1,3);
%     d=dist(positions(1:i,:)');
%     d2=d(i,1:i-1);
%     [d3,j] = min(d2); 
%     while (d3(1)-cellsize(i)-cellsize(j(1))) < MINSIZE
%         positions(i,:) = LIMIT*rand(1,3);
%         d=dist(positions(1:i,:)');
%         d2=d(i,1:i-1);
%         [d3,j] = min(d2); 
%     end;
% end;
 


matrix = zeros(n);
free_places = MAXPLACES * ones(n,1); % free places for synapses
%established_synapses = zeros(n,1); % number of outgoing synapses of a neuron
%finished = zeros(n,1); % neuron already has an incoming connection
connection=cell(n,1);% and can not establish another oned=dist(positions');
clength=cell(n,1);
connectible=cell(n,1);
d=dist(positions');

for i=1:n

        
       %phi=direc(i,1); theta=direc(i,2);
        %ray = positions(i,:) + [sin(theta)*cos (phi), sin(theta)*sin(phi), cos(theta)]; %allocate a direction for the ray
        z=para(:,1);theta=para(:,2);r=para(:,3); % for comparison we fix all directions for all algorithms 
        ray = [r(i)*cos(theta(i)), r(i)*sin(theta(i)), z(i)]; % for comparison
% allocate random directions       
%         z=-1+2*rand(1);
%         theta=2*pi*rand(1);
%         r=sqrt(1-z^2);
%         ray=[r*cos(theta), r*sin(theta), z];
         
        %nray=ray/norm(ray); %notmal vector of the given ray
        
        rv=positions-repmat(positions(i,:),n,1); 
        dis=d(:,i); 
        cosine=(rv*ray')./dis;       %    n*1
        pcosine=find(cosine>=0); 
        
        if flagdistance==0
            d_ray=dis(pcosine).*sqrt(1-cosine(pcosine).^2); %the distance btw ray and centers of spheres 
        else
            tempcellsize=cellsize*ones(n,1);
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
%                [unidis,ind]=unique(temp2(:,3));
%                connect=temp2(ind,:) ;
               
%                connect_ind=connect(:,1);
			% connect=temp2;
               connect_ind=temp2(:,1);
				connectible{i,1}=connect_ind; % the indices of connectible neurons with neuron i
            end;   
        
        %fc_abs_ori =fc_abs_ori+nnz(connect_ind); %original %fc + 1;  
        
            for j=1:length(connect_ind) %(1:MAXSYNAPSES)
                if connect_ind==0
                   connect=[0 0 0];
                   break;
                else fc_abs =fc_abs+1;
                %elseif
                     if (flagtaken == 0) ||(free_places(connect_ind(j)) ~= 0) 
                     matrix(i,connect_ind(j)) = 1;
                     %matrix(nn,i) = 1;
                     free_places(connect_ind(j)) =free_places(connect_ind(j))- 1;
                     %established_synapses(i) = established_synapses(i) + 1;
                     %fc_abs =fc_abs-(length(connect_ind)-j-1);
%                         if established_synapses(i) == MAXSYNAPSES
%                             break;
%                         end; % if 
                     end;
                 %break;
                 end;    
                
            end; 
        count= nnz(matrix(i,:));
        connection(i)={connect(1:count,1)}; %the index of connected neurons(arrival)
        clength(i)={connect(1:count,2)};% the distance btw the connected neurons, connection length
     
         
                 
end; % for i




w = nonzeros(matrix.*d);
fc = length(w) / fc_abs;
% 
% fprintf('Number of edges: %d\n', length(w));
% fprintf('Edge density: %f\n', length(w)/(n*(n-1)));
% fprintf('Percentage of filled slots: %f\n', n/(LIMIT*LIMIT*LIMIT));
%fprintf('Exponential parameter mu: %f\n', expfit(w));

%l=wirelengths2(celegans131matrix,celegans131positions);[y,x]=hist(l);figure(1);plot(x,y);figure(2);plot(x,length(l)*exppdf(x,expfit(l)));expfit(l)
% l=[];for i=1:20;[matrix,positions,w]=devolution(400,1);l=[l w];end;[y,x]=hist(nonzeros(l));cftool(x,y)
return


