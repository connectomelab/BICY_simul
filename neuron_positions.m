function [finalpos,check,cellsize]=neuron_positions(n,flagsize,size,limit)
%function [finalpos,l_index,n_step,t,check,cellsize]=neuron_positions(n,flagsize,size,limit)
%place neurons in parallel to shorten processing time in the given 3D space without overlapping
%n : the number of neurons
%size : the size (radius) of an neuron which includes its soma and
%dendrites
%flagsize: if flag is set, the size each neuron is assigned
%randomly between 0 and the size
%limit: limit of the embedding cube
%output finalpos: the centre position of a neuron (x, y, z)
%       check: whether all neurons are placed in the given space without
%       overlapping
%       cellsize: n x 1 matrix storing radii of neurons
%Author: Sol Lim

RandStream.setGlobalStream ...
     (RandStream('mt19937ar','seed',sum(100*clock)));


minsize=0;LIMIT=limit;n_step=0;

    
   
if flagsize == 0
  cellsize = size*ones(n,1); % unit cell size
else
  cellsize = size*rand(n,1); % random cell size in[0,2], average 1
end

   
    positions=[LIMIT*rand(n,1) LIMIT*rand(n,1) LIMIT*rand(n,1)];
    positions(1,:)=[LIMIT/2, LIMIT/2, LIMIT/2];
    index=0;l_index=n;
    cellsizesum=repmat(cellsize,1,n)+repmat(cellsize',n,1);
    
    while isempty(index)==0 
        d=dist(positions');% distances btw neurons  
        index=find((triu(d,1)-triu(cellsizesum,1))<minsize);%to find overlapping neurons
        dpos=[mod(index,n) ceil(index./n)]; % to find indices of  overlapping neurons
        zero_mod=find(dpos(:,1)==0);
        dpos(zero_mod,1)=n;
        positions(max(dpos,[],2),:)=LIMIT*rand(length(dpos(:,1)),3);%to assign a new position for the relevant neuron with bigger index


        l_index=[l_index; length(index)]; % to record the number of overlapping neurons
        n_step=n_step+1; % the number of steps needed
        if isempty(index)==1
            break
            else index=0;
        end
    end    
       
       
       finalpos=positions;
       d=dist(finalpos');% distances btw neurons

       
       index=find((triu(d,1)-triu(cellsizesum,1))<minsize);%to find ove
       check=isempty(index);
%end
