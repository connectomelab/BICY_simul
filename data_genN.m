%% data set generation --parallell embedding
% the number of neurons increasing leading to the spatial density increase
RandStream.setGlobalStream ...
     (RandStream('mt19937ar','seed',sum(100*clock)));

 nn=1000:500:2000; dataset=cell(100,3);process=cell(100,3);n_steps=zeros(100,3);tt=zeros(100,3);checkk=zeros(100,3);
for i=1:3
    n=nn(i);size=0.5;limit=34;
    for j=1:100
        
        [finalpos,l_index,n_step,t,check]=neuron_positions(n,size,limit);
        dataset{j,i}=finalpos;
        process{j,i}=l_index;
        n_steps(j,i)=n_step;
        tt(j,i)=t;
        checkk(j,i)=check;
        
    end
end
save 1000_1500_2000_size0.5

 
  nn=1000:500:2000; dataset=cell(100,3);process=cell(100,3);n_steps=zeros(100,3);tt=zeros(100,3);checkk=zeros(100,3);
for i=1:3

    n=nn(i);size=2^(1/3)*0.4+0.1;limit=34; %size=2^(1/3)*0.4+0.1 size = 0.6040
    for j=1:100
        
        [finalpos,l_index,n_step,t,check]=neuron_positions(n,size,limit);
        dataset{j,i}=finalpos;
        process{j,i}=l_index;
        n_steps(j,i)=n_step;
        tt(j,i)=t;
        checkk(j,i)=check;
        
    end
end
save (['1000_1500_2000_size' num2str(size) '.mat'])


clear all
 
 nn=1000:500:2000; dataset=cell(100,3);process=cell(100,3);n_steps=zeros(100,3);tt=zeros(100,3);checkk=zeros(100,3);
for i=1:3

    n=nn(i);size=4^(1/3)*0.4+0.1;limit=34; % size=4^(1/3)*0.4+0.1 size =0.7350

    for j=1:100
        
        [finalpos,l_index,n_step,t,check]=neuron_positions(n,size,limit);
        dataset{j,i}=finalpos;
        process{j,i}=l_index;
        n_steps(j,i)=n_step;
        tt(j,i)=t;
        checkk(j,i)=check;
        
    end
end
save (['1000_1500_2000_size' num2str(size) '.mat'])

 nn=1000:500:2000; dataset=cell(100,3);process=cell(100,3);n_steps=zeros(100,3);tt=zeros(100,3);checkk=zeros(100,3);
for i=1:3
    n=nn(i);size=0.9;limit=34; %size=2*a+0.1 a=0.4
    for j=1:100
        
        [finalpos,l_index,n_step,t,check]=neuron_positions(n,size,limit);
        dataset{j,i}=finalpos;
        process{j,i}=l_index;
        n_steps(j,i)=n_step;
        tt(j,i)=t;
        checkk(j,i)=check;
        
    end
end
save (['1000_1500_2000_size' num2str(size) '.mat'])
