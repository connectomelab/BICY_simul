% bidirectional connectivity in serial vs. parallel
clear all
load adjacencyMAT_discrete_parallel.mat


percent_bi_dp=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_dp=mat_dp{i,j}{k};
            

             percent_bi_dp(i,j,k)=nnz(temp_dp(temp_dp==temp_dp'))/nnz(temp_dp);
           

        end
    end
end
percent_bi_dpc=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_dpc=mat_dpc{i,j}{k};
%             temp_dpcc=mat_dpcc{i,j}{k};
%             bicount=0;
             percent_bi_dpc(i,j,k)=nnz(temp_dpc(temp_dpc==temp_dpc'))/nnz(temp_dpc);
            % percent_bi_dpcc(i,j,k)=nnz(temp_dpcc(temp_dpcc==temp_dpcc'))/nnz(temp_dpcc);

            
        end
    end
end


save bi_dp.mat percent_bi_dp percent_bi_dpc
%%

load adjacencyMAT_discrete_serial.mat  

percent_bi_ds=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_ds=mat_ds{i,j}{k};
%             temp_dsc=mat_dsc{i,j}{k};
%             bicount=0;
             percent_bi_ds(i,j,k)=nnz(temp_ds(temp_ds==temp_ds'))/nnz(temp_ds);
            % percent_bi_dsc(i,j,k)=nnz(temp_dsc(temp_dsc==temp_dsc'))/nnz(temp_dsc);

            
        end
    end
end
percent_bi_dsc=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_dsc=mat_dsc{i,j}{k};

             percent_bi_dsc(i,j,k)=nnz(temp_dsc(temp_dsc==temp_dsc'))/nnz(temp_dsc);
           
            
        end
    end
end
save bi_ds.mat percent_bi_ds percent_bi_dsc

%%
clear all
load adjacencyMAT_ray_serial.mat

percent_bi_rs=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_rs=mat_rs{i,j}{k};
%             temp_rsc=mat_rsc{i,j}{k};
%             bicount=0;
             percent_bi_rs(i,j,k)=nnz(temp_rs(temp_rs==temp_rs'))/nnz(temp_rs);
            % percent_bi_rsc(i,j,k)=nnz(temp_rsc(temp_rsc==temp_rsc'))/nnz(temp_rsc);

            
        end
    end
end
percent_bi_rsc=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_rsc=mat_rsc{i,j}{k};

             percent_bi_rsc(i,j,k)=nnz(temp_rsc(temp_rsc==temp_rsc'))/nnz(temp_rsc);
           
            
        end
    end
end
save bi_rs.mat percent_bi_rs percent_bi_rsc

%%
clear all
load adjacencyMAT_ray_parallel.mat
percent_bi_rp=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_rp=mat_rp{i,j}{k};
%             temp_rpc=mat_rpc{i,j}{k};
%             bicount=0;
             percent_bi_rp(i,j,k)=nnz(temp_rp(temp_rp==temp_rp'))/nnz(temp_rp);
            % percent_bi_rpc(i,j,k)=nnz(temp_rpc(temp_rpc==temp_rpc'))/nnz(temp_rpc);

            
        end
    end
end

percent_bi_rpc=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_rpc=mat_rpc{i,j}{k};

             percent_bi_rpc(i,j,k)=nnz(temp_rpc(temp_rpc==temp_rpc'))/nnz(temp_rpc);
           
            
        end
    end
end

save bi_rp.mat percent_bi_rp percent_bi_rpc

%% hypothesis: bidirectional connection would be more frequent in parallel growing than serial growing

load bi_ds
load bi_dp
load bi_rs
load bi_rp

% with out limiting incoming connections, serial and parallel schemes did
% not differ in the number of bidirectional connections but when
% competition was imposed, more bidirectional connections existed in
% parallel growing.

close all
figure('Position',[0 0 1500 300])
for i=1:3
    
    for j=1:4
        s1=percent_bi_dsc(i,j,:);
        s2=percent_bi_dpc(i,j,:);
        [ h p]=ttest(s1, s2, 0.05,'left');
        hh(i,j)=h;
        pp(i,j)=p;
        
        subplot(1,4,j),
        boxplot([squeeze(squeeze(s1));squeeze(squeeze(s2))],[repmat(0,50,1);repmat(1,50,1)],'labels',{'Serial','Parallel'})
        
    end
    set(gcf, 'Color','w')
    print(sprintf('%s_%d.ai','figure4D',i),'-dpsc2')
end

close all
for i=1:3
   %figure('Position',[0 0 1500 300]);
   %birsc=[];birpc=[];
    for j=1:4
        s1=squeeze(percent_bi_rsc(i,j,:));
        s2=squeeze(percent_bi_rpc(i,j,:));
        %[ h p]=ttest(s1, s2, 0.05,'left');
        [ h p]=signrank(s1, s2, 0.05,'tail','left');
        hh(i,j)=h;
        pp(i,j)=p;
%         subplot(1,4,j),
%         boxplot([squeeze(squeeze(s1));squeeze(squeeze(s2))],[repmat(0,50,1);repmat(1,50,1)],'labels',{'Serial','Parallel'})
%         birsc=[birsc s1];
%         birpc=[birpc s2];
    end
%     set(gcf, 'Color','w')
%     print(sprintf('%s_%d.ai','figure_ray4D',i),'-dpsc2')
end
         

% discrete : less consistent but in general discrete seial<discrete
% parallel in percentage of bidirectional connections after bonferroni
% correction (alpha 0.05/12=0.0042) 8/12, 66% 
% hh =
% 
%      1     1     0     1
%      0     1     1     1
%      1     1     1     1
% 
% 
% 
% pp =
% 
%     0.0430    0.0005    0.0751    0.0025
%     0.0811    0.0000    0.0000    0.0001
%     0.0000    0.0000    0.0001    0.0051            

% ray: In all conditions, parallel growing produce more bidirectional
% connections than serial growing except one case with 1000 neurons and 8
% incoming connections allowed
% hh =
% 
%      1     1     1     1
%      1     1     1     1
%      1     1     1     1
% 
% pp
% 
% pp =
% 
%    1.0e-13 *
% 
%     0.0000    0.0000    0.0000    0.2668
%     0.0000    0.0000    0.0000    0.0000
%     0.0000    0.0000    0.0000    0.0000
            
 
%%
clear all
tic
load rewired_rpc.mat

percent_bi_rpc=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_rpc=rewired_rpc{i,j}{k};

             percent_bi_rpc(i,j,k)=nnz(temp_rpc(temp_rpc==temp_rpc'))/nnz(temp_rpc);
           
            
        end
    end
end
toc

rewired_bi_rpc=percent_bi_rpc;
save bi_rewire_rpc.mat rewired_bi_rpc
clear all
load rewired_rsc.mat
tic
percent_bi_rsc=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_rsc=rewired_rsc{i,j}{k};

             percent_bi_rsc(i,j,k)=nnz(temp_rsc(temp_rsc==temp_rsc'))/nnz(temp_rsc);
           
            
        end
    end
end

rewired_bi_rsc=percent_bi_rsc;

save bi_rewire_rsc.mat rewired_bi_rsc
toc
tic
clear all
load rewired_dsc.mat

percent_bi_dsc=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_dsc=rewired_dsc{i,j}{k};

             percent_bi_dsc(i,j,k)=nnz(temp_dsc(temp_dsc==temp_dsc'))/nnz(temp_dsc);
           
            
        end
    end
end

rewired_bi_dsc=percent_bi_dsc;
save bi_rewire_dsc.mat rewired_bi_dsc
toc
tic
clear all
load rewired_dpc.mat

percent_bi_dpc=zeros(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            temp_dpc=rewired_dpc{i,j}{k};

             percent_bi_dpc(i,j,k)=nnz(temp_dpc(temp_dpc==temp_dpc'))/nnz(temp_dpc);
           
            
        end
    end
end

rewired_bi_dpc=percent_bi_dpc;
save bi_rewire_dpc.mat rewired_bi_dpc

toc
%% comparison between generated networks and rewired networks
% mean(percent_bi_rsc./rewired_bi_rsc,3)
% 
% ans =
% 
%        Inf       Inf   16.4143   10.1393
%        Inf       Inf   13.8434   10.9903
%        Inf       Inf   17.8061   11.6817
% 
% mean(percent_bi_rpc./rewired_bi_rpc,3)
% 
% ans =
% 
%        Inf       Inf   21.0073   11.3533
%        Inf       Inf   32.6273   11.6146
%        Inf       Inf   42.5981   14.4102
% 
% mean(percent_bi_dpc./rewired_bi_dpc,3)
% 
% ans =
% 
%        Inf       Inf       Inf       Inf
%        Inf       Inf       Inf   10.4147
%        Inf       Inf       Inf   10.1932
% 
% mean(percent_bi_dsc./rewired_bi_dsc,3)
% 
% ans =
% 
%        NaN       Inf       Inf    9.6211
%        Inf       Inf       Inf   10.0571
%        Inf       Inf       Inf    9.9175

load('bi_rewire_dpc.mat')
load('bi_rewire_dsc.mat')
load('bi_rewire_rpc.mat')
load('bi_rewire_rsc.mat')
close all
hh=[];pp=[];
for i=1:3
    figure('Position',[0 0 1500 300])
    for j=1:4
        s1=squeeze(percent_bi_dsc(i,j,:));
        s2=squeeze(rewired_bi_dsc(i,j,:));
        %[ h, p]=ttest(s1, s2, 0.05,'right');
        [ h p]=signrank(s1, s2, 0.05,'tail','right');
        hh(i,j)=h;
        pp(i,j)=p;
        
        %subplot(1,4,j),
        boxplot([squeeze(squeeze(s1));squeeze(squeeze(s2))],[repmat(0,50,1);repmat(1,50,1)],'labels',{'Serial','Parallel'})
        
    end
    set(gcf, 'Color','w')
    %print(sprintf('%s_%d.ai','figure4D',i),'-dpsc2')
end

close all
for i=1:3
    %figure('Position',[0 0 1500 300])
    for j=1:4
        s1=percent_bi_rsc(i,j,:);
        s2=rewired_bi_rsc(i,j,:);
        [ h p]=ttest(s1, s2, 0.05,'left');
        hh(i,j)=h;
        pp(i,j)=p;
        
        %subplot(1,4,j),
        %boxplot([squeeze(squeeze(s1));squeeze(squeeze(s2))],[repmat(0,50,1);repmat(1,50,1)],'labels',{'Serial','Parallel'})
        
    end
    %set(gcf, 'Color','w')
    %print(sprintf('%s_%d.ai','figure4D',i),'-dpsc2')
end