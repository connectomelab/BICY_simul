clear all
load adjacencyMAT_discrete_parallel.mat

dist_bi_dp=cell(3,4,50);
for i=1:3
    
       
    for j=1:4
        
        for k=1:50
             ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
           
            load (['number_neuron_increase_size' num2str(ss(j)) '.mat'], 'dataset') 
            tmp=dataset{k,2*i-1};
             distance=dist(tmp');
            temp_dp=mat_dp{i,j}{k};
            

             %percent_bi_dp(i,j,k)=nnz(temp_dp(temp_dp==temp_dp'))/nnz(temp_dp);
             biandzero=find(temp_dp==temp_dp');
             nozero=find(temp_dp);
             purebi=intersect(biandzero, nozero);
             dist_bi_dp{i,j,k}=distance(purebi);
             

        end
        clear biandzero nozero purebi
    end
end

dist_bi_dpc=cell(3,4,50);
for i=1:3
    
    for j=1:4
        
        for k=1:50
            ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
           
            load (['number_neuron_increase_size' num2str(ss(j)) '.mat'], 'dataset') 
            tmp=dataset{k,2*i-1};
             distance=dist(tmp');
            temp_dpc=mat_dpc{i,j}{k};
            

             %percent_bi_dpc(i,j,k)=nnz(temp_dpc(temp_dpc==temp_dpc'))/nnz(temp_dpc);
             biandzero=find(temp_dpc==temp_dpc');
             nozero=find(temp_dpc);
             purebi=intersect(biandzero, nozero);
             dist_bi_dpc{i,j,k}=distance(purebi);
             

        end
    end
end

save bidirec_dist_dp.mat dist_bi_dp dist_bi_dpc
clear all

%%
load adjacencyMAT_discrete_serial.mat
dist_bi_ds=cell(3,4,50);
for i=1:3
   
    for j=1:4
        
        for k=1:50
            ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
           
            load (['number_neuron_increase_size' num2str(ss(j)) '.mat'], 'dataset') 
            tmp=dataset{k,2*i-1};
             distance=dist(tmp');
            temp_ds=mat_ds{i,j}{k};
            

             %percent_bi_ds(i,j,k)=nnz(temp_ds(temp_ds==temp_ds'))/nnz(temp_ds);
             biandzero=find(temp_ds==temp_ds');
             nozero=find(temp_ds);
             purebi=intersect(biandzero, nozero);
             dist_bi_ds{i,j,k}=distance(purebi);
             

        end
    end
end
dist_bi_dsc=cell(3,4,50);
for i=1:3
    
    for j=1:4
        
        for k=1:50
            ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
           
            load (['number_neuron_increase_size' num2str(ss(j)) '.mat'], 'dataset') 
            tmp=dataset{k,2*i-1};
             distance=dist(tmp');
            temp_dsc=mat_dsc{i,j}{k};
            

             %percent_bi_dsc(i,j,k)=nnz(temp_dsc(temp_dsc==temp_dsc'))/nnz(temp_dsc);
             biandzero=find(temp_dsc==temp_dsc');
             nozero=find(temp_dsc);
             purebi=intersect(biandzero, nozero);
             dist_bi_dsc{i,j,k}=distance(purebi);
             

        end
    end
end
save bidirec_dist_ds.mat dist_bi_ds dist_bi_dsc


clear all
%%
clear all

load adjacencyMAT_ray_serial.mat
dist_bi_rs=cell(3,4,50);
for i=1:3
   
    for j=1:4
        
        for k=1:50
            ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
           
            load (['number_neuron_increase_size' num2str(ss(j)) '.mat'], 'dataset') 
            tmp=dataset{k,2*i-1};
             distance=dist(tmp');
            temp_rs=mat_rs{i,j}{k};
            

             %percent_bi_rs(i,j,k)=nnz(temp_rs(temp_rs==temp_rs'))/nnz(temp_rs);
             biandzero=find(temp_rs==temp_rs');
             nozero=find(temp_rs);
             purebi=intersect(biandzero, nozero);
             dist_bi_rs{i,j,k}=distance(purebi);
             

        end
        
    end
end
%%
dist_bi_rsc=cell(3,4,50);
for i=1:3
    
    for j=1:4
        
        for k=1:50
            ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
           
            load (['number_neuron_increase_size' num2str(ss(j)) '.mat'], 'dataset') 
            tmp=dataset{k,2*i-1};
             distance=dist(tmp');
            temp_rsc=mat_rsc{i,j}{k};
            

             %percent_bi_rsc(i,j,k)=nnz(temp_rsc(temp_rsc==temp_rsc'))/nnz(temp_rsc);
             biandzero=find(temp_rsc==temp_rsc');
             nozero=find(temp_rsc);
             purebi=intersect(biandzero, nozero);
             dist_bi_rsc{i,j,k}=distance(purebi);
             

        end
    end
end
save bidirec_dist_rs.mat dist_bi_rs dist_bi_rsc

clear all
%%
load adjacencyMAT_ray_parallel.mat

dist_bi_rp=cell(3,4,50);
for i=1:3
    
    for j=1:4
        
        for k=1:50
            ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
           
            load (['number_neuron_increase_size' num2str(ss(j)) '.mat'], 'dataset') 
            tmp=dataset{k,2*i-1};
             distance=dist(tmp');
            temp_rp=mat_rp{i,j}{k};
            

             %percent_bi_rp(i,j,k)=nnz(temp_rp(temp_rp==temp_rp'))/nnz(temp_rp);
             biandzero=find(temp_rp==temp_rp');
             nozero=find(temp_rp);
             purebi=intersect(biandzero, nozero);
             dist_bi_rp{i,j,k}=distance(purebi);
             

        end
    end
end

dist_bi_rpc=cell(3,4,50);
for i=1:3
    
    for j=1:4
        
        for k=1:50
            ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
           
            load (['number_neuron_increase_size' num2str(ss(j)) '.mat'], 'dataset') 
            tmp=dataset{k,2*i-1};
             distance=dist(tmp');
            temp_rpc=mat_rpc{i,j}{k};
            

             %percent_bi_rpc(i,j,k)=nnz(temp_rpc(temp_rpc==temp_rpc'))/nnz(temp_rpc);
             biandzero=find(temp_rpc==temp_rpc');
             nozero=find(temp_rpc);
             purebi=intersect(biandzero, nozero);
             dist_bi_rpc{i,j,k}=distance(purebi);
             

        end
    end
end
save bidirec_dist_rp.mat dist_bi_rp dist_bi_rpc
%%
clear all
load bidirec_dist_ds
load bidirec_dist_dp
load bidirec_dist_rs
load bidirec_dist_rp
% %%
% for i=1%:3
%     for j=1%:4
%         
%         for k=1:50
%             maxbirsc(k)=min(dist_bi_rsc{i,j,k});
%             maxbirpc(k)=min(dist_bi_rpc{i,j,k});
%             
%         end
%         maxbidist_rsc(i,j)=max(maxbirsc);
%         maxbidist_rpc(i,j)=max(maxbirpc);
%     end
% end
% example 1,1 case size 1 1000 neurons

load('adjacencyMAT_discrete_parallel.mat', 'ss')
load ('ray_serial_stat.mat','w_rsc')
load ('ray_parallel_stat.mat','w_rpc')
tic
for i=3
    for j=2%:4
        load (['number_neuron_increase_size' num2str(ss(j)) '.mat'], 'dataset')  
        for k=1:50
            ersc=2:1:18;erpc=2:1:13;% applicable for {1,1} case
            tempdata=dataset{k,2*i-1};
            tempdist=dist(tempdata');
            for nbin=2:18
                
                %distbirsc(k,:)=histc(dist_bi_rsc{i,j,k},ersc)./length(dist_bi_rsc);
                %distbirpc(k,:)=histc(dist_bi_rpc{i,j,k},erpc)./length(dist_bi_rsc);
                nomehist_rsc(1,nbin-1)= sum(dist_bi_rsc{i,j,k}<=nbin & dist_bi_rsc{i,j,k}>nbin-1);
                
                denomhist(1,nbin-1)=sum(tempdist(:)<=nbin & tempdist(:)>nbin-1);
            end
            distbirsc(k,:)=nomehist_rsc./denomhist;
            denomhist=[];
            for xx=2:13
              nomehist_rpc(1,xx-1)= sum(dist_bi_rpc{i,j,k}<=xx & dist_bi_rpc{i,j,k}>xx-1);  
              denomhist(1,xx-1)=sum(tempdist(:)<=xx & tempdist(:)>xx-1);
            end
            distbirpc(k,:)=nomehist_rpc./denomhist;
            
        end
    end
end
toc

mbirsc=mean(distbirsc);
% figure;loglog(ersc,mbirsc,'ro')
% figure;semilogy(ersc,mbirsc,'bo')

sdrsc=std(distbirsc);
mbirpc=mean(distbirpc);
sdrpc=std(distbirpc);

%%
close all
figure('Position',[0 0 600 500])

errorbar(ersc,mbirsc,sdrsc ,'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0],...
    'MarkerSize',3,...
    'Marker','o', 'LineStyle','-', 'Color', 'r', 'LineWidth',1.5)
hold on
errorbar(erpc,mbirpc,sdrpc ,'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor',[0 0 1],...
    'MarkerSize',3,...
    'Marker','o', 'LineStyle','-', 'Color', 'b', 'LineWidth',1.5)
set(gcf, 'Color','w')
xlabel('Distance between neurons (unit)')
ylabel('Proportion of bidirectional connections')
print ('bidirection.ai','-dpsc2')

%%
    figure('Position',[0 0 600 500])

loglog(ersc,mbirsc,'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0],...
    'MarkerSize',3,...
    'Marker','o', 'LineStyle','none')
hold on
loglog(erpc,mbirpc,'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor',[0 0 1],...
    'MarkerSize',3,...
    'Marker','o', 'LineStyle','none')


hold on
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Robust = 'LAR';
opts.Upper = [Inf Inf];


nzx=1:length(ersc(mbirsc>0));
[a b]=fit(log(ersc(mbirsc>0))',log(nonzeros(mbirsc)),ft,opts);
%loglog(nzx,nzx.^a.p1*exp(a.p2),'r', 'LineWidth',1.5)
loglog(ersc(mbirsc>0),ersc(mbirsc>0).^a.p1*exp(a.p2),'r', 'LineWidth',1.5)

nzx=1:length(erpc(mbirpc>0));
[c d]=fit(log(erpc(mbirpc>0))',log(nonzeros(mbirpc)),ft,opts);
loglog(erpc(mbirpc>0),erpc(mbirpc>0).^c.p1*exp(c.p2),'b', 'LineWidth',1.5)
xlabel('Distance between neurons (unit)')
ylabel('Connection probability')
set(gcf,'Color','w')
gtext(sprintf('%s%5.2f%s%s%5.2f','y = ',a.p1,' * x  ',a.p2))
gtext(sprintf('%s%5.2f%s%s%5.2f','y = ',c.p1,' * x  ',c.p2))
print('bidirecloglog.ai','-dpsc2')

%%
close all
d=2:60; r=ss;
figure('Position',[0 0 600 500])

for i=1:4
y=(asin((1+r(i))./d));
y=(y./sum(y));


plot(d,y,'o-','MarkerEdgeColor',[0.25*i,0 0],'Color','r')
hold on
end
xlabel('Distance between neurons (unit)')
ylabel('Connection probability')
set(gcf,'Color','w')
%print('arcsin1.ai','-dpsc2')
hold on

d=2:60; r=ss;
%figure('Position',[0 0 600 500])

for i=1:4
y=(asin((1+r(i))./d));
y=(y./sum(y)).^2;


plot(d,y,'o-','MarkerEdgeColor',[0 0 0.25*i])
hold on
end
xlabel('Distance between neurons (unit)')
ylabel('Connection probability')
set(gcf,'Color','w')
%print('arcsin_bi1.ai','-dpsc2')
%print('arcsin_overall_bi.ai','-dpsc2')
%%
%test pareto q-q
r=exprnd(1,1,length(y));
figure;qqplot(log(y),r)
%%
%%
ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
d=2:60;
  figure('Position',[0 0 600 500])
 for i=2%:4
y=(asin((1+ss(i))./d));
y=(y./sum(y));

loglog(d,y,'o','MarkerEdgeColor',[0.25*i,0 0],'Color','r')
hold on
x=(1+ss(i))./d;
y2=x;%+(1/6)*x.^3;%+(3/40)*x.^5+(55/112)*x.^7;
y2=(y2./sum(y2));
loglog(d,y2,'s')
figure;plot(d,y,'ro')
hold on 
plot(d,y2,'bs')
 end

 
 
 
xlabel('Distance between neurons (unit)')
ylabel('Connection probability')
set(gcf,'Color','w')
%print('arcsin2.ai','-dpsc2')
 hold on
 for i=1:4
y=(asin((1+ss(i))./d));
%y=(y./sum(y)).^2;

%loglog(d,y,'o','MarkerEdgey

loglog(d,y,'o','MarkerEdgeColor',[0 0 0.25*i])
hold on
end
xlabel('Distance between neurons (unit)')
ylabel('Connection probability')
set(gcf,'Color','w')
%print('arcsin_bi2.ai','-dpsc2')
%print('arcsin_bi_loglog.ai','-dpsc2')