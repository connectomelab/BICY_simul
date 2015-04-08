clear all
rs_out=cell(3,4,50);rs_in=cell(3,4,50);
rsc_out=cell(3,4,50);rsc_in=cell(3,4,50);rs_deg=cell(3,4,50);rsc_deg=cell(3,4,50);
rp_out=cell(3,4,50);rp_in=cell(3,4,50);
rpc_out=cell(3,4,50);rpc_in=cell(3,4,50);rp_deg=cell(3,4,50);rpc_deg=cell(3,4,50);
ds_out=cell(3,4,50);ds_in=cell(3,4,50);
dsc_out=cell(3,4,50);dsc_in=cell(3,4,50);ds_deg=cell(3,4,50);dsc_deg=cell(3,4,50);
dp_out=cell(3,4,50);dp_in=cell(3,4,50);
dpc_out=cell(3,4,50);dpc_in=cell(3,4,50);dp_deg=cell(3,4,50);dpc_deg=cell(3,4,50);


load ('adjacencyMAT_ray_serial.mat')
load ('adjacencyMAT_ray_parallel.mat')
load ('adjacencyMAT_discrete_serial.mat')
load ('adjacencyMAT_discrete_parallel.mat')

for j=1:4
    for k=1:3
        for i=1:50
            rs_in{k,j,i}=sum(mat_rs{k,j}{i});
            rs_out{k,j,i}=sum(mat_rs{k,j}{i},2)';
            rs_deg{k,j,i}=rs_in{k,j,i}+rs_out{k,j,i};
            rsc_in{k,j,i}=sum(mat_rsc{k,j}{i});
            rsc_out{k,j,i}=sum(mat_rsc{k,j}{i},2)';
            rsc_deg{k,j,i}=rsc_in{k,j,i}+rsc_out{k,j,i};
             
            rp_in{k,j,i}=sum(mat_rp{k,j}{i});
            rp_out{k,j,i}=sum(mat_rp{k,j}{i},2)';
             rp_deg{k,j,i}=rp_in{k,j,i}+rp_out{k,j,i};
            rpc_in{k,j,i}=sum(mat_rpc{k,j}{i});
            rpc_out{k,j,i}=sum(mat_rpc{k,j}{i},2)';
             rpc_deg{k,j,i}=rpc_in{k,j,i}+rpc_out{k,j,i};
             
            ds_in{k,j,i}=sum(mat_ds{k,j}{i});
            ds_out{k,j,i}=sum(mat_ds{k,j}{i},2)';
             ds_deg{k,j,i}=ds_in{k,j,i}+ds_out{k,j,i};
            dsc_in{k,j,i}=sum(mat_dsc{k,j}{i});
            dsc_out{k,j,i}=sum(mat_dsc{k,j}{i},2)';
             dsc_deg{k,j,i}=dsc_in{k,j,i}+dsc_out{k,j,i};
            
            dp_in{k,j,i}=sum(mat_dp{k,j}{i});
            dp_out{k,j,i}=sum(mat_dp{k,j}{i},2)';
             dp_deg{k,j,i}=dp_in{k,j,i}+dp_out{k,j,i};
            dpc_in{k,j,i}=sum(mat_dpc{k,j}{i});
            dpc_out{k,j,i}=sum(mat_dpc{k,j}{i},2)';
             dpc_deg{k,j,i}=dpc_in{k,j,i}+dpc_out{k,j,i};
            
        end
    end
end

save  deg_info.mat rs_in rs_out rsc_in rsc_out rs_deg rsc_deg ...
    ds_in ds_out dsc_in dsc_out ds_deg dsc_deg  ...
    rp_in rp_out rpc_in rpc_out rp_deg rpc_deg  ...
    dp_in dp_out dpc_in dpc_out dp_deg dpc_deg 
%%
for i=1:3
    for j=1:4
        for k=1:50
            den_ds(i,j,k) = density_dir(mat_ds{i,j}{k});den_dsc(i,j,k) = density_dir(mat_dsc{i,j}{k});
             den_dp(i,j,k) = density_dir(mat_dp{i,j}{k});den_dpc(i,j,k) = density_dir(mat_dpc{i,j}{k});
              den_rs(i,j,k) = density_dir(mat_rs{i,j}{k});den_rsc(i,j,k) = density_dir(mat_rsc{i,j}{k});
               den_rp(i,j,k) = density_dir(mat_rp{i,j}{k});den_rpc(i,j,k) = density_dir(mat_rpc{i,j}{k});
        end
    end
end
mean_den_ds=mean(den_ds,3);mean_den_dsc=mean(den_dsc,3);
mean_den_dp=mean(den_dp,3);mean_den_dpc=mean(den_dpc,3);
mean_den_rs=mean(den_rs,3);mean_den_rsc=mean(den_rsc,3);
mean_den_rp=mean(den_rp,3);mean_den_rpc=mean(den_rpc,3);


save edge_density.mat den_ds den_dsc den_dp den_dpc den_rs den_rsc den_rp den_rpc
%%

gp='dsdprsrp';
gpc='dscdpcrscrpc';
cc_ds=cell(3,4,50);cc_dp=cell(3,4,50);cc_rs=cell(3,4,50);cc_rp=cell(3,4,50);
cc_dsc=cell(3,4,50);cc_dpc=cell(3,4,50);cc_rsc=cell(3,4,50);cc_rpc=cell(3,4,50);
for i=1:3
    for j=1:4
        for k=1:50
            for x=1:4
                
                temp=clustering_coef_bd(eval(sprintf('mat_%s{%d,%d}{%d}',gp(2*x-1:2*x),i,j,k)));temp=full(temp);
               
                tempc=clustering_coef_bd(eval(sprintf('mat_%s{%d,%d}{%d}',gpc(3*x-2:3*x),i,j,k)));tempc=full(tempc);
                if x==1
                    cc_ds{i,j,k}=temp;cc_dsc{i,j,k}=tempc;
                
                elseif x==2
                      cc_dp{i,j,k}=temp;cc_dpc{i,j,k}=tempc;
                elseif x==3
                      
                          cc_rs{i,j,k}=temp;cc_rsc{i,j,k}=tempc;
                      
                elseif x==4
                          cc_rp{i,j,k}=temp;cc_rpc{i,j,k}=tempc;
                end
            end
        end
    end
end
save cc.mat cc_ds cc_dsc cc_dp cc_dpc cc_rs cc_rsc cc_rp cc_rpc

%%


%%
% compare ds vs dp
for i=1:3
    for j=1:4
        [pdsp(i,j) hdsp(i,j)]=ranksum(squeeze(squeeze(den_ds(i,j,:))),squeeze(squeeze(den_dp(i,j,:))));
        [prds(i,j) hrds(i,j)]=ranksum(squeeze(squeeze(den_ds(i,j,:))),squeeze(squeeze(den_rs(i,j,:))),'tail','left');
        [prsp(i,j) hrsp(i,j)]=ranksum(squeeze(squeeze(den_rs(i,j,:))),squeeze(squeeze(den_rp(i,j,:))));
        [prdp(i,j) hrdp(i,j)]=ranksum(squeeze(squeeze(den_dp(i,j,:))),squeeze(squeeze(den_rp(i,j,:))),'tail','left');
        
    end
end

%%
%%
for i=1:50
    figure;plot(dsc_out{1,1,i},'r*')
   
end
figure;plot(ds_in{1,1,1},'r*')
figure;plot(dp_in{1,1,1},'bo')
figure;plot(rs_in{1,1,1},'k*')
figure;plot(rp_in{1,1,1},'b*')
close all
figure;plot(ds_out{1,1,1},'r*')
figure;plot(dp_out{1,1,1},'bo')
figure;plot(rs_out{1,1,1},'k*')
figure;plot(rp_out{1,1,1},'b*')

close all

figure;plot(dsc_in{1,1,1},'r*')
figure;plot(dpc_in{1,1,1},'bo')
figure;plot(rsc_in{1,1,1},'k*')
figure;plot(rpc_in{1,1,1},'b*')
%%
clear all
load deg_info.mat
close all
ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];
num=1000:200:2000;m_zeroout=cell(4,1);sd_zeroout=cell(4,1);
%% for now
imgWidth = .42; imgHeight = .25;imgHMargin = .05; imgVMargin = .05;
for k=3;%4 
    figure;
    for i=6;%1:6
        subplot('Position', [0 + imgHMargin + mod(i - 1 , 2) * (imgWidth + imgHMargin) ... 
                       1 - imgVMargin - (floor((i - 1) / 2) + 1) * (imgHeight + imgVMargin ) ... 
                       imgWidth imgHeight]);
        for j=1:50
            plot(ds_out{j,i,k},'r*')
            set(gcf,'Color','w')
            ylabel({' Out-degree '},'fontsize',14)
            xlabel({' Indices of Neurons '},'fontsize',14)
            title([ num2str(num(i)),'   neurons with cell size ',num2str(ss(k)),'  '],'fontsize',14);
            %ylim([0 20])
            hold on;
        zero_out_dsc(j,i )=length(find(dsc_out{j,i,k}==0))/length(dsc_out{j,i,k});
        end
        
    end
    m_zeroout_dsc=mean(zero_out_dsc);
    sd_zeroout_dsc=std(zero_out_dsc);
    figure;
    for i=1:6
        subplot('Position', [0 + imgHMargin + mod(i - 1, 2) * (imgWidth + imgHMargin) ... 
                       1 - imgVMargin - (floor((i - 1) / 2) + 1) * (imgHeight + imgVMargin) ... 
                       imgWidth imgHeight]);
        for j=1:50
             plot(dpc_out{j,i,k},'bo')
            set(gcf,'Color','w')
        
            ylabel({' Out-degree '},'fontsize',14)
            xlabel({' Indices of Neurons '},'fontsize',14)
            title([ num2str(num(i)),'   neurons with cell size ',num2str(ss(k)),'  '],'fontsize',14)
            %ylim([0 10])
            hold on
            zero_out_dpc(j,i )=length(find(dpc_out{j,i,k}==0))/length(dpc_out{j,i,k});
        end
    end
    m_zeroout_dpc=mean(zero_out_dpc);
    sd_zeroout_dpc=std(zero_out_dpc);
    figure;
    for i=1:6
        subplot('Position', [0 + imgHMargin + mod(i - 1, 2) * (imgWidth + imgHMargin) ... 
                       1 - imgVMargin - (floor((i - 1) / 2) + 1) * (imgHeight + imgVMargin) ... 
                       imgWidth imgHeight]);
        for j=1:50
             plot(rsc_out{j,i,k},'k*')
            set(gcf,'Color','w')
        
            ylabel({' Out-degree '},'fontsize',14)
            xlabel({' Indices of Neurons '},'fontsize',14)
            title([ num2str(num(i)),'   neurons with cell size ',num2str(ss(k)),'  '],'fontsize',14)
            %ylim([0 25])
            hold on
            zero_out_rsc(j,i )=length(find(rsc_out{j,i,k}==0))/length(rsc_out{j,i,k});
        end
    end
    m_zeroout_rsc=mean(zero_out_rsc);
    sd_zeroout_rsc=std(zero_out_rsc);

    figure;
    for i=1:6
        subplot('Position', [0 + imgHMargin + mod(i - 1, 2) * (imgWidth + imgHMargin) ... 
                       1 - imgVMargin - (floor((i - 1) / 2) + 1) * (imgHeight + imgVMargin) ... 
                       imgWidth imgHeight]);
        for j=1:50
            plot(rpc_out{j,i,k},'b*')
            set(gcf,'Color','w')
        
            ylabel({' Out-degree '},'fontsize',14)
            xlabel({' Indices of Neurons '},'fontsize',14)
            title([ num2str(num(i)),'   neurons with cell size ',num2str(ss(k)),'  '],'fontsize',14)
            %ylim([0 8])
            hold on
            zero_out_rpc(j,i )=length(find(rpc_out{j,i,k}==0))/length(rpc_out{j,i,k});
        end
    end
    m_zeroout_rpc=mean(zero_out_rpc);
    sd_zeroout_rpc=std(zero_out_rpc);

    m_zeroout{k,1}=[m_zeroout_dsc;m_zeroout_dpc;m_zeroout_rsc;m_zeroout_rpc];
    sd_zeroout{k,1}=[sd_zeroout_dsc;sd_zeroout_dpc;sd_zeroout_rsc;sd_zeroout_rpc];

    figure; barweb(m_zeroout{k,1},sd_zeroout{k,1})
    set(gcf,'Color','w')
    ylabel({' % of zero-Outdegree '},'fontsize',14)
    ylim([0 1.0])
    set(gca,'XTickLabel',['DS';'DP';'RS';'RP'])
    title([ ' % of neurons with Zero-Outdegree      cell size ',num2str(ss(k)),'  '],'fontsize',14) 

%xlabel({'Indices of Neurons'},'fontsize',14)
%title([ num2str(num(i)),'neurons'],'fontsize',14);
    figure; barweb(m_zeroout{k,1}',sd_zeroout{k,1}')
    set(gcf,'Color','w')
    ylabel({' % of zero-Outdegree '},'fontsize',14)
    set(gca,'XTickLabel',['1000';'1200';'1400';'1600';'1800';'2000'])
    ylim([0 1.0])
    title([ ' % of neurons with Zero-Outdegree      cell size ',num2str(ss(k)),'  '],'fontsize',14) 
% xlabel({'Indices of Neurons'},'fontsize',14)
% title([ num2str(num(i)),'neurons'],'fontsize',14);
end
%%
num=1000:200:2000;
for i=1:4
    for j=1:6
        
        cum_ds_out=reshape(cell2mat(ds_out(:,j,i)),num(j),50);
        cum_dsc_out=reshape(cell2mat(dsc_out(:,j,i)),num(j),50);
        cum_dp_out=reshape(cell2mat(dp_out(:,j,i)),num(j),50);
        cum_dpc_out=reshape(cell2mat(dpc_out(:,j,i)),num(j),50);
        cum_rs_out=reshape(cell2mat(rs_out(:,j,i)),num(j),50);
        cum_rsc_out=reshape(cell2mat(rsc_out(:,j,i)),num(j),50);
        cum_rp_out=reshape(cell2mat(rp_out(:,j,i)),num(j),50);
        cum_rpc_out=reshape(cell2mat(rpc_out(:,j,i)),num(j),50);
    end
end
dl={'cum_ds_out(:)'  'cum_dp_out(:)' 'cum_rs_out(:)' 'cum_rp_out(:)'...
    'cum_dsc_out(:)'  'cum_dpc_out(:)' 'cum_rsc_out(:)'  'cum_rpc_out(:)' };
dl=reshape(dl,2,4);
figROW = 1;
figCOL = 1;

marginWIDTH = .03;
marginHEIGHT = .05;

imgWIDTH = (1 - marginWIDTH * 2) / figCOL - marginWIDTH;
imgHEIGHT = (1 - marginHEIGHT * 2) / figROW - marginHEIGHT;

picID = figure('Color', [1 1 1]);
for rowRUN = 1:1:figROW
       for colRUN = 1:1:figCOL
%                subplot('Position', [marginWIDTH + (colRUN - 1) * (imgWIDTH + marginWIDTH) ...
%                                        1 - rowRUN * (marginHEIGHT + imgHEIGHT) ...
%                                       imgWIDTH imgHEIGHT]);
             figure; 
               hold on;
               thisCOL = ['bd' ; 'rd'];
			for drawRUN = 1 :2
				[ R  Lambda x y NODES] = cumprobdist(eval(dl{drawRUN,figCOL}));
                
                hold on
                clear x y NODES
				%loglog(x, y./  NODES, thisCOL(drawRUN) );
			end	%drawRUN = 1:1:2
% 			xlabel('k', 'FontName',' Arial',' FontSize', 7, 'FontAngle', 'Italic');
% 			ylabel('Cumulative P(k)', 'FontName', 'Arial', 'FontSize', 7, 'FontAngle', 'Italic');

               %title([ '2000 neurons  0.7350 size  '],'fontsize',14) 
               hold off;
       end     %colRUN = 1:1:figCOL
end     %rowRUN = 1:1:figROW
