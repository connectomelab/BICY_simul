% edge density calculation

num=1000:400:1800;
deno=repmat(num.*(num-1),50,1);
edgedensity=zeros(50,3,4);edgedensity_c=zeros(50,3,4);
for i=1:4
edgedensity(:,:,i)=length_w(:,:,i)./deno;
edgedensity_c(:,:,i)=length_wc(:,:,i)./deno;
end

%%
% corrected # of potential synapses 
% discrete_serial
clear all

correct_cconnectible_ds=cell(50,3,4);
correct_cconnectible_dsc=cell(50,3,4);
 ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];num=1000:400:1800;
for i=1:4
    for j=1:3
         load(['discrete_serial_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'cconnectible','cconnectible_c')
        for m=1:50
            rep1=[];rep2=[];
            temp1=cell2mat(cconnectible{m,1});
            temp2=cell2mat(cconnectible_c{m,1});
            for k=1:length(temp1(:,1))-1
                
                if temp1(k,1:2)==temp1(k+1,1:2)
                    rep1=[rep1;k+1];
                end
            end
                temp1(rep1,:)=[];
                correct_cconnectible_ds{m,j,i}=temp1;
                
            for n=1:length(temp2(:,1))-1   
                if temp2(n,1:2)==temp2(n+1,1:2)
                    rep2=[rep2;n+1];
                end
              
            end
                temp2(rep2,:)=[];
                correct_cconnectible_dsc{m,j,i}=temp2;
                
        end
    end
end
clear cconnectible cconnectible_c

for i=1:4
    for j=1:3
        for m=1:50
            length_ds_corrected_connectible(m,j,i)=length(correct_cconnectible_ds{m,j,i});
            length_dsc_corrected_connectible(m,j,i)=length(correct_cconnectible_dsc{m,j,i});
        end
    end
end


for i=1:4
    temp_w=[];temp_wc=[];
    for j=1:3
        load(['discrete_serial_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'length_w','length_wc')
        temp_w=[temp_w length_w];
        temp_wc=[temp_wc length_wc];
    end
    tmp_w(:,:,i)=temp_w;
    tmp_wc(:,:,i)=temp_wc;
end
corrected_fc_ds=tmp_w./length_ds_corrected_connectible;
corrected_fc_dsc=tmp_wc./length_dsc_corrected_connectible;

m_ds_abs_corrected=squeeze(mean(length_ds_corrected_connectible))';
m_dsc_abs_corrected=squeeze(mean(length_dsc_corrected_connectible))';
sd_ds_abs_corrected=squeeze(std(length_ds_corrected_connectible))';
sd_dsc_abs_corrected=squeeze(std(length_dsc_corrected_connectible))';
m_ds_fc_corrected=squeeze(mean(corrected_fc_ds))';
m_dsc_fc_corrected=squeeze(mean(corrected_fc_dsc))';
sd_ds_fc_corrected=squeeze(std(corrected_fc_ds))';
sd_dsc_fc_corrected=squeeze(std(corrected_fc_dsc))';





save corrected_DS.mat correct_cconnectible_ds correct_cconnectible_dsc corrected_fc_ds corrected_fc_dsc m_ds_abs_corrected m_dsc_abs_corrected sd_ds_abs_corrected...
    sd_dsc_abs_corrected m_ds_fc_corrected m_dsc_fc_corrected sd_ds_fc_corrected sd_dsc_fc_corrected length_ds_corrected_connectible length_dsc_corrected_connectible


%%
% corrected # of potential synapses 
% discrete_parallel
clear all
 ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9];num=1000:400:1800;
correct_cconnectible_dp=cell(50,3,4);
correct_cconnectible_dpc=cell(50,3,4);
for i=1:4
    for j=1:3
        load(['discrete_parallel_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'cconnectible','cconnectible_c')
        for m=1:50
            rep1=[];rep2=[];
            temp1=sortrows(cell2mat(cconnectible{m,1}),2);
            temp2=sortrows(cell2mat(cconnectible_c{m,1}),2);
            for k=1:length(temp1(:,1))-1
                
                if temp1(k,1:2)==temp1(k+1,1:2)
                    rep1=[rep1;k+1];
                end
            end
                temp1(rep1,:)=[];
                correct_cconnectible_dp{m,j,i}=temp1;
                
            for n=1:length(temp2(:,1))-1   
                if temp2(n,1:2)==temp2(n+1,1:2)
                    rep2=[rep2;n+1];
                end
              
            end
                temp2(rep2,:)=[];
                correct_cconnectible_dpc{m,j,i}=temp2;
                
        end
    end
end

for i=1:4
    for j=1:3
        for m=1:50
            length_dp_corrected_connectible(m,j,i)=length(correct_cconnectible_dp{m,j,i});
            length_dpc_corrected_connectible(m,j,i)=length(correct_cconnectible_dpc{m,j,i});
%             
%              length_dp_corrected_connectible(m,j)=length(correct_cconnectible_dp{m,j});
%             length_dpc_corrected_connectible(m,j)=length(correct_cconnectible_dpc{m,j});
            
        end
    end
end

for i=1:4
    temp_w=[];temp_wc=[];
    for j=1:3
        load(['discrete_parallel_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'length_w','length_wc')
        temp_w=[temp_w length_w];
        temp_wc=[temp_wc length_wc];
    end
     tmp_w(:,:,i)=temp_w;
    tmp_wc(:,:,i)=temp_wc;
end
corrected_fc_dp=tmp_w./length_dp_corrected_connectible;
corrected_fc_dpc=tmp_wc./length_dpc_corrected_connectible;

m_dp_abs_corrected=squeeze(mean(length_dp_corrected_connectible))';
m_dpc_abs_corrected=squeeze(mean(length_dpc_corrected_connectible))';
sd_dp_abs_corrected=squeeze(std(length_dp_corrected_connectible))';
sd_dpc_abs_corrected=squeeze(std(length_dpc_corrected_connectible))';
m_dp_fc_corrected=squeeze(mean(corrected_fc_dp))';
m_dpc_fc_corrected=squeeze(mean(corrected_fc_dpc))';
sd_dp_fc_corrected=squeeze(std(corrected_fc_dp))';
sd_dpc_fc_corrected=squeeze(std(corrected_fc_dpc))';



save corrected_DP.mat correct_cconnectible_dp correct_cconnectible_dpc corrected_fc_dp corrected_fc_dpc m_dp_abs_corrected m_dpc_abs_corrected ... 
    sd_dp_abs_corrected sd_dpc_abs_corrected m_dp_fc_corrected m_dpc_fc_corrected sd_dp_fc_corrected sd_dpc_fc_corrected length_dp_corrected_connectible length_dpc_corrected_connectible

%%
rep=[];
for i=1:length(temp (:,1))-1
    if temp (i,1:2)==temp (i+1,1:2)
       rep=[rep;i+1];
    end
end
temp (rep,:)=[];

act2dpc=sortrows(act2dpc,2);
rep=[];
for i=1:length(act2dpc(:,1))-1
    if act2dpc(i,1:2)==act2dpc(i+1,1:2)
       rep=[rep;i+1];
    end
end

act2dpc(rep,:)=[];

rep=[];
for i=1:length(act2dp(:,1))-1
    if act2dp (i,1:2)==act2dp(i+1,1:2)
       rep=[rep;i+1];
    end
end

act2dp(rep,:)=[];
act2dsc=sortrows(act2dsc,2);
ind=nonzeros(act2dsc);

act2dsc(1:657,:)=[];
