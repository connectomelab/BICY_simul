clear all
num=1000:400:1800;ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9]; w_rp=cell(3,4);w_rpc=cell(3,4);
for i=1:4
    for j=1:3
         load (['ray_parallel_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'absfc','absfc_c','fcfc','fcfc_c','length_w','length_wc','ww','ww_c','tc','twc')
        m_rp_abs(j,i)=mean(absfc);
        m_rpc_abs(j,i)=mean(absfc_c);
        sd_rp_abs(j,i)=std(absfc);
        sd_rpc_abs(j,i)=std(absfc_c);
        m_rp_fc(j,i)=mean(fcfc);
        m_rpc_fc(j,i)=mean(fcfc_c);
        sd_rp_fc(j,i)=std(fcfc);
        sd_rpc_fc(j,i)=std(fcfc_c);
        m_rp_lengthw(j,i)=mean(length_w);
        m_rpc_lengthwc(j,i)=mean(length_wc);
        sd_rp_lengthw(j,i)=std(length_w);
        sd_rpc_lengthwc(j,i)=std(length_wc);
        w_rp{j,i}=ww;
        w_rpc{j,i}=ww_c;
        m_rpc_t(j,i)=mean(tc);
        m_rp_t(j,i)=mean(twc);
        sd_rpc_t(j,i)=std(tc);
        sd_rp_t(j,i)=std(twc);
    end
end
save ray_parallel_stat.mat m_rp_abs m_rpc_abs sd_rp_abs sd_rpc_abs m_rp_fc m_rpc_fc sd_rp_fc sd_rpc_fc m_rp_lengthw m_rpc_lengthwc...
    sd_rp_lengthw sd_rpc_lengthwc w_rp w_rpc m_rpc_t m_rp_t sd_rpc_t sd_rp_t
%%
clear all
num=1000:400:1800;ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9]; w_rs=cell(3,4);w_rsc=cell(3,4);
for i=1:4
    for j=1:3
         load (['ray_serial_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'absfc','absfc_c','fcfc','fcfc_c','length_w','length_wc','ww','ww_c','tc','twc')
        m_rs_abs(j,i)=mean(absfc);
        m_rsc_abs(j,i)=mean(absfc_c);
        sd_rs_abs(j,i)=std(absfc);
        sd_rsc_abs(j,i)=std(absfc_c);
        m_rs_fc(j,i)=mean(fcfc);
        m_rsc_fc(j,i)=mean(fcfc_c);
        sd_rs_fc(j,i)=std(fcfc);
        sd_rsc_fc(j,i)=std(fcfc_c);
        m_rs_lengthw(j,i)=mean(length_w);
        m_rsc_lengthwc(j,i)=mean(length_wc);
        sd_rs_lengthw(j,i)=std(length_w);
        sd_rsc_lengthwc(j,i)=std(length_wc);
        w_rs{j,i}=ww;
        w_rsc{j,i}=ww_c;
        m_rsc_t(j,i)=mean(tc);
        m_rs_t(j,i)=mean(twc);
        sd_rsc_t(j,i)=std(tc);
        sd_rs_t(j,i)=std(twc);
    end
end
save ray_serial_stat.mat m_rs_abs m_rsc_abs sd_rs_abs sd_rsc_abs m_rs_fc m_rsc_fc sd_rs_fc sd_rsc_fc m_rs_lengthw m_rsc_lengthwc...
    sd_rs_lengthw sd_rsc_lengthwc w_rs w_rsc m_rsc_t m_rs_t sd_rsc_t sd_rs_t

%%

clear all
num=1000:400:1800;ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9]; mat_rp=cell(3,4);mat_rpc=cell(3,4);
for i=1:4
    for j=1:3
         load (['ray_parallel_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'mat','mat_c')
         mat_rp{j,i}=mat;
         mat_rpc{j,i}=mat_c;
        
         clear mat mat_c
    end
end
save adjacencyMAT_ray_parallel
clear all
num=1000:400:1800;ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9]; mat_rs=cell(3,4);mat_rsc=cell(3,4);

for i=1:4
    for j=1:3
         load (['ray_serial_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'mat','mat_c')
         mat_rs{j,i}=mat;
         mat_rsc{j,i}=mat_c;
         
         clear mat mat_c
    end
end
save adjacencyMAT_ray_serial
%%
clear all
num=1000:400:1800;ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9]; mat_dp=cell(3,4);mat_dpc=cell(3,4);
for i=1:4
    for j=1:3
         load (['discrete_parallel_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'mat','mat_c')
         mat_dp{j,i}=mat;
         mat_dpc{j,i}=mat_c;
        
         clear mat mat_c
    end
end
save adjacencyMAT_discrete_parallel
clear all
num=1000:400:1800;ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9]; mat_ds=cell(3,4);mat_dsc=cell(3,4);

for i=1:4
    for j=1:3
         load (['discrete_serial_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'mat','mat_c')
         mat_ds{j,i}=mat;
         mat_dsc{j,i}=mat_c;
         
         clear mat mat_c
    end
end
save adjacencyMAT_discrete_serial

%%
clear all
num=1000:400:1800;ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9]; w_dp=cell(3,4);w_dpc=cell(3,4);
for i=1:4
    for j=1:3
         load (['discrete_parallel_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'absfc','absfc_c','fcfc','fcfc_c','length_w','length_wc','ww','ww_c','tc','twc')
        m_dp_abs(j,i)=mean(absfc);
        m_dpc_abs(j,i)=mean(absfc_c);
        sd_dp_abs(j,i)=std(absfc);
        sd_dpc_abs(j,i)=std(absfc_c);
        m_dp_fc(j,i)=mean(fcfc);
        m_dpc_fc(j,i)=mean(fcfc_c);
        sd_dp_fc(j,i)=std(fcfc);
        sd_dpc_fc(j,i)=std(fcfc_c);
        m_dp_lengthw(j,i)=mean(length_w);
        m_dpc_lengthwc(j,i)=mean(length_wc);
        sd_dp_lengthw(j,i)=std(length_w);
        sd_dpc_lengthwc(j,i)=std(length_wc);
        w_dp{j,i}=ww;
        w_dpc{j,i}=ww_c;
        m_dpc_t(j,i)=mean(tc);
        m_dp_t(j,i)=mean(twc);
        sd_dpc_t(j,i)=std(tc);
        sd_dp_t(j,i)=std(twc);
    end
end
save discrete_parallel_stat.mat m_dp_abs m_dpc_abs sd_dp_abs sd_dpc_abs m_dp_fc m_dpc_fc sd_dp_fc sd_dpc_fc m_dp_lengthw m_dpc_lengthwc...
    sd_dp_lengthw sd_dpc_lengthwc w_dp w_dpc m_dpc_t m_dp_t sd_dpc_t sd_dp_t
%%
clear all
num=1000:400:1800;ss=[0.5 2^(1/3)*0.4+0.1 4^(1/3)*0.4+0.1 0.9]; w_ds=cell(3,4);w_dsc=cell(3,4);
for i=1:4
    for j=1:3
         load (['discrete_serial_',num2str(num(j)),'_',num2str(ss(i)),'.mat'],'absfc','absfc_c','fcfc','fcfc_c','length_w','length_wc','ww','ww_c','tc','twc')
        m_ds_abs(j,i)=mean(absfc);
        m_dsc_abs(j,i)=mean(absfc_c);
        sd_ds_abs(j,i)=std(absfc);
        sd_dsc_abs(j,i)=std(absfc_c);
        m_ds_fc(j,i)=mean(fcfc);
        m_dsc_fc(j,i)=mean(fcfc_c);
        sd_ds_fc(j,i)=std(fcfc);
        sd_dsc_fc(j,i)=std(fcfc_c);
        m_ds_lengthw(j,i)=mean(length_w);
        m_dsc_lengthwc(j,i)=mean(length_wc);
        sd_ds_lengthw(j,i)=std(length_w);
        sd_dsc_lengthwc(j,i)=std(length_wc);
        w_ds{j,i}=ww;
        w_dsc{j,i}=ww_c;
        m_dsc_t(j,i)=mean(tc);
        m_ds_t(j,i)=mean(twc);
        sd_dsc_t(j,i)=std(tc);
        sd_ds_t(j,i)=std(twc);
    end
end
save discrete_serial_stat.mat m_ds_abs m_dsc_abs sd_ds_abs sd_dsc_abs m_ds_fc m_dsc_fc sd_ds_fc sd_dsc_fc m_ds_lengthw m_dsc_lengthwc...
    sd_ds_lengthw sd_dsc_lengthwc w_ds w_dsc m_dsc_t m_ds_t sd_dsc_t sd_ds_t
%%

%% corrected discrete routines
m_dp_abs_corrected=mean(length_dp_corrected_connectible);
m_dpc_abs_corrected=mean(length_dpc_corrected_connectible);
sd_dp_abs_corrected=std(length_dp_corrected_connectible);
sd_dpc_abs_corrected=std(length_dpc_corrected_connectible);
m_dp_fc_corrected=mean(corrected_fc_dp);
m_dpc_fc_corrected=mean(corrected_fc_dpc);
sd_dp_fc_corrected=std(corrected_fc_dp);
sd_dpc_fc_corrected=std(corrected_fc_dpc);
% m_dp_lengthw_corrected=mean(length_w);
% m_dpc_lengthwc_corrected=mean(length_wc);
% sd_dp_lengthw_corrected=std(length_w);
% sd_dpc_lengthwc_corrected=std(length_wc);
% m_dp_t=mean(twc);
% sd_dp_t=std(twc);
% m_dpc_t=mean(tc);
% sd_dpc_t=std(tc);
% ww_dp=ww;
% ww_dpc=ww_c;


m_ds_abs_corrected=mean(length_ds_corrected_connectible);
m_dsc_abs_corrected=mean(length_dsc_corrected_connectible);
sd_ds_abs_corrected=std(length_ds_corrected_connectible);
sd_dsc_abs_corrected=std(length_dsc_corrected_connectible);
m_ds_fc_corrected=mean(corrected_fc_ds);
m_dsc_fc_corrected=mean(corrected_fc_dsc);
sd_ds_fc_corrected=std(corrected_fc_ds);
sd_dsc_fc_corrected=std(corrected_fc_dsc);
% m_rs_lengthw_corrected=mean(length_w);
% m_rsc_lengthwc_corrected=mean(length_wc);
% sd_rs_lengthw_corrected=std(length_w);
% sd_rsc_lengthwc_corrected=std(length_wc);
% ww_rs=ww;
% ww_rsc=ww_c;
% m_rs_t=mean(twc);
% m_rsc_t=mean(tc);
% sd_rsc_t=std(tc);
% sd_rs_t=std(twc);

