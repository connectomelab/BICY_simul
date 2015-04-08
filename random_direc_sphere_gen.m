% %% angle generation
% clear all
% RandStream.setDefaultStream ...
%      (RandStream('mt19937ar','seed',sum(100*clock)));
% 
% for k=1:6
%     num=1000:200:2000;
%     angle=cell(100,1);
%     for i=1:100
%         dir=zeros(num(k),2);
%         for j=1:num(k)
%             phi = 2*pi*rand(1); % growth direction
%             theta  = 2*pi*rand(1);
%             dir(j,:)=[phi theta];
%         end;
%         angle{i}=dir;
%     end;
%     save (['angle' num2str(num(k)) '.mat'])
%     clear all
% end

%% random direction generation on a 3D surface
% http://www.math.niu.edu/~rusin/known-math/96/sph.rand (3) The trig method
clear all
RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',sum(100*clock)));

for k=1:6
    num=1000:200:2000;
    dir=cell(100,1);
    for i=1:100
        para=zeros(num(k),3);
        for j=1:num(k)
            z=-1+2*rand(1);
            phi=2*pi*rand(1);
            r=sqrt(1-z^2);
            para(j,:)=[z phi r];
           
        end;
        dir{i}=para;
    end;
    save (['rand_dir_sphere' num2str(num(k)) '.mat'])
    clear all
end