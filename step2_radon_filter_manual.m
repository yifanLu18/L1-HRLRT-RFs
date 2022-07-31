
% Perform filtering, inverse transformation, drawing, etc. on the result of forward transformation

clear;clc;close all;
set(0,'defaultfigurecolor','w');
km2deg = 6371*2*pi/360;

load('pre_processing.mat');

%% Calculate the theoretical position of the filter

% Ridge topography, station, source latitude and longitude information
lo_topo1 = -20/km2deg;
lo_topo2 = 20/km2deg;
la_topo1 = 0/km2deg;
la_topo2 = 0/km2deg;

la_sta1 = -50/km2deg;
lo_sta1 = 0;
la_sta2 = 50/km2deg;
lo_sta2 = 0;
[evla,evlo] = cal_event_location((la_sta1+la_sta2)/2,(lo_sta1+lo_sta2)/2,deg,baz,0.05,1);

depth_src = 0;
V_rayleigh = 3.36 * 0.9194;

[k_12 k_21] = cal_topo_visiualRF(la_topo1,lo_topo1,la_topo2,lo_topo2,la_sta1,lo_sta1,la_sta2,lo_sta2,evla,evlo,depth_src,V_rayleigh);

%% Figure 1: Input Data Plot

figure(1);
set(gcf,'Units','centimeter','Position',[5 4 30 6]);

subplot(1,2,1);
imagesc(xoff,t,shot,[-.5,.5]);
title('RF to be transformed     time lag = 2s');
xlabel('Distance (km)');ylabel('Time (s)');
set(gca,'ytick',0:2:18);ylim([-2,10]);
caxis([-0.5,0.5]);

subplot(1,2,2);
for i = 1:size(shot,2)
    plot(shot(:,i)/max(shot(:,i))*3 + xoff(i),t,'k');
    hold on;
end
title('RF to be transformed     time lag = 2s');
xlabel('Distance (km)');ylabel('Time (s)');
xlim([min(xoff),max(xoff)]);ylim([-2,10]);
set(gca,'ytick',min(t):2:max(t));
set(gca,'ydir','rev');

%% Figure 2: High Resolution Radon Transform

figure(2);
set(gcf,'Units','centimeter','Position',[5 6 30 6]);
subplot(1,2,2);
imagesc(p,t,R(:,1:length(t))',[-0.005,0.005]);colorbar;
title('high resolution LRT');xlabel('Slowness (s/km)');ylabel('Delay time (s)');
ylim([min(t),max(t)]);set(gca,'ytick',min(t):2:max(t));
hold on;
plot(-1/k_21,0,'rx','MarkerSize',10,'linewidth',3);

%% Figure 3: Filtering Operation

figure(3);
set(gcf,'Units','centimeter','Position',[5 8 30 6]);
subplot(1,2,1);
imagesc(R',[-0.005,0.005]);
ylim([1,length(t)]);colorbar;
xlabel('Slowness points');ylabel('tau points');
title('radon panel & select the filter window');

% ===============================================================   Filter window range  ========================================================
BW = zeros(size(R'));
windows_all_x = [];
windows_all_y = [];
windows_num = input('Please input the number of windows\n');
for j = 1:windows_num
    windows_number  = input('Please input the numbers of filter£º\n');
    [windows_x,windows_y] = ginput(windows_number);
    BW_temp = roipoly(R',windows_x,windows_y);
    
    for i = 1:windows_number-1
        hold on;
        plot([windows_x(i),windows_x(i+1)],[windows_y(i),windows_y(i+1)],'r','linewidth',2);
    end
    plot([windows_x(1),windows_x(end)],[windows_y(1),windows_y(end)],'r','linewidth',2);
    
    windows_all_x(j,:) = windows_x;
    windows_all_y(j,:) = windows_y;
    
    BW = BW+BW_temp;
end
% ===============================================================   Filter window range  ========================================================

R_filter = (BW') .* R;

subplot(1,2,2);
imagesc(R_filter',[-0.01,0.01]);colorbar;
ylim([1,length(t)]);
pause(0.01);

%% Inverse transformation

% Inverse transformation of the original data to check the reversibility of the transformation
ref_dist = 0;
line_model = 'linear';
M_org = Radon_inverse(t,p,R,xoff,ref_dist,line_model);
M_org_image = M_org';

% Inverse transform of filtered data
M_filter = Radon_inverse(t,p,R_filter,xoff,ref_dist,line_model);
M_filter_image = M_filter';
M_remain = M_org - M_filter;
M_remain_image = M_remain';

% Plot (contrast before and after filtering)
figure(4);
set(gcf,'Units','centimeter','Position',[5 10 30 6]);
subplot(1,2,1);
imagesc(xoff,t,shot,[-0.5,0.5]);
xlabel('Distance / km');ylabel('Time / s');title('Org RF');colorbar;ylim([-2,10]);
subplot(1,2,2);
imagesc(xoff,t,M_remain_image,[-0.5,0.5]);
xlabel('Distance / km');ylabel('Time / s');title('Filtered RF');colorbar;caxis([-0.5,0.5]);ylim([-2,10]);

% Drawing (showing the filtered part)
figure(5);
set(gcf,'Units','centimeter','Position',[5 12 30 6]);
subplot(1,2,1);
imagesc(xoff,t,shot,[-0.5,0.5]);
xlabel('Distance / km');ylabel('Time / s');title('Org RF');colorbar;ylim([-2,10]);
subplot(1,2,2);
imagesc(xoff,t,M_filter_image,[-0.5,0.5]);
xlabel('Distance / km');ylabel('Time / s');title('Filtered RF');colorbar;ylim([-2,10]);

% Drawing (compare the results of the forward and reverse transformation of the original data, test the selection of the 
% number of iterations and the discrete step size of p)
figure(6);
set(gcf,'Units','centimeter','Position',[5 14 30 6]);
for i = 1:1:size(shot,2)
    f1 = plot(xoff(i) + shot(:,i),t,'r');
    hold on;
    f2 = plot(xoff(i) + M_org_image(:,i),t,'b');
end
set(gca,'ydir','rev');
xlim([min(xoff),max(xoff)]);ylim([-2,10]);
title_name = strcat('iter =',32,num2str(maxiter),32,32,32,32,32,'dp =',32,num2str(dp));
title(title_name);
ylabel('Time / s');xlabel('Distance / km')
legend([f1 f2],'Org','iRadon');

% Plot (compare raw data with filtered results)
figure(7);
set(gcf,'Units','centimeter','Position',[5 16 30 6]);
for i = 1:1:size(shot,2)
    f1 = plot(xoff(i) + shot(:,i),t,'r');
    hold on;
    f2 = plot(xoff(i) + M_remain_image(:,i),t,'b');
end
set(gca,'ydir','rev');
xlim([min(xoff),max(xoff)]);ylim([-2,10]);
title_name = strcat('iter =',32,num2str(maxiter),32,32,32,32,32,'dp =',32,num2str(dp));
title(title_name);
ylabel('Time / s');xlabel('Distance / km')
legend([f1 f2],'Org','iRadon');

%% save the dataset

M_filter_image = M_filter_image(pad_npoint+1:end-pad_npoint,:);
M_org_image = M_org_image(pad_npoint+1:end-pad_npoint,:);
M_remain_image = M_remain_image(pad_npoint+1:end-pad_npoint,:);
RFs = shot(pad_npoint+1:end-pad_npoint,:);
t_addpad = t;
t = t(pad_npoint+1:end-pad_npoint);

save('process.mat','RFs','M_org_image','M_filter_image','M_remain_image','xoff','t');

R = R(:,1:length(t))';
R_filter = R_filter(:,1:length(t))';
save('radon.mat','R','R_filter','windows_all_x','windows_all_y','p','t');

save('backup.mat');



