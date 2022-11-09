% cruise with topography
clear all;close all;
files=dir('*.mat');
hfig1=figure(1)
colormap(m_colmap('blues'));  
m_proj('miller','lon',[-180 180],'lats',[-90 90])
% m_elev('contour',[-3500:500:-500],'edgecolor','b');
[CS,CH]=m_etopo2('contourf',[-7000:1000:-1000 -500 -200 0 ],'edgecolor','none');
clabel(CS,CH)
% colorbar
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none'); 
% m_grid
% m_grid('linest','none','tickdir','out','box','fancy','fontsize',16);
colors=[[0 0 0];[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.6350 0.0780 0.1840];[0.4940 0.1840 0.5560];[0.6350 0.0780 0.1840]];

Topo_all=[];Distance_all=[];KCO2_all=[];KCO2_all_to_distance=[];U10_all=[];Hs_all=[];
for ifile=1:length(files)
    load([files(ifile).name])
    distance=sw_dist(lon,lat,'km');
    
    Topo_all=[Topo_all;reshape(Topo,length(Topo),1)];
    Distance_all=[Distance_all;reshape(distance,length(distance),1)];
    KCO2_all=[KCO2_all;reshape(kCO2660_cm_hr,length(kCO2660_cm_hr),1)];
    KCO2_all_to_distance=[KCO2_all_to_distance;reshape(kCO2660_cm_hr(2:end),length(kCO2660_cm_hr(2:end)),1)];
    U10_all=[U10_all;u10n];
    Hs_all=[Hs_all;Hs_ECMWF_tot_m];
hold on;
% m_scatter(lon(2:end),lat(2:end),100,distance)
m_plot(lon,lat,'o','markersize',5,'markeredgecolor',colors(ifile,:))
hold on;
m_grid('linest','none','tickdir','out','box','fancy','fontsize',10);
caxis([-7000 000]);
[ax,h]=m_contfbar([.55 .75],.8,CS,CH,'endpiece','no','axfrac',.05);
title(ax,'meters')
set(gcf,'color','w');

end
% hfig1=figure(1);
set(hfig1,'paperunits','inches','paperposition',[0 0 10 8]);
exportgraphics(hfig1,['where_cruise_happen_deep_or_shallow.png'],'resolution',200,'ContentType','vector') %does not leaves white border
% 
%% statistical topography distribution at in-situ observation dataset
% The reason I did this check is based on Va in Romero 2019 is calculated
% from deep water dispersive relationship (keep in mind)
[loc,val]=find(Topo_all>=-20&~isnan(KCO2_all));
min(Distance_all)
% [loc,val]=find(Topo_all>=-20);
% 1/48 1/24 1/12 1/4 
% three bins
degree_bin=[1/48,1/24,1/12,1/4];
km_bin=degree_bin*110;
x_km=[km_bin+[0,km_bin(1:end-1)]]/2;
x_degree=[degree_bin+[0,degree_bin(1:end-1)]]/2;
[loc1,val1]=find(Distance_all<=km_bin(1)&~isnan(KCO2_all_to_distance));
[loc2,val2]=find(Distance_all>km_bin(1)&Distance_all<=km_bin(2)&~isnan(KCO2_all_to_distance));
[loc3,val3]=find(Distance_all>km_bin(2)&Distance_all<=km_bin(3)&~isnan(KCO2_all_to_distance));
[loc4,val4]=find(Distance_all>km_bin(3)&Distance_all<=km_bin(4)&~isnan(KCO2_all_to_distance));
[loc5,val5]=find(Distance_all>km_bin(4)&~isnan(KCO2_all_to_distance));

loc_tmp=find(~isnan(KCO2_all_to_distance));
Percentage_all=[length(loc1),length(loc2),length(loc3),length(loc4),length(loc5)]/length(loc_tmp);
% xticklabels={num2str(km_bin),'(km)'}
hfig2=figure(2)
t = tiledlayout(1,1);
ax1 = axes(t);
% plot(ax1,km_bin,Percentage_all(1:end-1)*100),'o','color',none)
ax1.XColor = 'k';
ax1.YColor = 'k';
hold on;
bar(ax1,[0,km_bin],Percentage_all*100);
grid minor
xlim([0 55])
xlabel('km', 'FontSize', 14)
ylabel('percentage of data samples (%)', 'FontSize', 14)
set(gca,'ylim',[0 100],'ytick',[0:10:100], 'FontSize', 14)
xtickformat('%.2f')


ax2 = axes(t);
plot(ax2,[0,degree_bin],Percentage_all*100,'o','color','none');
xlim([0 1/2])

ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
ax2.TickLabelInterpreter = 'latex';
set(gca,'ytick',[])
xlabel('^o (degree)')

% 
set(ax1,'xtick',[0,km_bin],'ytick',[0:10:100])
set(ax2,'xtick',[0,degree_bin],'xticklabel',{'0','$\frac{1}{48}$','$\frac{1}{24}$','$\frac{1}{12}$','$\frac{1}{4}$'}, 'FontSize', 14,'fontweight','bold')
set(hfig2,'paperunits','inches','paperposition',[0 0 10 8]);
exportgraphics(hfig2,['distribution_of_spatial_distribution.png'],'resolution',200,'ContentType','vector') %does not leaves white border
% 
%% check the condition when spatial distrubution is <1/48 degree
[loc1,val1]=find(Distance_all<=km_bin(1)&~isnan(KCO2_all_to_distance));
u10_tmp=U10_all(loc1);
hs_tmp=Hs_all(loc1);
wind_bin=[0:5:30];median_wind=[];ytxt=[];
pp_wind=[];
for ibin=1:length(wind_bin)-1
    median_wind(ibin)=(wind_bin(ibin)+wind_bin(ibin+1))/2;
    if ibin==1
pp_wind(ibin)=length(u10_tmp(u10_tmp<wind_bin(ibin+1)))/length(u10_tmp);
    ytxt(ibin)=length(u10_tmp(u10_tmp<wind_bin(ibin+1)));

    else
pp_wind(ibin)=length(u10_tmp(u10_tmp<wind_bin(ibin+1)&u10_tmp>=wind_bin(ibin)))/length(u10_tmp);
    ytxt(ibin)=length(u10_tmp(u10_tmp<wind_bin(ibin+1)&u10_tmp>=wind_bin(ibin)));

    end
        

end
%%
close all
hfig3=figure
subplot(1,2,1)
plot(u10_tmp,hs_tmp,'o','color',[0.5,0.5,0.5])
xlabel('u10(m/s)')
ylabel('Hs(m)')
grid on
set(gca,'FontSize', 10)


subplot(1,2,2)
bar(median_wind,pp_wind*100)
xlim([0 30])
ylim([0 100])
hold on;
for i=1:length(median_wind)  % iterate over number of bar objects
  text(median_wind(i)-2,pp_wind(i)*100+5, num2str(ytxt(i)));
  hold on;
end
grid on
xlabel('range of U10 (m/s)')
ylabel('percentage of data sample (%)')
set(gca,'xtick',[0:5:30],'ytick',[0:10:100],'FontSize', 10)
set(hfig3,'paperunits','inches','paperposition',[0 0 10 8]);
exportgraphics(hfig3,['distribution_of_high_resolution_wind.png'],'resolution',200,'ContentType','vector') %does not leaves white border
% 