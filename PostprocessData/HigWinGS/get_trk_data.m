clear all;close all;
% calculate Cg
file_glo='../ww3.2013.nc';
file_trk='/projects/DEIKE/XiaohuiZhou/ReevaluatedFluxData/StructureData/HigWinGS/CO2_Flux.mat'

%file_trk='../ww3.2019_trck.nc';
% calculate Cg
lat_glo_tmp=ncread(file_glo,'latitude');
lon_glo_tmp=ncread(file_glo,'longitude');
time_glo_tmp=ncread(file_glo,'time')+datenum('1990-01-01 00:00:00');

%lat_trk=ncread(file_trk,'latitude');
%lon_trk=ncread(file_trk,'longitude');
%time_trk=ncread(file_trk,'time');
load(file_trk);
lon_trk=CO2_Flux.lon;
lat_trk=CO2_Flux.lat;
time_trk=CO2_Flux.time;

xloc=find(lon_glo_tmp>=min(lon_trk)-2&lon_glo_tmp<=max(lon_trk)+2);
yloc=find(lat_glo_tmp>=min(lat_trk)-2&lat_glo_tmp<=max(lat_trk)+2);
tloc=find(time_glo_tmp>=min(time_trk)-1&time_glo_tmp<=max(time_trk)+1);
lon_glo=lon_glo_tmp(xloc);
lat_glo=lat_glo_tmp(yloc);
time_glo=time_glo_tmp(tloc);
count_x=length(xloc);
count_y=length(yloc);
count_t=length(tloc);
startLoc = [xloc(1) yloc(1) tloc(1)]; % Start location along each coordinate
count  = [count_x count_y count_t]; % Read until the end of each dimension
%% read global map 
uust_glo=ncread(file_glo,'uust',startLoc,count);
vust_glo=ncread(file_glo,'vust',startLoc,count);
hs_glo=ncread(file_glo,'hs',startLoc,count);
uwnd_glo=ncread(file_glo,'uwnd',startLoc,count);
vwnd_glo=ncread(file_glo,'vwnd',startLoc,count);
wnd_glo=sqrt(uwnd_glo.^2+vwnd_glo.^2);
wva_glo=ncread(file_glo,'wva',startLoc,count);
wcc_glo=ncread(file_glo,'wcc',startLoc,count);
hs_glo=ncread(file_glo,'hs',startLoc,count);
%%
ust=[];vst=[];hs=[];wva=[];u10=[];v10=[];wcc=[];
for ii=1:length(xloc)
    for jj=1:length(yloc)
        ust_tmp=double(squeeze(uust_glo(ii,jj,:)));
        vst_tmp=double(squeeze(vust_glo(ii,jj,:)));       
        hs_tmp=double(squeeze(hs_glo(ii,jj,:)));
        wva_tmp=double(squeeze(wva_glo(ii,jj,:)));
        wcc_tmp=double(squeeze(wcc_glo(ii,jj,:)));
        u10_tmp=double(squeeze(uwnd_glo(ii,jj,:)));
        v10_tmp=double(squeeze(vwnd_glo(ii,jj,:)));
        
        ust(ii,jj,:)=interp1(double(time_glo),ust_tmp,time_trk);
        vst(ii,jj,:)=interp1(double(time_glo),vst_tmp,time_trk);
        u10(ii,jj,:)=interp1(double(time_glo),u10_tmp,time_trk);
        v10(ii,jj,:)=interp1(double(time_glo),v10_tmp,time_trk);
        hs(ii,jj,:)=interp1(double(time_glo),hs_tmp,time_trk);
        wva(ii,jj,:)=interp1(double(time_glo),wva_tmp,time_trk);
        wcc(ii,jj,:)=interp1(double(time_glo),wcc_tmp,time_trk);
    end
end

ust_trk=[];vst_trk=[];
hs_trk=[];wva_trk=[];wcc_trk=[];uwnd_trk=[];vwnd_trk=[];
for it=1:length(time_trk)
    ust_trk(it)=interp2(lon_glo',lat_glo',ust(:,:,it)',lon_trk(it),lat_trk(it));
    vst_trk(it)=interp2(lon_glo',lat_glo',vst(:,:,it)',lon_trk(it),lat_trk(it));
    uwnd_trk(it)=interp2(lon_glo',lat_glo',u10(:,:,it)',lon_trk(it),lat_trk(it));
    vwnd_trk(it)=interp2(lon_glo',lat_glo',v10(:,:,it)',lon_trk(it),lat_trk(it));
    hs_trk(it)=interp2(lon_glo',lat_glo',hs(:,:,it)',lon_trk(it),lat_trk(it));
    wva_trk(it)=interp2(lon_glo',lat_glo',wva(:,:,it)',lon_trk(it),lat_trk(it));
    wcc_trk(it)=interp2(lon_glo',lat_glo',wcc(:,:,it)',lon_trk(it),lat_trk(it));
end
save('trk_data.mat','ust_trk','vst_trk','uwnd_trk','vwnd_trk','hs_trk','wva_trk','wcc_trk','lon_trk','lat_trk','time_trk');
