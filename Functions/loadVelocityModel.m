%% Load SAWum-NA2  (Yuan et al., 2011 - GJI)
% download from
% seismo.berkeley.edu/wiki_br/Regional.vectorial.tomography.download
clear all; close all; clc
cd ~/Scattered_Waves/Data/Velocity_Models/Mantle/Yuan2011GJI/

load F.interpolant.NA.mat
% Fs:   Vs
% Faz:  anisotropy direction
% Fx:   radial anisotropy
% Fg:   azimuthal anisotropy strength

lats=10:0.5:80; lons=(181:0.5:320)-360; deps=40:10:450;
[xq,yq]=meshgrid(lons+360,lats);
Vs_Model.Vs=zeros(length(lats),length(lons),length(deps));
for j=1:length(deps)
    Vs_Model.Vs(:,:,j)=Fs(yq,xq,deps(j)*ones(size(xq)))./1000;
    Vs_Model.AnisotropyDir(:,:,j)=Faz(yq,xq,deps(j)*ones(size(xq)));
    Vs_Model.RadAnisotropy(:,:,j)=Fx(yq,xq,deps(j)*ones(size(xq)));
    Vs_Model.AzAnisotropyStrength(:,:,j)=Fg(yq,xq,deps(j)*ones(size(xq)));
end

Vs_Model.Name='Yuan_etal_2011_GJI';
Vs_Model.Latitude=lats;
Vs_Model.Longitude=lons;
Vs_Model.Depth=deps;

save ../Vs/Yuan_etal_2011_GJI.mat Vs_Model


%% plot map sections
dep=60; type='Vs';

eval(['V=Vs_Model.' type ';']);

States = shaperead('usastatelo', 'UseGeoCoords', true); coast=load('coast');
[~,ind]=min(abs(deps-dep));

%lims=[-118 -94 36 52];
lims=[-87 -78 26 37]; 
figure; hold on; title([num2str(deps(ind)) 'km:  ' type])
contourf(lons,lats,V(:,:,ind),20); shading flat
for j = 1:length(States); plot(States(j).Lon, States(j).Lat,'k-'); end
plot(coast.long,coast.lat,'k-')
axis(lims); daspect([1/(111.16*distance(mean(lims(3:4)),0,mean(lims(3:4)),1)), 1/111.16, 1]);
c=colorbar; colormap(flipud(jet)); ylabel(c,type);


%% Load GYPSUM
clear all; close all; clc
cd ~/Scattered_Waves/Data/Velocity_Models/Mantle/Vp/GYPSUM/

info = ncinfo('GYPSUMP_kmps.nc');
lats=ncread('GYPSUMP_kmps.nc', 'latitude');
lons=ncread('GYPSUMP_kmps.nc', 'longitude');
deps = ncread('GYPSUMP_kmps.nc', 'depth');
vp = ncread('GYPSUMP_kmps.nc','vp');
Vp_Model.Name='GYPSUM';
Vp_Model.Latitude=lats;
Vp_Model.Longitude=lons;
Vp_Model.Depth=deps;

vpt=nan(length(lats),length(lons),length(deps));
for i =1:length(deps)
    vpt(:,:,i)=vp(:,:,i)';
end

[~,i]=min(abs(deps-100));
if min(lons)<0; ax=[-180 180 -90 90]; else ax=[0 360 -90 90]; end
figure; hold on; axis(ax); daspect([1 0.6 1]); 
contourf(lons, lats, vpt(:,:,i),'linestyle','none')
coast=load('coast'); colorbar; colormap(flipud(jet));
plot(coast.long+360,coast.lat,'k'); plot(coast.long, coast.lat,'k');

Vp_Model.Vp = vpt;

save ../GYPSUM.mat Vp_Model

%% Load TX 2011 (Grand, 2011) 
clear all; close all; clc

name='TX2011'; %'S362WMANIM';
str='_percent'; %'_kmps';
cd(['~/Scattered_Waves/Data/Velocity_Models/Mantle/Vs/' name]);

info = ncinfo([name str '.nc']);
lats=ncread([name str '.nc'], 'latitude');
lons=ncread([name str '.nc'], 'longitude');
deps = ncread([name str '.nc'], 'depth');
%vs = ncread([name '.nc'],'vs');

dvs = ncread([name str '.nc'],'dvs');
load rem.mat
vs=zeros(size(dvs));
for id=1:length(deps)
    [~,irem]=min(abs(deps(id)-rem(:,1)));
    vs(:,:,id)=(1+dvs(:,:,id)./100).*rem(irem,2);
end

Vs_Model.Name=name;
Vs_Model.Latitude=lats;
Vs_Model.Longitude=lons;
Vs_Model.Depth=deps;

vst=nan(length(lats),length(lons),length(deps));
for i =1:length(deps)
    vst(:,:,i)=vs(:,:,i)';
end

[~,i]=min(abs(deps-100));
if min(lons)<0; ax=[-180 180 -90 90]; else ax=[0 360 -90 90]; end
figure; hold on; axis(ax); daspect([1 0.6 1]); 
contourf(lons, lats, vst(:,:,i),'linestyle','none')
coast=load('coast'); colorbar; colormap(flipud(jet));
plot(coast.long+360,coast.lat,'k'); plot(coast.long, coast.lat,'k');

Vs_Model.Vs = vst;

eval(['save ../' name '.mat Vs_Model']);


%%
% load SEMUM velocity model
% downloaded as netCDF file from http://www.iris.edu/dms/products/emc-semum/
% chose version expressed as absolute velocities
clear all; close all; clc
cd ~/Scattered_Waves/Data/Velocity_Models/Mantle/Vs/SEMUM/

info = ncinfo('SEMUM_kmps.nc');

lats=ncread('SEMUM_kmps.nc', 'latitude');
lons=ncread('SEMUM_kmps.nc', 'longitude');
deps = ncread('SEMUM_kmps.nc', 'depth');
vs = ncread('SEMUM_kmps.nc','vs');

Vs_Model.Name='SEMUM';
Vs_Model.Latitude=lats;
Vs_Model.Longitude=lons;
Vs_Model.Depth=deps;

vst=nan(length(lats),length(lons),length(deps));
for i =1:length(deps)
    vst(:,:,i)=vs(:,:,i)';
end

if min(lons)<0; ax=[-180 180 -90 90]; else ax=[0 360 -90 90]; end

figure; hold on; axis(ax); daspect([1 0.6 1]); daspect([1 0.6 1]); 
contourf(lons, lats, vst(:,:,2),'linestyle','none')
coast=load('coast'); colorbar; colormap(flipud(jet));
plot(coast.long+360,coast.lat,'k'); plot(coast.long, coast.lat,'k');

Vs_Model.Vs = vst;

save ../SEMUM.mat Vs_Model

%%
clear all; %figure;
Data='~/Scattered_Waves';
load([Data '/Data/Velocity_Models/Mantle/Vs/Wagner10.mat'])
[ob_x, ob_y, ob_z] = meshgrid(Vs_Model.Longitude,Vs_Model.Latitude,Vs_Model.Depth);
hold on
[~,ilat]=min(abs(Vs_Model.Latitude-43)); [~,ilon]=min(abs(Vs_Model.Longitude+108.2));
plot(squeeze(Vs_Model.Vs(ilat,ilon,:)), Vs_Model.Depth,'b'); axis ij
xlabel('Vs (km/s)'); ylabel('Depth (km)'); ylim([0 400])
