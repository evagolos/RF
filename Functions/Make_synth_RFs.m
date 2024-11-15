% Make synthetic receiver functions for a given velocity model.
% Written by Eva Golos, 7/2021. Based off 'RFcomp' functionality from Junlin Hua/Isabella Gama's code.

basedir= '/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin/';
javaaddpath([basedir 'Functions/taup/lib/TauP-1.1.7.jar']);
javaaddpath([basedir 'Functions/taup/lib/log4j-1.2.8.jar']);
javaaddpath([basedir 'Functions/taup/lib/seisFile-1.0.1.jar']);

% Load synthetic velocity model
VsModelname= 'MITPS_10kmZ'; VpModelname= 'MITPS_10kmZ';
tmp_vfile_s = dir(['/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin/Data/Velocity_Models/Mantle/Vs/' VsModelname '.mat']);
tmp_vfile_p = dir(['/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin/Data/Velocity_Models/Mantle/Vp/' VpModelname '.mat']);
orimodel_s = load([tmp_vfile_s.folder '/' tmp_vfile_s.name]);
orimodel_p = load([tmp_vfile_p.folder '/' tmp_vfile_p.name]);

pdir= pwd;
% Event and station information
STA.Latitude= 32.31; STA.Longitude= -110.785; % Use TUC info for now.
EQ.Depth= 33; EQ.Latitude= 39.0300; EQ.Longitude= 143.4400; % Use parameters from EQ_1998_150_18_18 event for now. (SNR=4)
Phase= 'P'; % Must be P or S for TauP
REM= 'ak135';


% Get ray parameter
cd('/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin/Functions/Matlab_TauP');
%addpath('/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin/Functions/taup');
Taup_output = Matlab_TauP('Path',REM,EQ.Depth,Phase,...
            'sta',[STA.Latitude STA.Longitude],...
            'evt',[EQ.Latitude EQ.Longitude]);
rayp= Taup_output.rayParam*pi/180/(6371*pi/180);

cd(pdir);


%rp_ref= [rp_min,rp_med,rp_max];
rp_ref= rayp;
depth = 0:450;

[LonG,LatG,DepG] = meshgrid(orimodel_s.Vs_Model.Longitude,orimodel_s.Vs_Model.Latitude,orimodel_s.Vs_Model.Depth);
latmin = orimodel_s.Vs_Model.Latitude(1);
latmax = orimodel_s.Vs_Model.Latitude(end);
lonmin = orimodel_s.Vs_Model.Longitude(1);
lonmax = orimodel_s.Vs_Model.Longitude(end);
lat_model = latmin:1:latmax;
lon_model = lonmin:1:lonmax;
[LonM,LatM,DepM] = meshgrid(lon_model,lat_model,depth');
Vs_M = interp3(double(LonG),double(LatG),double(DepG),...
    double(orimodel_s.Vs_Model.Vs),double(LonM),double(LatM),double(DepM));
Vp_M = interp3(double(LonG),double(LatG),double(DepG),...
    double(orimodel_p.Vp_Model.Vp),double(LonM),double(LatM),double(DepM));
max_V = nanmax(nanmax(nanmax(Vp_M)));
max_rp = sind(85)/max_V;
ind_ex = find(rp_ref>max_rp);
if ~isempty(ind_ex)
    rp_ref(ind_ex) = [];
    rp_ref = [rp_ref,max_rp];
end

Period = 13;
recordid='RFcomp';

Phase= 'Sp';
cd([basedir,'NewFunctions/PROPMAT/']) % Can we add as path instead??

ilat= 15; ilon=40; % For now. -- which one is close to TUC?
Vp_ray = permute(Vp_M(ilat,ilon,:),[2,3,1]);
Vs_ray = permute(Vs_M(ilat,ilon,:),[2,3,1]);
dV_ray = gradient(Vs_ray)./gradient(depth);
dV_ray(isnan(dV_ray)) = 0;
indv1 = find(~isnan(Vs_ray),1);
indv2 = find(~isnan(Vs_ray),1,'last');
if isempty(indv1)
    Vs_ray(isnan(Vs_ray)) = 3.2;
    Vp_ray(isnan(Vp_ray)) = 5.6;
    Moho_ray = 100000;
else
    Vs_ray(1:indv1) = Vs_ray(indv1);
    Vp_ray(1:indv1) = Vp_ray(indv1);
    if indv2>length(Vs_ray)
        Vs_ray(indv2:end) = Vs_ray(indv2);
        Vp_ray(indv2:end) = Vp_ray(indv2);
    end
    Moho_ray = depth(dV_ray==max(dV_ray(depth>30 & depth<60)));
    if length(Moho_ray)>1
        Moho_ray = Moho_ray(1);
    end
end
vs = Vs_ray(1);
vp = Vp_ray(1);
rho = 2.8;
zbot = [];
ztop = 0;
for idep = 1:length(depth)
    if abs(Vs_ray(idep)-vs(end))/vs(end)>0.005
        vs = [vs;Vs_ray(idep)];
        vp = [vp;Vp_ray(idep)];
        if depth(idep)<Moho_ray
            rho = [rho;2.8];
        else
            rho = [rho;3.3];
        end
        zbot = [zbot;depth(idep)];
        ztop = [ztop;depth(idep)];
    end
end
if ztop(end) ~= depth(end)
    zbot = [zbot;depth(end)];
else
    zbot = [zbot;depth(end)+1];
end
for irp = 1:length(rp_ref)
    write_propmat_syn(ztop,zbot,vp,vs,rho,rp_ref(irp),Period,[recordid,'.mat'],Phase);
    SynWave = load([basedir,'NewFunctions/PROPMAT/',recordid,'.mat']);
    if isnan(SynWave.R(100))
     %   rf_model(irp,ilat,ilon,:) = nan;
        continue;
    end
    [SynRF,SynTime] = getrf(SynWave.R,SynWave.Z,SynWave.dt,vp(1),vs(1),rp_ref(irp));
    [SynDepth,Synindts,Synindte] = getdepth(Phase,SynTime,zbot,vp,vs,rp_ref(irp));
    SynRF_dep = interp1(SynDepth,SynRF(Synindts:Synindte),depth);
    SynRF_dep(isnan(SynRF_dep)) = 0;
    %rf_model(irp,ilat,ilon,:) = permute(SynRF_dep,[1,4,2,3]);
end


figure;
plot(SynRF_dep,depth); 
set(gca,'Ydir','reverse');
xlabel('RF Amplitude'); ylabel('Depth (km)');
title(['Synthetic RF for ',Phase]);