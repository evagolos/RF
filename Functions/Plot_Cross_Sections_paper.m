clear all
% So these are the data coverage matrices that I saved after the fact with
% Ved's script.  They're the same size as plane_RF.
Data = '/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin'; 
fol = 'Data/Projects/SWUS_big/CCP';
cd(Data)
addpath([Data '/Functions'])
Ps= 0;
if Ps==1
    file='CCP_Ps_13_SWUS_Ps_4_100s_05.10.20_1_780km_deg_1.mat';
else
    %file='CCP_Sp_13_SWUS_Sp_4_100s_05.10.20_1_2_450km_deg_1.mat';
    %file='CCP_Sp_13_SWUS_big_Sp_4_100s_18.05.21_450km_deg_1.mat';
    %file='CCP_Sp_13_SWUS_big_Sp_2_100s_22.04.22_450km_deg_1.mat';
    file='CCP_Sp_13_SWUS_big_Sp_2_100s_22.04.22_450km_deg_2.mat';
end
    
load([Data '/' fol '/' file]);
junk=strsplit(file,[phase '_']); junk=strsplit(junk{end},'_'); T=junk{1};
string = ['pics_' T 's_'];
VsModelName='MITPS_10kmZ';
%VsModelName='DNA13';
%VsModelName='FichtnerSV13';
load([Data '/Data/Velocity_Models/Mantle/Vs/' VsModelName '.mat']);

[~,~,Vs_Model.Vs] = gradient(Vs_Model.Vs,1,1,Vs_Model.Depth(2)-Vs_Model.Depth(1));
% 
% 
plane_RF=real(plane_RF);  %CHANGE
%% Interpolate data
min_dep = 0; 
if Ps==1
    max_dep = 780; 
else
    max_dep= 450;
end
DepNpts = (max_dep-min_dep)/5;
rplane_rf=reshape(plane_RF,[size(tmplat), length(cp_depths)]);
rplane_std=reshape(std_RF,[size(tmplat), length(cp_depths)]);
%weighted_n_events = weighted_n_events./repmat(nansum(weighted_n_events,1),size(weighted_n_events,1),1)*nansum(nansum(weighted_n_events,1))/size(weighted_n_events,2);
weighted_n_events = weighting_norm;
rweighted_n_events=reshape(weighted_n_events,[size(tmplat), length(cp_depths)]);
rplane_rf(rweighted_n_events<0.05) = NaN;
rplane_rf(rplane_std>0.2 & 3*abs(rplane_std)>abs(rplane_rf)) = NaN;
%rweighted_n_events=reshape(n_events,[size(tmplat), length(cp_depths)]);
%rstd_RF=reshape(std_RF,[size(tmplat),length(cp_depths)]);
deps_tmp = linspace(min_dep,max_dep,DepNpts); deps=double(cp_depths)';
[glo,gla,gdep]=ndgrid(tmplon(:,1),tmplat(1,:)',deps);
[plo,pla]=ndgrid(tmplon(:,1),tmplat(1,:)');
RF=griddedInterpolant(glo,gla,gdep,rplane_rf);
%STD=griddedInterpolant(glo,gla,gdep,rplane_std);
COV=griddedInterpolant(glo,gla,gdep,rweighted_n_events);



% [vellat,vellon]=meshgrid(Vs_Model.Latitude,Vs_Model.Longitude);
% [vlo,vla,vdep]=ndgrid(vellon(:,1),vellat(1,:)',Vs_Model.Depth);
% vsmo=Vs_Model.Vs;
% P=[2,1,3];
% vlo=permute(vlo, P);
% vla=permute(vla, P);
% vdep=permute(vdep, P);
% velu=interp3(vlo,vla,vdep,vsmo,glo,gla,gdep);
% VEL=griddedInterpolant(glo,gla,gdep,velu);




%% set up cross section lines

%close all;
%sects='randomW';%'NW-SE'; %'NW-SW'; % 'SE-NW'; %  'W-E'; %
% nsects=30; xsl=5;
% iftrans=1;
% dlat=1;
% dlon=0;
% direc=[dlat,dlon];
% ntrans=5;
% direc=direc/ntrans;

%latlons=[33 -115 33 -105; 35 -115 35 -105; 37 115 37 -105; 39 -115 39 -105; 41 -115 41 -105];
latlons=[34 -117 34 -104; 36 -117 36 -104; 38 117 38 -104; 42 -117 42 -104];
k= size(latlons,1);
kj=round(k/2); sc = 0.04;%0.15;%0.06;

Npts=ceil(111.16*distance(latlons(kj,1),latlons(kj,2), ...
    latlons(kj,3), latlons(kj,4))/10);
latsAll=zeros(size(latlons,1),Npts); lonsAll=latsAll;
for k=1:size(latlons,1);
    lonsAll(k,:)= linspace(latlons(k,2),latlons(k,4),Npts);
    latsAll(k,:)=linspace(latlons(k,1),latlons(k,3),Npts);
end

PLOTNOW=1:size(latlons,1);

%% Generate the figure


fig1= figure();  
ppos= get(fig1,'Position');
set(fig1,'Position',[ppos(1)-500 ppos(2) ppos(3)*2 ppos(4)*2.3]);
    
% Load text files with boundaries of geologic provinces.
BNR= importdata('/Users/evagolos/Research/Misc/GeoBounds/BNR.txt');
Col= importdata('/Users/evagolos/Research/Misc/GeoBounds/Colorado.txt');
RGR= importdata('/Users/evagolos/Research/Misc/GeoBounds/RioGrande.txt');
SNB= importdata('/Users/evagolos/Research/Misc/GeoBounds/Sierra.txt');



for kj =1:length(PLOTNOW)
%for kj=18
    clear x y z;
    clear x y z; clear discmdl n_eventsINTERP weighted_eventsINTERP Overlay;
    lons = lonsAll(kj,:); lats = latsAll(kj,:);
    [gdep1,glo1]=ndgrid(deps,lons); [~,gla1]=ndgrid(deps,lats);

     latlima=max(latlons(kj,1),latlons(kj,3));
     latlimi=min(latlons(kj,1),latlons(kj,3));
     lonlima=max(latlons(kj,2),latlons(kj,4));
     lonlimi=min(latlons(kj,2),latlons(kj,4));
     coplat=99;
     coplon=[];

     

    %interpolate plane_RF onto grid of cross sections
     discmdl=RF(glo1,gla1,gdep1);
%     %discmdlstd=STD(glo1,gla1,gdep1);
%     discmdlvel=VEL(glo1,gla1,gdep1);
    weighted_eventsINTERP=COV(glo1,gla1,gdep1);
    weighted_eventsINTERP(weighted_eventsINTERP<0)= NaN; % Otherwise get errors.
    
    
    dalj = 111.16*distance(latlons(kj,1),latlons(kj,2),lats,lons);
     if coplat~=99
         dalcop= 111.16*distance(latlons(kj,1),latlons(kj,2),coplat,coplon);
     end
    
    Cutoff=50;
    weighted_or_not = weighted_eventsINTERP;
    
    fontsz=20;
    
    % plot cross sections and coverage
    %close all
    savefig='no';
    min_dep = 0;   
    max_dep= 230;
    maxlat=max(Vs_Model.Latitude)-0.1; %sc=0.25;
    clear x y z %close all; clear x y z;
    
    %f=figure(kj);set(f,'position',[1 41 1600 783], 'PaperPositionMode','auto');
    Title=axes('position',[0 0 1 1],'visible','off');
    switch (lower(phase)); case 'sp'; cmap=fliplr(colormap(jet)); case 'ps'; cmap=colormap(jet); end
    text(0.35, 0.95,['\fontsize{14} \bf ' num2str(latlons(kj,1)) 'N, ' num2str(latlons(kj,2)) 'E   to   ' ...
        num2str(latlons(kj,3)) 'N, ' num2str(latlons(kj,4)) 'W'])
    
    if max_dep==75; fpos=[60 350 1480 750]; f2pos=[60 120 1480 300];
    else f2pos=[985 555 920 307]; fpos=[55 265 920 820];
    end
    
    % Set up coordinates for profile line
    if lats(1)>lats(end)
        dalj2=lats(1):-0.5:lats(end);
    else
        dalj2=lats(1):0.5:lats(end);
    end
    halfdeg=abs((dalj(end)-dalj(1))/diff(dalj2([1 end]))/2);
    lims=[min(dalj) max(dalj)];
    
    
    % This section finds intersection points of the geologic provinces to
    % plot.
%     P_BNR= InterX(BNR',[lons; lats]); nBNR= size(P_BNR,2);
%     PointPlot= [];
%     for ii=1:nBNR
%         dlatlon= sqrt(((lons-P_BNR(1,ii)).^2 - (lats-P_BNR(2,ii)).^2));
%         [~,minD]= min(dlatlon);
%         PointPlot= [PointPlot; minD];
%     end
%     P_Col= InterX(Col',[lons; lats]); nCol= size(P_Col,2);
%         for ii=1:nCol
%         dlatlon= sqrt(((lons-P_Col(1,ii)).^2 - (lats-P_Col(2,ii)).^2));
%         [~,minD]= min(dlatlon);
%         PointPlot= [PointPlot; minD];
%     end
%     P_RGR= InterX(RGR',[lons; lats]); nRGR= size(P_RGR,2);
%         for ii=1:nRGR
%         dlatlon= sqrt(((lons-P_RGR(1,ii)).^2 - (lats-P_RGR(2,ii)).^2));
%         [~,minD]= min(dlatlon);
%         PointPlot= [PointPlot; minD];
%     end
%     P_SNB= InterX(SNB',[lons; lats]); nSNB= size(P_SNB,2);
%         for ii=1:nSNB
%         dlatlon= sqrt(((lons-P_SNB(1,ii)).^2 - (lats-P_SNB(2,ii)).^2));
%         [~,minD]= min(dlatlon);
%         PointPlot= [PointPlot; minD];
%     end
    
% Plot CCP stack profile
    a1=axes('position',[0.03 0.49-(kj-1)*0.15 0.8 0.3]); hold on
    % Why are these getting smaller each time?
    h = imagesc(dalj,deps,discmdl); caxis(sc*[-1.3 1.3]);
    set(h,'alphadata',~isnan(discmdl));
    ylabel('Depth (km)','fontsize',10);
    %plot(dalj,5*topo_prof,'color',[0 0 0],'linewidth',2);
    % Line markers at surface
    plot(a1,dalj(1:10:end),0.*dalj(1:10:end)-20,'ko','markersize',6,'markerfacecolor','g'); 
     if coplat ~= 99
         plot(a1,dalcop,0,'ko','markersize',6,'markerfacecolor','r');
     end
     % Intersection points
%      for jj=1:length(PointPlot)
%         plot(a1,[dalj(PointPlot(jj)) dalj(PointPlot(jj))],[-10 max(deps)],'--k','linewidth',2);
%      end
    set(a1,'ydir','reverse','ytick',0:50:max_dep,'yticklabel',0:50:max_dep); daspect([1 1 1]);
    set(a1,'xtick',dalj(1):halfdeg:dalj(end),'xticklabel',dalj2)
    set(a1,'ygrid','on','xgrid','on','layer','top','fontsize',10,'box','on');
    axis([lims -25 max_dep]);   ylabel('Depth (km)','fontsize',10);
    %title(['a) CCP Stack at ' T '-100s Filter Band']);
    xlabel('X Distance (degrees)');
    colormap(cmap); c=colorbar('location', 'southoutside');
    c.Label.String= 'Stack Amplitude'; c.Label.FontSize=12;
    set(a1,'FontSize',12);
    
end
    
figure();
    % Plot map view, colored by weighting.
    a3=axes('position',[0.72 0.49 0.25 0.4]); hold on
    lon=[-120 -100]; lat=[28 42];
    fig2= figure;
    daspect([111.16, 111.16*distance(mean(lat),0,mean(lat),1), 1]);
    box on;  axis([lon, lat]);%[-90 -77, lat]);
    contourf(tmplon(:,1),tmplat(1,:),log(rweighted_n_events(:,:,200))','linestyle','none');
    grid on; %set(a3,'xgrid','on','ygrid','on','layer','top')
    
    % Add state boundaries, geologic boundaries to map
        S=shaperead('usastatelo','UseGeoCoords',true);
        for i=3:51; plot(S(i).Lon,S(i).Lat,'k'); end
        plot(Col(:,1),Col(:,2),'k','linewidth',2);
        text(248-360,37,'Colo','FontWeight','bold','FontSize',12);
        
        plot(RGR(:,1),RGR(:,2),'k','linewidth',2);
        text(256-360,34,'RG','Fontweight','bold');
        plot([255 256]-360,[33 33.5],'k');
        plot(BNR(:,1),BNR(:,2),'k','linewidth',2);
        text(241-360,40.1,'BNR','FontWeight','bold');
        plot(SNB(:,1),SNB(:,2),'k','linewidth',2);
        text(-120.6,38.1,'SNB','FontWeight','bold');
    
    for ii=1:size(latlons,1)
        lon= [min([latlons(ii,2) latlons(ii,4)])-2, max([latlons(ii,2) latlons(ii,4)+2])]; 
        lat= [min([latlons(ii,1) latlons(ii,3)])-2, max([latlons(ii,1) latlons(ii,3)+2])]; 
        %lat=[min(lats)-2 max(lats)+2]; lon=[min(lons)-2 max(lons)+2];%[26 40];
        if abs(diff(lat))<10; lat=mean(lat)+[-5 5]; end
        if abs(diff(lon))<6; lon=mean(lon)+[-5 5]; end
        plot(lons,lats,'k-','linewidth',2); hold on;
        plot(lons(1:10:end),lats(1:10:end),'ko','markerfacecolor','c','markersize',8);
    end
    title('Map View, Weighting');
    coast=load('coast');
    hold on
    plot(coast.long,coast.lat)
    xlabel('Longitude (degrees)','FontSize',12); ylabel('Latitude (degrees)','FontSize',12);
    c=colorbar('location', 'southoutside'); 
    set(a3,'FontSize',12);
    c.Label.String= 'log(weighting) at 200 km'; c.Label.FontSize=12;
    xlim([lon(1) lon(2)]); ylim([lat(1) lat(2)]);
    
    % Plot weighting in vertical profile view.
%     a4=axes('position',[0.72 0.04 0.25 0.4]); hold on
%     imagesc(dalj,deps,log(weighted_eventsINTERP));% caxis([0 2.5]);
%     ylabel('Depth (km)','fontsize',10);
%     %plot(dalj,5*topo_prof,'color',[0 0 0],'linewidth',2);
%     plot(a4,dalj(1:10:end),0.*dalj(1:10:end)-20,'ko','markersize',6,'markerfacecolor','g');
%     if coplat ~= 99
%          plot(a2,dalcop,0,'ko','markersize',6,'markerfacecolor','r');
%      end
%     set(a4,'ydir','reverse','ytick',0:50:max_dep,'yticklabel',0:50:max_dep); daspect([1 1 1]);
%     set(a4,'xtick',dalj(1):halfdeg:dalj(end),'xticklabel',dalj2,'fontsize',11)
%     set(a4,'ygrid','on','xgrid','on','layer','top','fontsize',10,'box','on');
%     axis([lims -25 max_dep]);   ylabel('Depth (km)','fontsize',10);
%     title(['Depth View, Weighting']);
%     xlabel('X Distance (degrees)');
%     c=colorbar('location', 'southoutside');
%     c.Label.String='log(weighting)'; c.Label.FontSize=12;
%     caxis([0,2.5])
%     caxis([0,4])
%     
%     if strcmp(savefig, 'yes')
%         
%         pause(1); figure(f); fext='png';
%         if strcmp(fext,'ps'); pstr='epsc'; elseif strcmp(fext,'jpg'); pstr='jpeg'; else pstr=fext; end
%         eval(['print -d' pstr ' ' figname '.' fext]);
%         if strcmp(fext,'ps'); system(['ps2pdf -dEPSCrop ' figname '.' fext]);
%             system(['rm ' figname '.' fext]);   end
%     end
%end

