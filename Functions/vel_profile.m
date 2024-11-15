function vel = vel_profile(PS,dep,lat,lon,elev,moh,sed)
% vel = vel_profile(PS,dep,lat,lon,elev,moh,sed)
%
% For elev?0 (land), use crust 1.0 for the crust.
% For elev<0 (ocean) use standard oceanic crust, interp sediments.
% For the mantle velociites, it uses iasp91. 

if nargin < 3
    model = iasp91('depths',dep(:));
    if PS==1, vel = model.vp; end
    if PS==2, vel = model.vs; end
    return
end

if nargin < 5 || isempty(elev)
    elev = 0;
end
if nargin < 6 || isempty(moh)
    moh = false;
end
if nargin < 7 || isempty(sed)
    sed = 0;
end

dep=dep(:);
N = length(lat);

% ocean profile
zo  = [0    1.15    1.16	2.32    2.33    6.86    6.87];
vpo = [4.7  4.7     6.3     6.3     6.9     6.9     7.55];
vso = [2.54 2.54    3.59    3.59    3.94    3.94    4.35];

% seds 
sed_vp = 1.6;
sed_vs = 0.3;

vel = nan(length(dep),N);

% mantle
model = iasp91('depths',dep(dep>0),'crust',false);



for ii = 1:N

if elev >= 0 % land station
    %% CRUST
    crust = getcrust(lat(ii),lon(ii));
    % massage crustal layers so they have constant V in layers and doublets
    % a metre apart at the discontinuities, and NO sub-moho V
    glay = (crust.thk > 0) & [0 0 1 1 1 1 1 1];
    zc = cumsum(crust.thk(glay));
    % to include CRUST1.0 submoho velocity
%     zc = [0,zc(1:end-1)+1e-3,zc,zc(end)+1e-3];
%     vpc = [crust.vp(glay),crust.vp(glay),crust.vp(end)];
%     vsc = [crust.vs(glay),crust.vs(glay),crust.vs(end)];
    % to get the sub-moho velocity from the mantle model
    zc = [0,zc(1:end-1)+1e-3,zc];
    vpc = [crust.vp(glay),crust.vp(glay)];
    vsc = [crust.vs(glay),crust.vs(glay)];
    [~,ic] = sort(zc);
    
    zc = zc(ic)';
    vpc = vpc(ic)';
    vsc = vsc(ic)';
    
    zc = zc-elev; % reset top as the topog.
    if moh
        zc(end) = moh;
    end
       
elseif elev < 0 % obs sta
    
    %% SEDS
    zs = [0 sed]; 
    vps = [sed_vp sed_vp];
    vss = [sed_vs sed_vs];
    
    %% CRUST
    zc = [zs zo+sed]'-elev;
    vpc = [vps vpo]';
    vsc = [vss vso]';    
end

%% CAT
zz = [zc; model.depth(model.depth>max(zc))];
vp = [vpc; model.vp(model.depth>max(zc))];
vs = [vsc; model.vs(model.depth>max(zc))];

% silly thing to avoid repeated depths
while any(diff(zz)==0)
    zz = zz + [0;1e-10*double(diff(zz)==0)];
end 
    

if PS==1
    V = vp;
elseif PS==2
    V = vs;
end

% interp to dep vals
vel(:,ii) = interp1(zz,V,dep);


end % loop on points

end








