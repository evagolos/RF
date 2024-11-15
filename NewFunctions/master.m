%% Master code to make CCP stacks
clear all; close all; clc
basedir='/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin/'; % Change this for your file system.
cd(basedir);
javaaddpath([basedir 'Functions/taup/lib/TauP-1.1.7.jar']);
javaaddpath([basedir 'Functions/taup/lib/log4j-1.2.8.jar']);
javaaddpath([basedir 'Functions/taup/lib/seisFile-1.0.1.jar']);
addpath([basedir 'Functions/Matlab_TauP']); addpath([basedir 'NewFunctions']);
addpath([basedir 'Functions']);
addpath([basedir 'NewFunctions/SACfun']);
addpath([basedir 'NewFunctions/PROPMAT']);
addpath(genpath([basedir 'seizmo']),'-end');


Project = 'Alaska_test';

%% PREP: Request data
%First, make sure NEIC Catalogue is up to date
% neic='no';
% if strcmp(neic, 'yes'); Auto_Prep(Project,basedir,'neic'); end; clear neic

% Auto_Prep(Project,basedir,'request') % Request workflow needs to be
% overhauled, don't use.

%% PREP: Unpack SAC Files
%Auto_Prep(Project, basedir, 'unpack'); % Similarly, don't use until fix is pushed.
% check for and remove repeat data
Auto_Prep(Project,basedir,'multiples'); % This generates Network.mat file.

%% Get all file names

Proj_Dir = [basedir 'Data/Projects/' Project '/'];
load([Proj_Dir,'Networks.mat']);
Nets = fieldnames(Networks);
cps = {'H1','H2','HZ','HN','HE'};
for in = 1:length(Nets)
   disp(Nets{in});
   Stas = fieldnames(Networks.(Nets{in}));
   for is = 1:length(Stas)
       sac_files = dir([Proj_Dir,Nets{in},'/',Stas{is},'/SAC_Files/*SAC']);
       SAC_Filenames.(Nets{in}).(Stas{is}).mark = cell(length(sac_files),1);
       SAC_Filenames.(Nets{in}).(Stas{is}).name = cell(length(sac_files),1);
       for ifile = 1:length(sac_files)
           sac_comp = sac_files(ifile).name(end-9:end);
           ic = [];
           for icp = 1:length(cps)
               if contains(sac_comp,cps{icp})
                   ic = [ic,icp];
               end
           end
           if length(ic)~=1
               disp(['error ',Nets{in},' ',Stas{is},' ',sac_comp])
               continue;
           end
           SAC_Filenames.(Nets{in}).(Stas{is}).mark(ifile) = {[sac_files(ifile).name(1:11),cps{ic}]};
           SAC_Filenames.(Nets{in}).(Stas{is}).name(ifile) = {sac_files(ifile).name};
       end
       for ifile = length(sac_files):-1:1
           if isempty(SAC_Filenames.(Nets{in}).(Stas{is}).mark{ifile})
               SAC_Filenames.(Nets{in}).(Stas{is}).mark(ifile) = [];
               SAC_Filenames.(Nets{in}).(Stas{is}).name(ifile) = [];
           end
       end   
       [SAC_Filenames.(Nets{in}).(Stas{is}).mark,inds] = unique(SAC_Filenames.(Nets{in}).(Stas{is}).mark);
       SAC_Filenames.(Nets{in}).(Stas{is}).name = SAC_Filenames.(Nets{in}).(Stas{is}).name(inds);
   end
end
save([Proj_Dir,'SAC_Filenames.mat'],'SAC_Filenames');

%% PREP: Update pre-existing project
% First, make sure NEIC Catalogue is up to date
neic='no';
if strcmp(neic, 'yes'); Auto_Prep(Project,basedir,'neic'); end; clear neic


%% PREP: Run Array picker
% This takes FOREVER so specify a network to run in parallel (by hand)
% If array picker not used, signal-to-noise ratio picker is used instead.

lowT = 1; % s - period of high frequency corner
highT = 33; % You can change lowT, highT
%nets = 'all'; % Default: run all networks sequentially
nets= {'YV'};

Auto_Prep_RT(Project,basedir,'arrayPick', nets,lowT,highT);

%% PREP: Once all networks have been array picked, amalgamate and reformat
Auto_Prep_RT(Project,basedir,'arrayPickReformat',lowT,highT);

%% PREP: Prep Waveforms
% Hardwire in networks to prep, or can input by running code if these
% variables are set to empty
NETs={'YV'};
%NETs = []; % all


% Set Free Surface Velocity correctioxn as Fixed (set Vp, Vs) or Best
% Doesn't matter if we use RT
FS='Best'; % 'Fixed';
if strcmp(FS,'Fixed'); Vp=5.7; Vs=3.2; else Vp=[]; Vs=[]; end

% I don't use the following line.
%Auto_Prep_free_surface(Project,basedir,'prep', NETs, FS, Vp, Vs, 'para');

% Use this line.
Auto_Prep_RT(Project,basedir,'prep', NETs, FS, Vp, Vs,lowT,highT);
%% RFs:  Set variables for RF calculation
% Hardwire in networks to generate RFs for, or can input by running code if
% these variables are set to empty
%NETs=[];% % 
NETs= {'YV'};

% Set high frequency corner and migration models
lowT = 1; % s - period of high frequency corner
highT = 33;
CrustModel  = 'None'; % - Migration models
MantleVpModel = 'AKAN2020'; % Model for Alaska
MantleVsModel = 'AKAN2020'; 
SvSh='SV'; % 'SH':  P/SH coupling;  'S  V':  P/SV coupling
isRPdiff = 'no';

version=['names_' num2str(lowT) '_' num2str(highT)];
phase= {'Ps'}; % or {'Sp'}
%% RFs:  Generate RFs
Generate_RFs(Project, basedir, NETs, phase, 'generate', lowT,highT,SvSh,CrustModel,MantleVpModel,MantleVsModel,isRPdiff,version);
% You can preliminarily examine RFs after this step.
% % To fully complete this step, you'll need to use the PROPMAT software, in directory NewFunctions. Might not be
% % able to compile this on new Mac chips.
%% Generate Single Station Stack in Depth
% Again, will need PROPMAT to complete this step.

% lowt=lowT;
% %NETs=[];  
% NETs= {'YV'};
% ifplot='yes'; % Spit out plots of single-station stacks
% ifbaz='no'; % If you want to make stacks at a particular backazimuth range (need a lot of events)
% Generate_RFs(Project, basedir, NETs, phase, 'SingleStaRF', lowT,highT,ifplot,ifbaz,isRPdiff,version);


%% CCP:  Calculate depth/lat/lon matrix
version='names';
phase = {'Ps'};
isRPdiff='no';

if strcmp(isRPdiff,'no')
    CCP_stack_slopeangle(Project, basedir, phase, 'calc_deplatlon',version);
    
    % CCP: Calculate stack
    CCP_stack_slopeangle(Project, basedir, phase, 'stack',version)
else
    vmodel.Vp = 'DNA13';
    vmodel.Vs = 'DNA13';
    vmodel.Cr = 'None';
    CCP_stack_slopeangle_RPdiff(Project, basedir, phase, 'calc_deplatlon',vmodel,version);
    % CCP: Calculate stack
    CCP_stack_slopeangle_RPdiff(Project, basedir, phase, 'stack',vmodel,version)
end


