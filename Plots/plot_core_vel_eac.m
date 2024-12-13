%%% This script tries to identify the core velocity of the eac by shifting 
%%% the geogrpahical coordinate to jet cooordinate following the
%%% method of Archer et al., 2017: doi: 10.1002/2017JC013286

clear;

%% Read COF radar data 
datDir = '../DATA/radar-data';
datIn  = fullfile( datDir, ...
    "TUV_COF_20120101T003000Z_20240201T003000Z_2dVar_1-hour-avg.nc" );

glon = ncread( datIn, 'LONGITUDE' );
glat = ncread( datIn, 'LATITUDE' );

[glon,glat] = meshgrid(glon,glat);
glon = glon';glat = glat';

t = ncread( datIn, 'TIME' );
t = t + datenum(1950,01,01);

U = ncread( datIn, 'UCUR' );
V = ncread( datIn, 'VCUR' );

%% Set timeframe 
tStart = datenum(2019,01,01);
tEnd   = datenum(2021,01,01);

%% Cut data within timeframe
indD = find( t >= tStart & t <= tEnd );
t = t(indD);
U = U(:,:,indD);V = V(:,:,indD);
% 
tax = datenum(year(t(1)),month(t(1)),15);
while tax(end) < t(end)
    tax = horzcat(tax,addtodate(tax(end),6,"month"));
end
%% Plot one vel filed for test 
% Aux data for plot
load("../DATA/aux_data/EAC_bathy.mat")
load("../DATA/aux_data/EAC_coastline.mat")
load("../DATA/aux_data/EAC_radarstations.mat")

% For COF
LonLims = [152.9 154.2];
LatLims = [-31.2 -29.8];

s1 = [-33.0109167000000,151.726800000000]; % RHED;
s2 = [-32.4414667000000,152.539083300000]; % SEAL;
s3 = [-29.9838888888889,153.231111111111]; % RRK;
s4 = [-30.6241666666667,153.011388888889]; % NNB;
%%

m_proj('lambert','lon',LonLims,'lat',LatLims);

%% Flags 
meshgrid_flag = 0;

% For COF
vecScl = 40.5;
vecHeadangle=40;
vecShaftwidth=0.025;
vecHeadwidth=0.15;
vecHeadlength=0.15;
dx = 1;

%% Identify jet core 
load_data_flag = 1; % Running the code for identifying jet core takes time. 
                    % The processed results can be loaded from mat file. 

if ~load_data_flag
    lonCoreVel = nan(size(glat,2),length(t));
    latCoreVel = nan(size(glat,2),length(t));
    lonCoreWidth = nan(size(glat,2),2,length(t));
    latCoreWidth = nan(size(glat,2),2,length(t));
    
    coreVel = nan(size(glat,2),length(t));
    coreAngle = nan(size(glat,2),length(t));
    coreWidth = nan(size(glat,2),length(t));
    
    jetxAx = [-1*35*1.5:1.5:35*1.5];
    jetAxVel = nan(size(glat,2),length(jetxAx),length(t));
    
    nbJetAx = 35; % Nb. grid points in one side of the jet axis
    
    for ii=1:length(t)
        disp(['Checking data at ' datestr(t(ii),'yyyy/mm/dd HH:MM')]);
        
        Utmp = squeeze(U(:,:,ii))';Vtmp = squeeze(V(:,:,ii))';
    
        [lonCoreVel(:,ii),latCoreVel(:,ii),lonCoreWidth(:,:,ii),latCoreWidth(:,:,ii),...
            coreVel(:,ii),coreAngle(:,ii),coreWidth(:,ii),jetxAx,jetAxVel(:,:,ii)] = ...
            func_get_core_vel_eac(glon',glat',Utmp,Vtmp,1.5,30,nbJetAx,4);
    end
else
    load('../DATA/processed-data/EAC_core.mat');
end

%% Plot radar vel 
% Time for plotting
tInd = 1824;

h=figure(1);clf;hold on;
set(h, 'Position',[533 403 380 261]);

% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier');
% Coastline
m_plot(clon,clat,'k','Linewidth',1.2);
%
xr = reshape(glon,numel(glon),1);
yr = reshape(glat,numel(glat),1);
uPlot = reshape(U(:,:,tInd),numel(glon),1);
vPlot = reshape(V(:,:,tInd),numel(glon),1);
colormap(cm_delta(32));
m_vec( vecScl*2,xr,yr,uPlot,vPlot,sqrt( uPlot.^2 + vPlot.^2 ),...
            'edgeclip', 'k', 'headangle',vecHeadangle,...
            'shaftwidth',vecShaftwidth,'headwidth',vecHeadwidth,...
            'headlength',vecHeadlength );

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

%
m_grid();clim([0 1.2]);colorbar;
title( [ 'Surface currents at ' datestr(t(tInd), ...
    'dd/mm/yyyy HH:MM') ' UTC' ] );
% % Plot core vel 
m_plot(lonCoreVel(:,tInd),latCoreVel(:,tInd),'rx');

%% 
jetAxVel99 = jetAxVel;
for ix=1:size(jetAxVel99,1)
    for iy=1:size(jetAxVel99,2)

        % if isnan(jetAxVel99(ix,iy,66277))
        %     continue
        % end

        meanVel = nanmean( squeeze(jetAxVel99(ix,iy,:)) );
        meanVelStd = nanstd( squeeze(jetAxVel99(ix,iy,:)) );
        indRem = find( abs(squeeze(jetAxVel99(ix,iy,:))) > meanVel+3*meanVelStd );
        jetAxVel99( ix,iy,indRem ) = NaN;
    end
end

%% 
tmonthStrs = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
%% Linear interpolate to jet coordinate 
h=figure(2);clf;hold on;
set(h,'Position',[564 149 679 620]);

colmap = cm_balance(12);

%
tiledlayout(3,4,"TileSpacing","compact")
for ii=1:12
   nexttile;hold on;
    % Mean
    indT = find(month(t) == ii);
    spd = nanmean(nanmean(squeeze(jetAxVel99(:,:,indT)),3),1);
    spdStd  = nanstd(nanmean(squeeze(jetAxVel99(:,:,indT)),1),[],3);

    jetxAxL = jetxAx(1:nbJetAx);jetxAxR = jetxAx(nbJetAx+2:end);
    indBL = find( spd(1:nbJetAx) < spd(nbJetAx+1) * 0.5 );
    if ~isempty(indBL);bl = jetxAxL(indBL(end));else; bl = NaN;end;
    indBR = find( spd(nbJetAx+2:end) < spd(nbJetAx+1) * 0.5 );
    if ~isempty(indBR);br = jetxAxR(indBR(1));else; br = NaN;end;

    a1 = patch([jetxAx fliplr(jetxAx)], ...
        [spd-spdStd/2 fliplr(spd+spdStd/2)]*-1, ...
        'r');
    set(a1,'FaceAlpha',0.2,'EdgeColor','none');
    plot(jetxAx,spd*-1,'r','LineWidth',1.5);
    plot([bl bl],[-6 3],'k--','LineWidth',2);
    plot([br br],[-6 3],'k--','LineWidth',2);
    plot(jetxAx,spd,'Color',colmap(ii,:),'LineWidth',2.0);
    plot([bl bl],[-6 3],'--','Color',colmap(ii,:),'LineWidth',2);
    plot([br br],[-6 3],'--','Color',colmap(ii,:),'LineWidth',2);

    ylim([-1.8 0]);xlim([-40 50]);

    txt = sprintf('Distance from core (km)');
    xlabel(txt,'FontWeight','bold');
    ylabel('m s^-^1','FontWeight','bold');

    grid on;

    text(10,-1.6,tmonthStrs(ii,:),'FontWeight','bold','FontSize',13);
end

%% Estimate spectrum 
spdOrig = squeeze(nanmean(jetAxVel99(:,nbJetAx+1,:),1));
spdOrig = spdOrig - nanmean(spdOrig);
spdOrig(isnan(spdOrig)) = 0;

nx = length(spdOrig);
nfft = floor(nx/16);
nooverlap = nfft/2;
window = hanning(nx);
p = 0.95;
fs = 24;

vtsRotComp1 = spdOrig; 
vtsRotComp1(isnan(vtsRotComp1)) = 0;
splOrig1   = spectrelisse2(vtsRotComp1,fs,12,p);

errOrig1 = nanmean( splOrig1(3:end,3) );
%% 
% Semi-diurnal Signal
k1 = 1/24*24;
m2 = 1/12*24;
f1 = 23.05/24;
synoptic = [1/20 1/3];
b_monthly = [1/40 1/25];
meso = [1/120 1/60];

h=figure(4);clf;hold on;
set(h,'Position', [680 543 259 260]);

%% Core-vel
plot(splOrig1(3:end-5,1),splOrig1(3:end-5,2),'b','LineWidth',1.5);
% loglog(freq,pxxc,'k--');

a2 = patch([splOrig1(3:end-5,1);flipud(splOrig1(3:end-5,1))], ...
    [splOrig1(3:end-5,2)-splOrig1(3:end-5,3);flipud(splOrig1(3:end-5,2)+splOrig1(3:end-5,3))], ...
    'b');
set(a2,'FaceAlpha',0.2,'EdgeColor','none');

%% Indicators
plot([k1 k1],[2e-2 1e5],'--','Color',[.6 .6 .6]);
plot([m2 m2],[1e-2 1e5],'--','Color',[.6 .6 .6]);
plot([f1 f1],[1e-5 6e-3],'--','Color',[.6 .6 .6]);
plot([2e-1:0.01:1e0],3e-2*[2e-1:0.01:1e0].^(-5/3),'--','Color',[.6 .6 .6]);

% Synotpic 
a = patch([synoptic(1);synoptic(1);synoptic(2);synoptic(2)],[[1e-5;4e1];flipud([1e-5;4e1])] ...
    ,[0.50,0.50,0.50]);
set(a,'FaceAlpha',0.1,'EdgeColor','none');

% Monthly/bi-monthly 
a = patch([b_monthly(1);b_monthly(1);b_monthly(2);b_monthly(2)],[[1e-5;4e1];flipud([1e-5;4e1])] ...
    ,[0.50,0.50,0.50]);
set(a,'FaceAlpha',0.1,'EdgeColor','none');

% Mesoscale
a = patch([meso(1);meso(1);meso(2);meso(2)],[[1e-5;4e1];flipud([1e-5;4e1])] ...
    ,[0.50,0.50,0.50]);
set(a,'FaceAlpha',0.1,'EdgeColor','none');

% Set text label 
text(0.0033,20.7113,sprintf('60 - \n100d'), ...
    'FontWeight','bold','FontSize',9);
text(0.0166,43.2422,sprintf('25 - \n40d'), ...
    'FontWeight','bold','FontSize',9);
text(0.0894,15.4910,sprintf('3 - \n20d'), ...
    'FontWeight','bold','FontSize',9);
text(1.0897,1.5819,sprintf('Diurnal'), ...
    'FontWeight','bold','FontSize',9);
text(2.0885,0.0917,sprintf('Semi-\ndiurnal'), ...
    'FontWeight','bold','FontSize',9);
text(1.2673,0.0003,sprintf('f'), ...
    'FontWeight','bold','FontSize',9);
text(0.3387,0.2551,sprintf('-5/3'), ...
    'FontWeight','bold','FontSize',9);
%%
legend({'Core velocity'},'box','off','Location','southwest');
ylim([5e-5 1e2]);ylabel('Vel. variance (cm^2 s^-^2 cpd^-^1)', ...
    'FontWeight','bold')
xlim([1e-4 1e2]);xlabel('cpd','FontWeight','bold');
set(gca,'XScale','log','YScale','log');