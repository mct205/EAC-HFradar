clear;
%% Read NEWC radar data 
datDir = 'DATA/radar-data';
datIn1  = fullfile( datDir, ...
    "TUV_NEWC_20171101T000000Z_20240101T000000Z_2dVar_1-hour-avg.nc" );

glon1 = ncread( datIn1, 'LONGITUDE' );
glat1 = ncread( datIn1, 'LATITUDE' );

[glon1,glat1] = meshgrid(glon1,glat1);
glon1 = glon1';glat1 = glat1';

t1 = ncread( datIn1, 'TIME' );
t1 = t1 + datenum(1950,01,01);

u1 = ncread( datIn1, 'UCUR' );
v1 = ncread( datIn1, 'VCUR' );
%% Read COF radar data 
datDir = 'radar-data';
datIn2  = fullfile( datDir, ...
    "TUV_COF_20120101T003000Z_20240201T003000Z_2dVar_1-hour-avg.nc" );

glon2 = ncread( datIn2, 'LONGITUDE' );
glat2 = ncread( datIn2, 'LATITUDE' );

[glon2,glat2] = meshgrid(glon2,glat2);
glon2 = glon2';glat2 = glat2';

t2 = ncread( datIn2, 'TIME' );
t2 = t2 + datenum(1950,01,01);

u2 = ncread( datIn2, 'UCUR' );
v2 = ncread( datIn2, 'VCUR' );

%% Load mooring data 
load('aux_data/Mooring-data/CH100_uv.mat','tM','uM','vM');
%% Cut data 
tStart1 = datenum(2019,09,01);
tEnd1 = datenum(2021,09,01);

tStart2 = datenum(2019,01,01);
tEnd2 = datenum(2021,01,01);

% For NEWC
indT1 = find(t1 >= tStart1 & t1 <= tEnd1);
u1 = u1(:,:,indT1); v1 = v1(:,:,indT1);
t1 = t1(indT1);

% For COF
indT2 = find(t2 >= tStart2 & t2 <= tEnd2);
u2 = u2(:,:,indT2); v2 = v2(:,:,indT2);
t2 = t2(indT2);

% For mooring
indTm = find(tM >= tStart2 & tM <= tEnd2);
uM = uM(indTm); vM = vM(indTm);
tM = tM(indTm);
%% 
x1 =  [152.5342  -33.0019]; % NEWC mid point
x2 =  [153.3953  -30.2652]; % CH100

% Find indexes 
% NEWC
d = deg2km( sqrt( (x1(1)-glon1).^2 + (x1(2)-glat1).^2 ) );
[ii,jj] = find( d == min(min(d)) );
u1mean = squeeze(nanmean(nanmean(u1(ii-2:ii+2,jj-2:jj+2,:),1),2));
v1mean = squeeze(nanmean(nanmean(v1(ii-2:ii+2,jj-2:jj+2,:),1),2));

% COF
d = deg2km( sqrt( (x2(1)-glon2).^2 + (x2(2)-glat2).^2 ) );
[ii,jj] = find( d == min(min(d)) );
u2mean = squeeze(nanmean(nanmean(u2(ii:ii,jj:jj,:),1),2));
v2mean = squeeze(nanmean(nanmean(v2(ii:ii,jj:jj,:),1),2));

%% Compute domain-average current vel 
% Speed 
spd1mean = sqrt( u1mean.^2 + v1mean.^2 );
spd2mean = sqrt( u2mean.^2 + v2mean.^2 );
spdM     = sqrt( uM.^2 + vM.^2 );

spd1mean = spd1mean - nanmean(spd1mean);
spd2mean = spd2mean - nanmean(spd2mean);
spdM = spdM - nanmean(spdM);
%% X labels
tax1 = datenum(year(tStart1(1)),month(tStart1(1)),day(tStart1(1)));
while tax1(end)<tEnd1
    tax1 = horzcat( tax1,addtodate(tax1(end),1,'month') );
end

tax2 = datenum(year(tStart2(1)),month(tStart2(1)),day(tStart2(1)));
while tax2(end)<tEnd2
    tax2 = horzcat( tax2,addtodate(tax2(end),1,'month') );
end
%% Plot time series 
h=figure(2);clf; hold on;
set(h,'Position', [455 379 580 375]);
colmap = cm_balance(10);

% NEWC
subplot(3,1,1); hold on;
plot(t1,spd1mean,'Color',colmap(3,:));
ylim([-1.5 1.5]);ylabel('NEWC (m s^-^1)','FontWeight','bold');
xticks(tax1);xticklabels(datestr(tax1,'mm/yyyy'));

% COF
subplot(3,1,2); hold on;
plot(t2,spd2mean,'Color',colmap(8,:));
ylim([-1.5 1.5]);ylabel('COF (m s^-^1)','FontWeight','bold');
xticks(tax2);xticklabels(datestr(tax2,'mm/yyyy'));

% Mooring
subplot(3,1,3); hold on;
plot(tM,spdM,'Color',colmap(8,:));
ylim([-1.5 1.5]);ylabel('COF (m s^-^1)','FontWeight','bold');
xticks(tax2);xticklabels(datestr(tax2,'mm/yyyy'));
%% 
spd1 = spd1mean;spd1(isnan(spd1)) = 0;
spd2 = spd2mean;spd2(isnan(spd2)) = 0;
spdM(isnan(spdM)) = 0;
%% Test spectral analysis
nx = length(spd1);
nfft = nx;
window = hamming(4096);
nooverlap = floor(length(window)/2);
p = 0.95;
fs = 24;

splOrig1   = spectrelisse2(spd1,fs,10,p);
splOrig2   = spectrelisse2(spd2,fs,10,p);
splOrigM   = spectrelisse2(spdM,fs,10,p);

errOrig1 = nanmean( splOrig1(3:end,3) );
errOrig2 = nanmean( splOrig2(3:end,3) );
errOrigM = nanmean( splOrigM(3:end,3) );
%% 
k1 = 1/24*24;
m2 = 1/12*24;

f1 = 22.1/24; % f for NEWC
f2 = 23.05/24;% f for COF

synoptic = [1/20 1/3];
b_monthly = [1/60 1/25];
meso = [1/200 1/70];

h=figure(3);clf;hold on;
set(h,'Position', [680 543 518 260]);

% NEWC
subplot(1,2,2);hold on;

% NEWC
plot(splOrig1(2:end-5,1),splOrig1(2:end-5,2),'Color',colmap(3,:),'LineWidth',2);
a = patch([splOrig1(2:end-5,1);flipud(splOrig1(2:end-5,1))], ...
    [splOrig1(2:end-5,2)-splOrig1(2:end-5,3);flipud(splOrig1(2:end-5,2)+splOrig1(2:end-5,3))], ...
    'b');
set(a,'FaceAlpha',0.2,'EdgeColor','none');

% Indicators

% f, tides
plot([k1 k1],[2e-2 1e5],'--','Color',[.6 .6 .6]);
plot([m2 m2],[1e-2 1e5],'--','Color',[.6 .6 .6]);
plot([f2 f2],[1e-5 6e-3],'--','Color',[.6 .6 .6]);
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
text(0.0023,20.7113,sprintf('70 - \n200d'), ...
    'FontWeight','bold','FontSize',9);
text(0.0136,43.2422,sprintf('25 - \n60d'), ...
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

legend({'NEWC'},'box','off','Location','southwest');

ylim([5e-5 1e2]);ylabel('Vel. variance (cm^2 s^-^2 cpd^-^1)', ...
    'FontWeight','bold')
xlim([1e-4 1e2]);xlabel('cpd','FontWeight','bold');
set(gca,'XScale','log','YScale','log');

% COF
subplot(1,2,1);hold on;
% Radar
plot(splOrig2(2:end-5,1),splOrig2(2:end-5,2),'Color',colmap(3,:),'LineWidth',2);
a = patch([splOrig2(2:end-5,1);flipud(splOrig2(2:end-5,1))], ...
    [splOrig2(2:end-5,2)-splOrig2(2:end-5,3);flipud(splOrig2(2:end-5,2)+splOrig2(2:end-5,3))], ...
    colmap(3,:));
set(a,'FaceAlpha',0.2,'EdgeColor','none');

% Mooring
plot(splOrigM(2:end-5,1),splOrigM(2:end-5,2),'Color',colmap(8,:),'LineWidth',2);
a = patch([splOrigM(2:end-5,1);flipud(splOrigM(2:end-5,1))], ...
    [splOrigM(2:end-5,2)-splOrigM(2:end-5,3);flipud(splOrigM(2:end-5,2)+splOrigM(2:end-5,3))], ...
    colmap(8,:));
set(a,'FaceAlpha',0.2,'EdgeColor','none');

% Indicators

% f, tides
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
text(0.0023,20.7113,sprintf('70 - \n200d'), ...
    'FontWeight','bold','FontSize',9);
text(0.0136,43.2422,sprintf('25 - \n60d'), ...
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

legend({'COF','','CH100'},'box','off','Location','southwest');

ylim([5e-5 1e2]);ylabel('Vel. variance (cm^2 s^-^2 cpd^-^1)', ...
    'FontWeight','bold')
xlim([1e-4 1e2]);xlabel('cpd','FontWeight','bold');
set(gca,'XScale','log','YScale','log');