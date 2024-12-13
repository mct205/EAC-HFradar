clear;
%% Read NEWC radar data 
datDir = '../DATA/radar-data';
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
%
maskNEWC = nan(size(glon1)); % Masking grid with data
indnotnan = find( ~isnan(u1(:,:,40656)) );
maskNEWC(indnotnan) = 1;
%% Read COF radar data 
datDir = '../DATA/radar-data';
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
%
maskCOF = nan(size(glon2)); % Masking grid with data
indnotnan = find( ~isnan(u2(:,:,75865)) );
maskCOF(indnotnan) = 1;

%% Load mooring data 
load('../DATA/aux_data/Mooring-data/CH100_uv.mat','tM','uM','vM');
%% Cut data 
tStart = datenum(2012,01,01);
tEnd = datenum(2024,01,01);

% For NEWC
indT1 = find(t1 >= tStart & t1 <= tEnd);
u1 = u1(:,:,indT1); v1 = v1(:,:,indT1);
t1 = t1(indT1);

% For COF
indT2 = find(t2 >= tStart & t2 <= tEnd);
u2 = u2(:,:,indT2); v2 = v2(:,:,indT2);
t2 = t2(indT2);

% For mooring
indTm = find(tM >= tStart & tM <= tEnd);
uM = uM(indTm); vM = vM(indTm);
tM = tM(indTm);

% For tax 
tax = datenum( year(tStart),month(tStart),01,00,00,00 );
while tax(end) < tEnd
    tax = horzcat(tax,addtodate(tax(end),6,'month'));
end

%% Mooring info
% For NEWC
rot_angle1 = 30;
p200  = [152.5342  -33.0019];
p1000 = [152.7309  -33.0013];

% For COF
ch70_loc  = [153.3,-30.275];
ch100_loc = [153.397,-30.268];
c200_loc  = [153.4631,-30.2693];
rot_angle2 = 18;

rotmat1 = [cosd(rot_angle1) -sind(rot_angle1);sind(rot_angle1) cosd(rot_angle1)];
rotmat2 = [cosd(rot_angle2) -sind(rot_angle2);sind(rot_angle2) cosd(rot_angle2)];

%% For filter
fs = 1; % in hours

fc = [1/9360 ]; % 1y low pass
ft = 'low';

%% Extract cross-section mid-point at NEWC
disp('Extracting cross-section at NEWC ...')
indLat1 = 29;indLon1 = 19:37;
plot(glon1(indLon1,indLat1),squeeze(v1(indLon1,indLat1,30234)))

uCross1 = squeeze(u1(indLon1,indLat1,:));
vCross1 = squeeze(v1(indLon1,indLat1,:));

uRotCross1 = nan(size(vCross1));
vRotCross1 = nan(size(vCross1));

% Rotate vectors
for ii=1:size(vCross1,1)
    uts = squeeze(uCross1(ii,:));
    vts = squeeze(vCross1(ii,:));

    velRot = rotmat1 * [uts;vts];

    uRotCross1(ii,:) = velRot(1,:)';
    vRotCross1(ii,:) = velRot(2,:)';
end

%% Extract cross-section along CH100
disp('Extracting cross-section at COF ...')
indLat2 = 100;indLon2 = 20:62;
plot(glon2(indLon2,indLat2),squeeze(v2(indLon2,indLat2,30234)))

uCross2 = squeeze(u2(indLon2,indLat2,:));
vCross2 = squeeze(v2(indLon2,indLat2,:));

uRotCross2 = nan(size(vCross2));
vRotCross2 = nan(size(vCross2));

% Rotate vectors
for ii=1:size(vCross2,1)
    uts = squeeze(uCross2(ii,:));
    vts = squeeze(vCross2(ii,:));

    velRot = rotmat2 * [uts;vts];

    uRotCross2(ii,:) = velRot(1,:)';
    vRotCross2(ii,:) = velRot(2,:)';
end

%% Fitler vel timeseries 
disp('Filtering Velocity timeseries ...')
uRotCrossFilt1 = nan(size(vCross1));
vRotCrossFilt1 = nan(size(vCross1));
uRotCrossFilt2 = nan(size(vCross2));
vRotCrossFilt2 = nan(size(vCross2));

% NEWC
for ii=1:size(vRotCross1,1)
    uRotCrossFilt1(ii,:) = butter_filt(squeeze(uRotCross1(ii,:)),fc,fs,ft);
    vRotCrossFilt1(ii,:) = butter_filt(squeeze(vRotCross1(ii,:)),fc,fs,ft);
end
% COF
for ii=1:size(vRotCross2,1)

    uRotCrossFilt2(ii,:) = butter_filt(squeeze(uRotCross2(ii,:)),fc,fs,ft);
    vRotCrossFilt2(ii,:) = butter_filt(squeeze(vRotCross2(ii,:)),fc,fs,ft);
end

%% Extract velocity in COF and NEWC
% For CH100 
disp('Getting data at CH100....')
d = deg2km( sqrt( (ch100_loc(1)-glon2).^2 + (ch100_loc(2)-glat2).^2 ) );
[ii,jj] = find( d == min(min(d)) );
u2e = squeeze(nanmean(nanmean(u2(ii:ii,jj:jj,:),1),2));
v2e = squeeze(nanmean(nanmean(v2(ii:ii,jj:jj,:),1),2));

velRot2  = rotmat2 * [u2e';v2e'];
velRotM = rotmat2 * [uM';vM'];
ueRot2 = squeeze(velRot2(1,:));veRot2 = squeeze(velRot2(2,:));
ueRotM = squeeze(velRotM(1,:));veRotM = squeeze(velRotM(2,:));

% For mid-point in NEWC 
disp('Getting data at NEWC....')
d = deg2km( sqrt( (p200(1)-glon1).^2 + (p200(2)-glat1).^2 ) );
[ii,jj] = find( d == min(min(d)) );
u1e = squeeze(nanmean(nanmean(u1(ii:ii,jj:jj,:),1),2));
v1e = squeeze(nanmean(nanmean(v1(ii:ii,jj:jj,:),1),2));

velRot1  = rotmat1 * [u1e';v1e'];
ueRot1 = squeeze(velRot1(1,:));veRot1 = squeeze(velRot1(2,:));

% Filtering
ueRot1Filt = butter_filt(ueRot1,fc,fs,ft);
veRot1Filt = butter_filt(veRot1,fc,fs,ft);
ueRot2Filt = butter_filt(ueRot2,fc,fs,ft);
veRot2Filt = butter_filt(veRot2,fc,fs,ft);
ueRotMFilt = butter_filt(ueRotM,fc,fs,ft);
veRotMFilt = butter_filt(veRotM,fc,fs,ft);
%% Test plot 
h=figure(1);clf;hold on;
set(h, 'Position', [371 136 804 662]);
tiledlayout(1,3,"TileSpacing","compact");

% Plot COF
nexttile;hold on;
val2 = vRotCrossFilt2';

yax2 = repmat(t2',size(val2,2),1);
xax2 = repmat(glon2(indLon2,indLat2),1,length(t2));

% Cross-section
pcolor(xax2,yax2,val2');
shading flat;
caxis([-1.2 0.0]);
colormap((flipud(cm_delta(32))));colorbar;
% Bathy location
plot([ch70_loc(1) ch70_loc(1)],[tStart tEnd],'--','Color',[.3 .3 .3],'LineWidth',0.3);
plot([ch100_loc(1) ch100_loc(1)],[tStart tEnd],'--','Color',[.3 .3 .3],'LineWidth',0.3);
plot([c200_loc(1) c200_loc(1)],[tStart tEnd],'--','Color',[.3 .3 .3],'LineWidth',0.3);
% plot([p100(1) p100(1)],[tStart tEnd],'--','Color',[.3 .3 .3],'LineWidth',0.3);
plot([p1000(1) p1000(1)],[tStart tEnd],'--','Color',[.3 .3 .3],'LineWidth',0.3);
%
set(gca,'YDir','reverse');

yticks(tax);yticklabels( datestr(tax,'mmm/yy') );
ylim( [tStart tEnd] );
xlim([153.2 154]);

% Plot NEWC
nexttile;hold on;
val1 = vRotCrossFilt1';

yax1 = repmat(t1',size(val1,2),1);
xax1 = repmat(glon1(indLon1,indLat1),1,length(t1));

% Cross-section
pcolor(xax1,yax1,val1');
shading flat;
caxis([-1.2 0.0]);
colormap((flipud(cm_delta(32))));colorbar;
% Bathy location
plot([p200(1) p200(1)],[tStart tEnd],'--','Color',[.3 .3 .3],'LineWidth',0.3);
plot([p1000(1) p1000(1)],[tStart tEnd],'--','Color',[.3 .3 .3],'LineWidth',0.3);
%
set(gca,'YDir','reverse');

yticks(tax);yticklabels( datestr(tax,'mmm/yy') );
ylim( [tStart tEnd] );
xlim([ 151.8 153]);

% Time-series 
veRot2Filt = veRot2Filt - nanmean(veRot2Filt); % Anomaly
veRot1Filt = veRot1Filt - nanmean(veRot1Filt); % Anomaly
veRotMFilt = veRotMFilt - nanmean(veRotMFilt); % Anomaly

nexttile;hold on;

plot([0 0],[tax(1) tax(end)],'k--');
plot(veRot2Filt,t2,'-o','LineWidth',2,'MarkerSize',2);
plot(veRotMFilt,tM,'LineWidth',3);
plot(veRot1Filt,t1,'-o','LineWidth',2,'MarkerSize',2);
xlim([-0.4 0.6]);
ylim([tStart tEnd]);
yticks(tax);yticklabels(datestr(tax,'mmm/yy'));
xlabel('Vel. (m s^-^1)');
legend({'','COF','CH100','NEWC'},'Orientation','vertical','Box','off');
set(gca,'YDir','reverse');


%% Function
% 4th order Butterworth lowpass filtering
function [vel_filt,meanVel,trndVel]=butter_filt(vel,fc,fs,ft)

% De-mean
meanVel = nanmean(vel);
vel = vel - meanVel;
% Set nan vel to zeros
indnan = find(isnan(vel));
vel(indnan) = 0;
% De-trend
veldt = detrend(vel);

trndVel = vel - veldt;

% Start filtering 
[b,a] = butter(2,fc/(fs/2),ft);
vel_filt = filtfilt(b,a,veldt);

%% Add mean and trend
vel_filt = vel_filt + trndVel + meanVel;

% Set back NaN
vel_filt(indnan) = NaN;

return 
end