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
%% 
tax = datenum(2013,01,01):1:datenum(2013,12,31); % Abitrarily choosing so 
                                               % we can have daily 
                                               % timestamp for a random year
taxMonthly = tax(1);
while taxMonthly(end) < tax(end)
    taxMonthly = vertcat(taxMonthly,addtodate(taxMonthly(end),1,'month'));
end
taxMonthly(end) = [];
%%
tmonths = [1:12];
tmonthStrs = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

nDat1 = squeeze( sum( sum( isfinite(u1), 2 ), 1 ) );
nDat2 = squeeze( sum( sum( isfinite(u2), 2 ), 1 ) );

iDatAvail1 = zeros(1,length(t1));
iDatAvail1(nDat1~=0) = 1;

iDatAvail2 = zeros(1,length(t2));
iDatAvail2(nDat2~=0) = 1;

nDatMonths1 = zeros(1,length(tmonths));
nDatMonths2 = zeros(1,length(tmonths));

for ii=tmonths
    ind = find(month(t1) == ii);
    nDatMonths1(ii) = sum( iDatAvail1(ind) );

    ind = find(month(t2) == ii);
    nDatMonths2(ii) = sum( iDatAvail2(ind) );
end

%% Plot nb of data each month
h=figure(1);clf;hold on;
set(h,'Position',[455 453 696 301]);

% NEWC
subplot(2,1,1);hold on;
bar(tmonths,nDatMonths1);
xticks(1:12);
ylabel('NEWC','FontWeight','bold');

% COF
subplot(2,1,2);hold on;
bar(tmonths,nDatMonths2);
xticks(1:12);
ylabel('COF','FontWeight','bold');

%% Calculate radar vel each month 
Umonthly1 = nan(size(glon1,1),size(glon1,2),length(tmonths));
Vmonthly1 = nan(size(glon1,1),size(glon1,2),length(tmonths));

Umonthly2 = nan(size(glon2,1),size(glon2,2),length(tmonths));
Vmonthly2 = nan(size(glon2,1),size(glon2,2),length(tmonths));

SpdMonthly1 = nan(size(glon1,1),size(glon1,2),length(tmonths));
SpdMonthly2 = nan(size(glon2,1),size(glon2,2),length(tmonths));

fmaj1 = nan(size(glon1,1),size(glon1,2),length(tmonths));
fmin1 = nan(size(glon1,1),size(glon1,2),length(tmonths));
finc1 = nan(size(glon1,1),size(glon1,2),length(tmonths));

fmaj2 = nan(size(glon2,1),size(glon2,2),length(tmonths));
fmin2 = nan(size(glon2,1),size(glon2,2),length(tmonths));
finc2 = nan(size(glon2,1),size(glon2,2),length(tmonths));

for ii=tmonths

    disp(['Computing NEWC mean vel for month ' num2str(ii) '...']);

    % For NEWC 
    ind = find(month(t1) == ii & iDatAvail1' == 1);
    Umonthly1(:,:,ii) = nanmean(u1(:,:,ind),3).*maskNEWC;
    Vmonthly1(:,:,ii) = nanmean(v1(:,:,ind),3).*maskNEWC;
    
    SpdMonthly1(:,:,ii) = sqrt( Umonthly1(:,:,ii).^2 + Vmonthly1(:,:,ii).^2 );

    % Compute var ellipses for NEWC
    [nx,ny] = size(glon1);
    for xind=1:nx
        for yind=1:ny
            if isnan(Umonthly1(xind,yind,ii))
                continue
            end
            % 
            [fmaj1(xind,yind,ii),fmin1(xind,yind,ii),finc1(xind,yind,ii)] = ...
                get_PCA( squeeze(u1(xind,yind,ind)), squeeze(v1(xind,yind,ind)) );
        end
    end
    
    disp(['Computing COF mean vel for month ' num2str(ii) '...']);

    % For COF 
    ind = find(month(t2) == ii & iDatAvail2' == 1);
    Umonthly2(:,:,ii) = nanmean(u2(:,:,ind),3).*maskCOF;
    Vmonthly2(:,:,ii) = nanmean(v2(:,:,ind),3).*maskCOF;
    
    SpdMonthly2(:,:,ii) = sqrt( Umonthly2(:,:,ii).^2 + Vmonthly2(:,:,ii).^2 );

    % Compute var ellipses for COF
    [nx,ny] = size(glon2);
    for xind=1:nx
        for yind=1:ny
            if isnan(Umonthly2(xind,yind,ii))
                continue
            end
            % 
            [fmaj2(xind,yind,ii),fmin2(xind,yind,ii),finc2(xind,yind,ii)] = ...
                get_PCA( squeeze(u2(xind,yind,ind)), squeeze(v2(xind,yind,ind)) );
        end
    end
end

%% Now for plot 
% Aux data for plot
load("../DATA/aux_data/EAC_bathy.mat")
load("../DATA/aux_data/EAC_coastline.mat")
load("../DATA/aux_data/EAC_radarstations.mat")

% For both sites
LonLims = [150.5 154.3];
LatLims = [-34.0 -29.8];
% For NEWC
LonLims1 = [151.6 153.1];
LatLims1 = [-33.6 -32.3];
% For COF
LonLims2 = [152.7 154.2];
LatLims2 = [-31.2 -29.8];



%% For plotting monthly NEWC
m_proj('lambert','lon',LonLims1,'lat',LatLims1);

vecScl = 8.5;
vecHeadangle=50;
vecShaftwidth=0.025;
vecHeadwidth=0.25;
vecHeadlength=0.33;
dx = 1;

% For Ellipses
dEll = 3;

ellScl = 42000;

% For plotting NEWC
h=figure(2);clf;hold on;
set(h, 'Position',[316 182 1090 659]);
tiledlayout(3,4,"TileSpacing","compact");

% Plot radar vel 
for indPlot=1:12
    
    switch indPlot
        % Austral summer 
        case 12
            figPos = 1;
        case 1 
            figPos = 5;
        case 2
            figPos = 9;

        % Austral autumn 
        case 3 
            figPos = 2;
        case 4 
            figPos = 6;
        case 5 
            figPos = 10;

        % Austral winter 
        case 6 
            figPos = 3;
        case 7 
            figPos = 7;
        case 8 
            figPos = 11;

        % Austral spring 
        case 9 
            figPos = 4;
        case 10 
            figPos = 8;
        case 11 
            figPos = 12;

    end
    ax=nexttile(figPos);hold on;
    % Bathy lines 
    [C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
        'ShowText','off','Color',[.6 .6 .6]);
    clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier');
    % Coastline
    m_plot(clon,clat,'k','Linewidth',1.2);

    % Vectors
    xr = reshape(glon1(1:dx:end,1:dx:end),numel(glon1(1:dx:end,1:dx:end)),1);
    yr = reshape(glat1(1:dx:end,1:dx:end),numel(glat1(1:dx:end,1:dx:end)),1);

    uPlot = reshape(Umonthly1(1:dx:end,1:dx:end,indPlot),numel(glon1(1:dx:end,1:dx:end)),1);
    vPlot = reshape(Vmonthly1(1:dx:end,1:dx:end,indPlot),numel(glon1(1:dx:end,1:dx:end)),1);
    colormap(cm_delta(32));
    m_vec( vecScl*2,xr,yr,uPlot,vPlot,sqrt( uPlot.^2 + vPlot.^2 ),...
                'edgeclip', 'k', 'headangle',vecHeadangle,...
                'shaftwidth',vecShaftwidth,'headwidth',vecHeadwidth,...
                'headlength',vecHeadlength );
    
    % Ellipses
    xre2 = reshape(glon1(1:dEll:end,1:dEll:end),numel(glon1(1:dEll:end,1:dEll:end)),1);
    yre2 = reshape(glat1(1:dEll:end,1:dEll:end),numel(glon1(1:dEll:end,1:dEll:end)),1);
    zr12= reshape(fmaj1(1:dEll:end,1:dEll:end,indPlot),numel(glon1(1:dEll:end,1:dEll:end)),1);
    zr22= reshape(fmin1(1:dEll:end,1:dEll:end,indPlot),numel(glon1(1:dEll:end,1:dEll:end)),1);
    zr32= reshape(finc1(1:dEll:end,1:dEll:end,indPlot),numel(glon1(1:dEll:end,1:dEll:end)),1);
    
    for iii=1:length(zr12)
        [ll1,ll2] = m_ellipse(xre2(iii),yre2(iii),90-zr32(iii),zr12(iii),zr22(iii),100,ellScl);

        a = m_plot(ll1,ll2,'Color',[234 95 148]/256,'LineWidth',2);
        a.Color(4) = [0.65];
    end

    % For plotting scale ellipse 
    [ll1,ll2] = m_ellipse(152.73,-33.44,90-0,0.25,0.1,100,ellScl);
    m_plot(ll1,ll2,'Color',[234 95 148]/256,'LineWidth',2.0);
    m_text(152.55,-33.52,'0.25 m^2 s^-^2');
    
    % m_plot(xy(:,1),xy(:,2),'k','LineWidth',1.3);
    
    %
    m_grid('GridLineStyle','none');clim([0 1.2]);colorbar;
    titStr = sprintf('%s',tmonthStrs(indPlot,:));
    % title( titStr );
    m_text(151.67,-32.45,titStr,'FontWeight','bold','FontSize',17);
    
    % Plot radar sites
    m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
        'gs','MarkerSize',11,'MarkerFaceColor','g');
    m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
        'gs','MarkerSize',11,'MarkerFaceColor','g');
    m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
        'rs','MarkerSize',11,'MarkerFaceColor','r');
    m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
        'rs','MarkerSize',11,'MarkerFaceColor','r');
end

% Set title 
ax=nexttile(1);title('Austral summer','FontSize',17);
ax=nexttile(2);title('Austral autumn','FontSize',17);
ax=nexttile(3);title('Austral winter','FontSize',17);
ax=nexttile(4);title('Austral spring','FontSize',17);
%% For plotting monthly COF
m_proj('lambert','lon',LonLims2,'lat',LatLims2);

vecScl = 10.5;
vecHeadangle=50;
vecShaftwidth=0.025;
vecHeadwidth=0.25;
vecHeadlength=0.33;
dx = 4;

% For Ellipses
dEll = 9;

ellScl = 40000;

% For plotting COF
h=figure(3);clf;hold on;
set(h, 'Position',[316 182 1090 659]);
tiledlayout(3,4,"TileSpacing","compact");

% Plot radar vel 
for indPlot=1:12
    
    switch indPlot
        % Austral summer 
        case 12
            figPos = 1;
        case 1 
            figPos = 5;
        case 2
            figPos = 9;

        % Austral autumn 
        case 3 
            figPos = 2;
        case 4 
            figPos = 6;
        case 5 
            figPos = 10;

        % Austral winter 
        case 6 
            figPos = 3;
        case 7 
            figPos = 7;
        case 8 
            figPos = 11;

        % Austral spring 
        case 9 
            figPos = 4;
        case 10 
            figPos = 8;
        case 11 
            figPos = 12;

    end
    ax=nexttile(figPos);hold on;
    % Bathy lines 
    [C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
        'ShowText','off','Color',[.6 .6 .6]);
    clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier');
    % Coastline
    m_plot(clon,clat,'k','Linewidth',1.2);

    % Vectors
    xr = reshape(glon2(1:dx:end,1:dx:end),numel(glon2(1:dx:end,1:dx:end)),1);
    yr = reshape(glat2(1:dx:end,1:dx:end),numel(glat2(1:dx:end,1:dx:end)),1);

    uPlot = reshape(Umonthly2(1:dx:end,1:dx:end,indPlot),numel(glon2(1:dx:end,1:dx:end)),1);
    vPlot = reshape(Vmonthly2(1:dx:end,1:dx:end,indPlot),numel(glon2(1:dx:end,1:dx:end)),1);
    colormap(cm_delta(32));
    m_vec( vecScl*2,xr,yr,uPlot,vPlot,sqrt( uPlot.^2 + vPlot.^2 ),...
                'edgeclip', 'k', 'headangle',vecHeadangle,...
                'shaftwidth',vecShaftwidth,'headwidth',vecHeadwidth,...
                'headlength',vecHeadlength );
    
    % Ellipses
    xre2 = reshape(glon2(1:dEll:end,1:dEll:end),numel(glon2(1:dEll:end,1:dEll:end)),1);
    yre2 = reshape(glat2(1:dEll:end,1:dEll:end),numel(glon2(1:dEll:end,1:dEll:end)),1);
    zr12= reshape(fmaj2(1:dEll:end,1:dEll:end,indPlot),numel(glon2(1:dEll:end,1:dEll:end)),1);
    zr22= reshape(fmin2(1:dEll:end,1:dEll:end,indPlot),numel(glon2(1:dEll:end,1:dEll:end)),1);
    zr32= reshape(finc2(1:dEll:end,1:dEll:end,indPlot),numel(glon2(1:dEll:end,1:dEll:end)),1);
    
    for iii=1:length(zr12)
        [ll1,ll2] = m_ellipse(xre2(iii),yre2(iii),90-zr32(iii),zr12(iii),zr22(iii),100,ellScl);

        a = m_plot(ll1,ll2,'Color',[234 95 148]/256,'LineWidth',1.5);
        a.Color(4) = [0.65];
    end

    % For plotting scale ellipse 
    [ll1,ll2] = m_ellipse(153.90,-31.0,90-0,0.25,0.1,100,ellScl);
    m_plot(ll1,ll2,'Color',[234 95 148]/256,'LineWidth',2.0);
    m_text(153.68,-31.08,'0.25 m^2 s^-^2');
%
    % m_plot(xy(:,1),xy(:,2),'k','LineWidth',1.3);
    
    %
    m_grid('GridLineStyle','none');clim([0 1.2]);colorbar;
    titStr = sprintf('%s',tmonthStrs(indPlot,:));
    % title( titStr );
    m_text(152.80,-30.0,titStr,'FontWeight','bold','FontSize',17);
    
    % Plot radar sites
    m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
        'gs','MarkerSize',11,'MarkerFaceColor','g');
    m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
        'gs','MarkerSize',11,'MarkerFaceColor','g');
    m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
        'rs','MarkerSize',11,'MarkerFaceColor','r');
    m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
        'rs','MarkerSize',11,'MarkerFaceColor','r');
end

% Set title 
ax=nexttile(1);title('Austral summer','FontSize',17);
ax=nexttile(2);title('Austral autumn','FontSize',17);
ax=nexttile(3);title('Austral winter','FontSize',17);
ax=nexttile(4);title('Austral spring','FontSize',17);
%% Now computing domain averaged current velocity
u1mean = squeeze( nanmean(nanmean(u1.*maskNEWC,2),1)) ;
v1mean = squeeze( nanmean(nanmean(v1.*maskNEWC,2),1));
u2mean = squeeze( nanmean(nanmean(u2.*maskCOF,2),1) );
v2mean = squeeze( nanmean(nanmean(v2.*maskCOF,2),1) );
%% Rotate vectors 
rot_angle1 = 30;
rot_angle2 = 18;

rotmat1 = [cosd(rot_angle1) -sind(rot_angle1);sind(rot_angle1) cosd(rot_angle1)];
rotmat2 = [cosd(rot_angle2) -sind(rot_angle2);sind(rot_angle2) cosd(rot_angle2)];

vel1Rot = rotmat1 * [u1mean';v1mean'];
vel2Rot = rotmat2 * [u2mean';v2mean'];

u1meanRot = vel1Rot(1,:);v1meanRot = vel1Rot(2,:);
u2meanRot = vel2Rot(1,:);v2meanRot = vel2Rot(2,:);

%% Sort intra-annual data 
nb1 = nan(length(taxMonthly),1);
nb2 = nan(length(taxMonthly),1);
u1Monthly = nan(length(taxMonthly),1);v1Monthly = nan(length(taxMonthly),1);
u2Monthly = nan(length(taxMonthly),1);v2Monthly = nan(length(taxMonthly),1);
u1MonthlyStd = nan(length(taxMonthly),1);v1MonthlyStd = nan(length(taxMonthly),1);
u2MonthlyStd = nan(length(taxMonthly),1);v2MonthlyStd = nan(length(taxMonthly),1);

for ii=1:length(taxMonthly)
    % For NEWC
    ind = find(month(taxMonthly(ii)) == month(t1));
    nb1(ii) = length( find(~isnan(u1meanRot(ind))) );

    u1Monthly(ii) = nanmean( u1meanRot(ind) );
    u1MonthlyStd(ii) = nanstd( u1meanRot(ind) );
    v1Monthly(ii) = nanmean( v1meanRot(ind) );
    v1MonthlyStd(ii) = nanstd( v1meanRot(ind) );

    % For COF
    ind = find(month(taxMonthly(ii)) == month(t2));
    nb2(ii) = length( find(~isnan(u2meanRot(ind))) );
    u2Monthly(ii) = nanmean( u2meanRot(ind) );
    u2MonthlyStd(ii) = nanstd( u2meanRot(ind) );
    v2Monthly(ii) = nanmean( v2meanRot(ind) );
    v2MonthlyStd(ii) = nanstd( v2meanRot(ind) );
end

%% Plot domain averaged data 
colmap = cm_balance(10);
% 
h=figure(4);clf;hold on;
set(h,'Position', [480 432 745 200]);
tiledlayout(1,2,"TileSpacing",'compact');

% Cross-shore
nexttile;hold on;
plot([tax(1)-30 tax(end)],[0 0],'--','Color',[.8 .8 .8],'LineWidth',0.4);
a = patch([taxMonthly;flipud(taxMonthly)],[u1Monthly+u1MonthlyStd;flipud(u1Monthly-u1MonthlyStd)], ...
    colmap(3,:),'EdgeColor','none');
set(a,'FaceAlpha',0.2);

a = patch([taxMonthly;flipud(taxMonthly)],[u2Monthly+u2MonthlyStd;flipud(u2Monthly-u2MonthlyStd)], ...
    colmap(8,:),'EdgeColor','none');
set(a,'FaceAlpha',0.2);

plot(taxMonthly,u1Monthly,'-o','LineWidth',2,'Color',colmap(3,:));
plot(taxMonthly,u2Monthly,'-o','LineWidth',2,'Color',colmap(8,:));
xticks(taxMonthly);xticklabels(datestr(taxMonthly,'mmm'));
xlim([tax(1)-30 tax(end)]);ylim([-0.3 0.3]);ylabel('Vel. (m s^-^1)');
text(tax(1)-15,0.25,'Cross-shore','FontWeight','bold','FontSize',12);
text(tax(340),-0.20,'a','FontWeight','bold','FontSize',12);
legend({'','','','NEWC','COF'},'Orientation','horizontal','Box','off');

% Along-shore
nexttile;hold on;
plot([tax(1)-30 tax(end)],[0 0],'--','Color',[.8 .8 .8],'LineWidth',0.4);
a = patch([taxMonthly;flipud(taxMonthly)],[v1Monthly+v1MonthlyStd;flipud(v1Monthly-v1MonthlyStd)], ...
    colmap(3,:),'EdgeColor','none');
set(a,'FaceAlpha',0.2);

a = patch([taxMonthly;flipud(taxMonthly)],[v2Monthly+v2MonthlyStd;flipud(v2Monthly-v2MonthlyStd)], ...
    colmap(8,:),'EdgeColor','none');
set(a,'FaceAlpha',0.2);

plot(taxMonthly,v1Monthly,'-o','LineWidth',2,'Color',colmap(3,:));
plot(taxMonthly,v2Monthly,'-o','LineWidth',2,'Color',colmap(8,:));
xticks(taxMonthly);xticklabels(datestr(taxMonthly,'mmm'));
xlim([tax(1)-30 tax(end)]);ylim([-1.2 0.5]);ylabel('Vel. (m s^-^1)');
text(tax(1)-15,0.30,'Along-shore','FontWeight','bold','FontSize',12);
text(tax(340),-1.0,'b','FontWeight','bold','FontSize',12);
legend({'','','','NEWC','COF'},'Orientation','horizontal','Box','off');


%% Lastly calculate multi-year averaged radar vel
disp(['Computing NEWC mean vel ...']);

fmaj1total = nan(size(glon1,1),size(glon1,2));
fmin1total = nan(size(glon1,1),size(glon1,2));
finc1total = nan(size(glon1,1),size(glon1,2));

fmaj2total = nan(size(glon2,1),size(glon2,2));
fmin2total = nan(size(glon2,1),size(glon2,2));
finc2total = nan(size(glon2,1),size(glon2,2));

% For NEWC 
ind = find(iDatAvail1' == 1);
UmeanTotal1 = nanmean(u1(:,:,ind),3).*maskNEWC;
VmeanTotal1 = nanmean(v1(:,:,ind),3).*maskNEWC;

SpdmeanTotal1 = sqrt( UmeanTotal1.^2 + VmeanTotal1.^2 );

% Compute var ellipses for NEWC
[nx,ny] = size(glon1);
for xind=1:nx
    for yind=1:ny
        if isnan(UmeanTotal1(xind,yind))
            continue
        end
        % 
        [fmaj1total(xind,yind),fmin1total(xind,yind),finc1total(xind,yind)] = ...
            get_PCA( squeeze(u1(xind,yind,ind)), squeeze(v1(xind,yind,ind)) );
    end
end

disp(['Computing COF mean vel ...']);

% For COF 
ind = find(iDatAvail2' == 1);
UmeanTotal2 = nanmean(u2(:,:,ind),3).*maskCOF;
VmeanTotal2 = nanmean(v2(:,:,ind),3).*maskCOF;

SpdmeanTotal2 = sqrt( UmeanTotal2.^2 + VmeanTotal2.^2 );

% Compute var ellipses for COF
[nx,ny] = size(glon2);
for xind=1:nx
    for yind=1:ny
        if isnan(UmeanTotal2(xind,yind))
            continue
        end
        % 
        [fmaj2total(xind,yind),fmin2total(xind,yind),finc2total(xind,yind)] = ...
            get_PCA( squeeze(u2(xind,yind,ind)), squeeze(v2(xind,yind,ind)) );
    end
end
%% For plotting mean current vel
m_proj('lambert','lon',LonLims,'lat',LatLims);

vecScl = 2.5;
vecHeadangle=60;
vecShaftwidth=0.055;
vecHeadwidth=0.55;
vecHeadlength=0.73;
ellScl = 35000;

% For plotting COF
h=figure(5);clf;hold on;
set(h, 'Position',[661 155 467 511]);

% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier');
% Coastline
m_plot(clon,clat,'k','Linewidth',1.2);

% For NEWC
dx = 1;
xr = reshape(glon1(1:dx:end,1:dx:end),numel(glon1(1:dx:end,1:dx:end)),1);
yr = reshape(glat1(1:dx:end,1:dx:end),numel(glat1(1:dx:end,1:dx:end)),1);

uPlot = reshape(UmeanTotal1(1:dx:end,1:dx:end),numel(glon1(1:dx:end,1:dx:end)),1);
vPlot = reshape(VmeanTotal1(1:dx:end,1:dx:end),numel(glon1(1:dx:end,1:dx:end)),1);
colormap(cm_delta(32));
m_vec( vecScl*2,xr,yr,uPlot,vPlot,sqrt( uPlot.^2 + vPlot.^2 ),...
            'edgeclip', 'k', 'headangle',vecHeadangle,...
            'shaftwidth',vecShaftwidth,'headwidth',vecHeadwidth,...
            'headlength',vecHeadlength );

% Ellipses
% NEWC
dEll = 2;
xre2 = reshape(glon1(1:dEll:end,1:dEll:end),numel(glon1(1:dEll:end,1:dEll:end)),1);
yre2 = reshape(glat1(1:dEll:end,1:dEll:end),numel(glon1(1:dEll:end,1:dEll:end)),1);
zr12= reshape(fmaj1total(1:dEll:end,1:dEll:end),numel(glon1(1:dEll:end,1:dEll:end)),1);
zr22= reshape(fmin1total(1:dEll:end,1:dEll:end),numel(glon1(1:dEll:end,1:dEll:end)),1);
zr32= reshape(finc1total(1:dEll:end,1:dEll:end),numel(glon1(1:dEll:end,1:dEll:end)),1);

for iii=1:length(zr12)
    [ll1,ll2] = m_ellipse(xre2(iii),yre2(iii),90-zr32(iii),...
        zr12(iii),zr22(iii),100,ellScl);

    a = m_plot(ll1,ll2,'Color',[0.87,0.53,0.29],'LineWidth',1.3);
    a.Color(4) = [0.65];
end

% For COF
dx = 3;
% Vectors
xr = reshape(glon2(1:dx:end,1:dx:end),numel(glon2(1:dx:end,1:dx:end)),1);
yr = reshape(glat2(1:dx:end,1:dx:end),numel(glat2(1:dx:end,1:dx:end)),1);

uPlot = reshape(UmeanTotal2(1:dx:end,1:dx:end),numel(glon2(1:dx:end,1:dx:end)),1);
vPlot = reshape(VmeanTotal2(1:dx:end,1:dx:end),numel(glon2(1:dx:end,1:dx:end)),1);
colormap(cm_delta(32));
m_vec( vecScl*2,xr,yr,uPlot,vPlot,sqrt( uPlot.^2 + vPlot.^2 ),...
            'edgeclip', 'k', 'headangle',vecHeadangle,...
            'shaftwidth',vecShaftwidth,'headwidth',vecHeadwidth,...
            'headlength',vecHeadlength );

% Ellipses
% COF
dEll = 8;
xre2 = reshape(glon2(1:dEll:end,1:dEll:end),numel(glon2(1:dEll:end,1:dEll:end)),1);
yre2 = reshape(glat2(1:dEll:end,1:dEll:end),numel(glon2(1:dEll:end,1:dEll:end)),1);
zr12= reshape(fmaj2total(1:dEll:end,1:dEll:end),numel(glon2(1:dEll:end,1:dEll:end)),1);
zr22= reshape(fmin2total(1:dEll:end,1:dEll:end),numel(glon2(1:dEll:end,1:dEll:end)),1);
zr32= reshape(finc2total(1:dEll:end,1:dEll:end),numel(glon2(1:dEll:end,1:dEll:end)),1);

for iii=1:length(zr12)
    [ll1,ll2] = m_ellipse(xre2(iii),yre2(iii),90-zr32(iii), ...
        zr12(iii),zr22(iii),100,ellScl);

    a = m_plot(ll1,ll2,'Color',[0.87,0.53,0.29],'LineWidth',1.3);
    a.Color(4) = [0.65];
end

% For plotting scale ellipse 
[ll1,ll2] = m_ellipse(150.90,-31.0,90-0,0.25,0.1,100,ellScl);
m_plot(ll1,ll2,'Color',[0.87,0.53,0.29],'LineWidth',2.0);
m_text(150.68,-31.18,'0.25 m^2 s^-^2');
%
% m_plot(xy(:,1),xy(:,2),'k','LineWidth',1.3);

%
m_grid('GridLineStyle','none');clim([0 1.2]);colorbar;

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

%% Location 
ll1 = [-33.856765, 151.229640]; % Sydney
ll2 = [-32.4132,   152.3170  ]; % Newcastle radar site 
ll3 = [-32.706386, 152.125112]; % Nelson Bay
ll4 = [-34.429123, 150.909155]; % Wollongong
ll5 = [-31.467922, 152.901629]; % Port Macquarie
ll6 = [-31.879287, 152.686935]; % Harrington
ll7 = [-30.655346, 153.046770]; % Nambucca Heads
ll8 = [-30.288152, 153.143107]; % Coffs Harbour radar site
m_text(ll1(2)-0.68,ll1(1),sprintf('Sydney'),'Color','k', ...
    'FontWeight','bold','FontSize',12);
m_text(ll2(2)-1.08,ll2(1),sprintf('Newcastle \n  (NEWC)'),'Color', ...
    'k','FontWeight','bold','FontSize',12);
m_text(ll3(2)-0.90,ll3(1),sprintf('Nelson Bay'),'Color','k');
% m_text(ll4(2)-0.75,ll4(1),sprintf('Wollongong'),'Color','k');
m_text(ll5(2)-0.95,ll5(1),sprintf('Port Macquarie'),'Color','k');
% m_text(ll6(2)-0.70,ll6(1),sprintf('Harrington'),'Color','k');
% m_text(ll7(2)-1.10,ll7(1),sprintf('Nambucca Heads'),'Color','k');
m_text(ll8(2)-0.85,ll8(1),sprintf('Coffs Habour \n       (COF)'), ...
    'Color','k','FontWeight','bold');

%% Function 
% PCA the vel 
function [fmaj,fmin,finc]=get_PCA(x,y)

    fmaj=NaN;fmin=NaN;finc=NaN;
    c=cov(x,y,'omitrows');
    [v,d] = svd(c);
    
    fmaj = d(1,1);
    fmin = d(2,2);
    finc = atan2d(v(2,1),v(1,1));
    
    return
end