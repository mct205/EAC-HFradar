%% Require to run the script get_tidal_analysis_Utide.m first!! 
clear;
%% Load data
load(['../DATA/processed-data/' ...
    'Tidal_utide_20200101_20210101.mat']);
%% 
tidalMethod = 'utide';
% tidalMethod = 't_tide';
conLists = ['K1';'M2';'O1';'S2'];

[nx1,ny1] = size(glon1);[nx2,ny2] = size(glon2);
%% tt1
% K1
[k1_fmaj_t1,k1_fmin_t1,k1_fpha_t1,k1_finc_t1] = ...
    get_tidal_val(tt1,conLists(1,:),tidalMethod,nx1,ny1);
% M2
[m2_fmaj_t1,m2_fmin_t1,m2_fpha_t1,m2_finc_t1] = ...
    get_tidal_val(tt1,conLists(2,:),tidalMethod,nx1,ny1);
% O1
[o1_fmaj_t1,o1_fmin_t1,o1_fpha_t1,o1_finc_t1] = ...
    get_tidal_val(tt1,conLists(3,:),tidalMethod,nx1,ny1);
% S2
[s2_fmaj_t1,s2_fmin_t1,s2_fpha_t1,s2_finc_t1] = ...
    get_tidal_val(tt1,conLists(4,:),tidalMethod,nx1,ny1);
%% tt2
% K1
[k1_fmaj_t2,k1_fmin_t2,k1_fpha_t2,k1_finc_t2] = ...
    get_tidal_val(tt2,conLists(1,:),tidalMethod,nx2,ny2);
% M2
[m2_fmaj_t2,m2_fmin_t2,m2_fpha_t2,m2_finc_t2] = ...
    get_tidal_val(tt2,conLists(2,:),tidalMethod,nx2,ny2);
% O1
[o1_fmaj_t2,o1_fmin_t2,o1_fpha_t2,o1_finc_t2] = ...
    get_tidal_val(tt2,conLists(3,:),tidalMethod,nx2,ny2);
% S2
[s2_fmaj_t2,s2_fmin_t2,s2_fpha_t2,s2_finc_t2] = ...
    get_tidal_val(tt2,conLists(4,:),tidalMethod,nx2,ny2);
%
xe1 = [152.1508, 152.3077, 152.4424];
ye1 = [-32.9132, -33.0556, -33.1882];
%
xe2 = [153.3282, 153.5728, 153.7828];
ye2 = [-30.2354, -30.2858, -30.3711];
%%
load('../DATA/aux_data/EAC_coastline.mat');
load('../DATA/aux_data/EAC_bathy.mat');
% load radar site locations 
load ../DATA/aux_data/EAC_radarstations.mat

ch70_pos = [153.30 -30.2758];
ch100_pos= [153.3978 -30.2688];
%% 

%% Domain
LonLims = [151.5 154.0];
LatLims = [-33.5 -29.8];
% NEWC
LonLims1 = [151.5 153.5];
LatLims1 = [-34.0 -32.0];
% COF
LonLims2 = [152.9 154.2];
LatLims2 = [-30.8 -29.8];
%%
m_proj('lambert','lon',LonLims,'lat',LatLims);

h=figure(1);clf;hold on;

% Plot K1
ax=subplot(1,2,1);hold on;
set(h, 'Position', [455 167 748 587]);

% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,k1_fmaj_t1);
m_pcolor(glon2,glat2,k1_fmaj_t2);
colormap(ax,cm_speed(32));clim([0 0.1]);colorbar;

% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
%
title( [ 'Tidal analysis on HF-derived surface velocity, tidal const ' conLists(1,:)] );

% Plot M2
ax=subplot(1,2,2);hold on;

% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,m2_fmaj_t1);
m_pcolor(glon2,glat2,m2_fmaj_t2);
colormap(ax,cm_curl(32));clim([0 0.1]);colorbar;

% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
%
title( [ 'Tidal analysis on HF-derived surface velocity, tidal const ' conLists(2,:)] );
%% Plot for each site 
% Plot NEWC region 
m_proj('lambert','lon',LonLims1,'lat',LatLims1);

h=figure(2);clf;
set(h,'Position', [563 363 946 456]);

tiledlayout(2,4,"TileSpacing",'tight');

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,k1_fmaj_t1);colorbar
colormap(ax,cm_speed(32));clim([0 0.1]);
% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(1,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,m2_fmaj_t1);colorbar;
colormap(ax,cm_speed(32));clim([0 0.1]);
% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(2,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,o1_fmaj_t1);colorbar;
colormap(ax,cm_speed(32));clim([0 0.1]);
% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(3,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,s2_fmaj_t1);colorbar;
colormap(ax,cm_speed(32));clim([0 0.1]);
% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(4,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,k1_fpha_t1);colorbar;
colormap(ax,cm_curl(32));clim([0 360]);
% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

title([conLists(1,:)]);
m_grid('GridLineStyle','none');

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,m2_fpha_t1);colorbar;
colormap(ax,cm_curl(32));clim([0 360]);
% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(2,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,o1_fpha_t1);colorbar;
colormap(ax,cm_curl(32));clim([0 360]);
% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(3,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon1,glat1,s2_fpha_t1);colorbar;
colormap(ax,cm_curl(32));clim([0 360]);
% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(4,:)]);

%% For COF site 
m_proj('lambert','lon',LonLims2,'lat',LatLims2);

h=figure(3);clf;hold on;
set(h,'Position', [563 363 946 456]);
tiledlayout(2,4,"TileSpacing",'compact');

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon2,glat2,k1_fmaj_t2);colorbar
colormap(ax,cm_speed(32));clim([0 0.1]);
% Plot mooring site 
m_plot(ch70_pos(1),ch70_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');
m_plot(ch100_pos(1),ch100_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');

% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(1,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon2,glat2,m2_fmaj_t2);colorbar;
colormap(ax,cm_speed(32));clim([0 0.1]);
% Plot mooring site 
m_plot(ch70_pos(1),ch70_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');
m_plot(ch100_pos(1),ch100_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');

% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(2,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon2,glat2,o1_fmaj_t2);colorbar;
colormap(ax,cm_speed(32));clim([0 0.1]);
% Plot mooring site 
m_plot(ch70_pos(1),ch70_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');
m_plot(ch100_pos(1),ch100_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');

% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(3,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon2,glat2,s2_fmaj_t2);colorbar;
colormap(ax,cm_speed(32));clim([0 0.1]);
% Plot mooring site 
m_plot(ch70_pos(1),ch70_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');
m_plot(ch100_pos(1),ch100_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');

% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(4,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon2,glat2,k1_fpha_t2);colorbar;
colormap(ax,cm_curl(32));clim([0 360]);
% Plot mooring site 
m_plot(ch70_pos(1),ch70_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');
m_plot(ch100_pos(1),ch100_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');

title([conLists(1,:)]);
% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon2,glat2,m2_fpha_t2);colorbar;
colormap(ax,cm_curl(32));clim([0 360]);
% Plot mooring site 
m_plot(ch70_pos(1),ch70_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');
m_plot(ch100_pos(1),ch100_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');

% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(2,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon2,glat2,o1_fpha_t2);colorbar;
colormap(ax,cm_curl(32));clim([0 360]);
% Plot mooring site 
m_plot(ch70_pos(1),ch70_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');
m_plot(ch100_pos(1),ch100_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');

% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(3,:)]);

ax=nexttile;hold on;
% Bathy lines 
[C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
    'ShowText','off','Color',[.6 .6 .6]);
clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier')

m_pcolor(glon2,glat2,s2_fpha_t2);colorbar;
colormap(ax,cm_curl(32));clim([0 360]);
% Plot mooring site 
m_plot(ch70_pos(1),ch70_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');
m_plot(ch100_pos(1),ch100_pos(2),'^','MarkerFaceColor', ...
    'b','MarkerEdgeColor','k');

% Plot coast line 
m_plot(clon,clat,'k','Linewidth',1.2);

% Plot radar sites
m_plot(sites(1).SiteOrigin(2),sites(1).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(2).SiteOrigin(2),sites(2).SiteOrigin(1), ...
    'gs','MarkerSize',11,'MarkerFaceColor','g');
m_plot(sites(3).SiteOrigin(2),sites(3).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');
m_plot(sites(4).SiteOrigin(2),sites(4).SiteOrigin(1), ...
    'rs','MarkerSize',11,'MarkerFaceColor','r');

m_grid('GridLineStyle','none');
title([conLists(4,:)]);

%%
function [fmaj,fmin,fpha,finc] = get_tidal_val(tt,const,tidalMethod,nx,ny)

fmaj = nan(nx,ny);
fmin = nan(nx,ny);
fpha = nan(nx,ny);
finc = nan(nx,ny);

switch tidalMethod
    case 't_tide'
        for jj=1:length(tt.nameu)
            if strcmpi( strtrim(tt.nameu(jj,:)),const )
                ind = jj;
                break;
            end
        end
        %
        fmaj = tt.fmaj(:,:,ind);
        fmin = tt.fmin(:,:,ind);
        fpha = tt.fpha(:,:,ind);
        finc = tt.finc(:,:,ind);
    case 'utide'
        for ix=1:nx
            for iy=1:ny
                if isempty(tt.nameu{ix,iy})
                    continue
                end
                indC = find(ismember(tt.nameu{ix,iy},const));
                if ~isempty(indC)
                    ind = indC;
                else
                    ind = NaN;
                end
                %
                fmaj(ix,iy) = tt.fmaj(ix,iy,ind);
                fmin(ix,iy) = tt.fmin(ix,iy,ind);
                fpha(ix,iy) = tt.fpha(ix,iy,ind);
                finc(ix,iy) = tt.finc(ix,iy,ind);
            end
        end
    otherwise
        disp('Nothing there!');
end

return
end