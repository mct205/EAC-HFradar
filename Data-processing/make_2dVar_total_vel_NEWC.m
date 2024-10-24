clear;
%% Load total grid 
load('../DATA/aux_data/NEWC_grid.mat');

[nx,ny] = size(glon);
refPoint = [151,-34.50];
[x0,y0] = deal(refPoint(1),refPoint(2));

%% Data info 
siteNames = {'RHED','SEAL'};
rPrefix = 'IMOS_ACORN_RV';
rType = 'seasonde';
%% Load sample data
load('sample_data/Radial_NEWC_sample.mat')
%% Time to load data
t = datenum(2021,11,01,00,00,00):1/24:datenum(2021,11,01,00,00,00);

%% Data out infomation 
make_plot_flag = 1;
indPlot = 1;

%% If plot, load some aux data
if make_plot_flag
    load('../DATA/aux_data/EAC_coastline.mat');
    load('../DATA/aux_data/EAC_bathy.mat');
    load('../DATA/aux_data/EAC_radarstations.mat');
end
%% Mask file 
mask = load(['/Users/z3541616/Library/CloudStorage/OneDrive-UNSW/' ...
    'aus-data/2dVar/mask11.dat']);
mask = flipud(mask);
nbGrid = find( mask == 1 );

%% Calc radial stats
for nsite=1:length(siteNames)
    stat(nsite).nd = squeeze( sum( isfinite(R(nsite).RadComp), 1 ) );
end
%% 
critNdVal = 5;
%% Make 2dVar gap-filling surface current maps
U = nan(nx,ny,length(t));
V = nan(nx,ny,length(t));

cDir = pwd;
pDir = [pwd,'/2dVar-NEWC'];
cd(pDir);
fOut = [pDir,'/out12.dat'];

for indt = 1:length(t)
    % Do not continue if coverage of both radar is lower than critical
    % value 
    % if stat(1).nd(indt) + stat(2).nd(indt) < critNdVal
    %     continue;
    % end

    nsites = 0;

    xr1=[];yr1=[];x1=[];y1=[];sn1=[];cs1=[];rcomp1=[];

    for nsite=1:length(siteNames)
        % Reshape first
        if isempty( R(nsite).SiteOrigin )
            continue
        else
            llOrig = [R(nsite).SiteOrigin(1),R(nsite).SiteOrigin(2)];
        end
        %
        xr = R(nsite).Lon(:,indt);
        yr = R(nsite).Lat(:,indt);
        br = R(nsite).rBear(:,indt);
        rr = R(nsite).Range(:,indt);
        rcomp = R(nsite).RadComp(:,indt);
        % Calculate distance to the reference point
        % Since we work on Cartesian coordinates, assume the points are
        % lying on a horizontal flat plane
        uu0 = deg2km(llOrig(1)-x0);
        vv0 = deg2km(llOrig(2)-y0);
        ar0 = mod(atan2d(vv0,uu0),360);
        xr0 = deg2km( (xr-llOrig(1)) );
        yr0 = deg2km( (yr-llOrig(2)) );
        % Calculate distance of each radial points to the reference point
        x = uu0 + xr0;
        y = vv0 + yr0;
        %
        sn= sind( br );
        cs= cosd( br );
        % Remove NaN
        ind=find(isnan(rcomp));
        if ~isempty(ind)
            xr(ind)  = [];
            yr(ind)  = [];
            br(ind)  = [];
            rr(ind)  = [];
            x(ind)  = [];
            y(ind)  = [];
            sn(ind) = [];
            cs(ind) = [];
            rcomp(ind) = [];
        end
        if ~isempty(rcomp)
            nsites = nsites+1;
        end

        % Concate
        xr1= vertcat(xr1,xr);
        yr1= vertcat(yr1,yr);
        x1 = vertcat(x1,x);
        y1 = vertcat(y1,y);
        sn1= vertcat(sn1,sn);
        cs1= vertcat(cs1,cs);
        rcomp1= vertcat(rcomp1,rcomp);
    end 
    % Write to file 
    if exist(fOut,'file')
        delete(fOut);
    end

    if isempty(sn1) || nsites < 2
        continue
    end
        
    dlmwrite( fOut,[sn1 cs1 x1 y1 -1*rcomp1],'delimiter',...
    ' ','precision','%8.3f %8.3f %8.3f %8.3f %8.3f');

    % Now try with 2dVar 
    % Remove previous running to make sure we have the latest data
    if exist([pDir '/u_opt.dat'],'file')
        delete([pDir '/u_opt.dat']);
    end
    if exist([pDir '/v_opt.dat'],'file')
        delete([pDir '/v_opt.dat']);
    end
    % Then execute
    try
        eval(['!' pDir '/2dVar'])   %interpolation 2dvar
    catch ex
        disp(['Encounter error at ' datestr(t(indt))]);
        disp(ex)
        return
    end 

    % Finally read back the data
    U(:,:,indt) = dlmread([pDir '/u_opt.dat']);
    V(:,:,indt) = dlmread([pDir '/v_opt.dat']);
end
%
disp('Finish processing total current velocities!!');
cd(cDir);
%% Set NaN
U(U==0) = NaN;
V(V==0) = NaN;

%% For test plot
if make_plot_flag
    LonLims = [151 153.5];
    LatLims = [-34.2 -32];
    
    % Flag for plotting radar vectors
    vecScl = 10.5;
    vecHeadangle=60;
    vecShaftwidth=0.105;
    vecHeadwidth=0.5;
    vecHeadlength=1.05;
    dx = 1;
    %
    h=figure(1);clf;hold on;
    %
    m_proj('lambert','lon',LonLims,'lat',LatLims);
    
    % Bathy lines 
    [C,h] = m_contour(xbath,ybath,zbath,[-100 -200 -1000 -2000 -3000 -4000],...
        'ShowText','off','Color',[.6 .6 .6]);
    clabel(C,h,'FontSize',8,'Color',[.6 .6 .6],'FontName','Courier');
    % Coastline
    m_plot(clon,clat,'k','Linewidth',1.2);
    
        xr = reshape(glon(1:dx:end,1:dx:end),numel(glon(1:dx:end,1:dx:end)),1);
    yr = reshape(glat(1:dx:end,1:dx:end),numel(glat(1:dx:end,1:dx:end)),1);
    
    uPlot = reshape(U(1:dx:end,1:dx:end,indPlot),numel(glon(1:dx:end,1:dx:end)),1);
    vPlot = reshape(V(1:dx:end,1:dx:end,indPlot),numel(glon(1:dx:end,1:dx:end)),1);
    
    colormap(cm_delta(32));
    m_vec( vecScl,xr,yr,uPlot,vPlot,sqrt( uPlot.^2 + vPlot.^2 ),...
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
    
    m_grid('GridLineStyle','none');clim([0 1]);colorbar;
    
    m_text(151.2735,-32.1794,datestr(t(indPlot),'yyyy-mm-ddTHH:MM:SSZ'), ...
        'FontWeight','bold','FontSize',12);
end