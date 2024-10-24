clear;
%% Load total grid 
load('../DATA/aux_data/COF_grid.mat');

[nx,ny] = size(glon);
refPoint = [152.94,-31.5];
[x0,y0] = deal(refPoint(1),refPoint(2));
%% Data info 
siteNames = {'NNB','RRK'};
rPrefix = 'IMOS_ACORN_RV';
rType = 'wera';
%% Load sample data
load('sample_data/Radial_COF_sample.mat')
%% Time to load data
t = datenum(2015,01,01,03,00,00):1/24:datenum(2015,01,01,03,00,00);

t1(1,:) = [t(1)-30/60/24:10/60/24:t(end)+30/60/24];
t1(2,:) = [t(1)-30/60/24:10/60/24:t(end)+30/60/24]+5/60/24; % Since NNB time 
                                                            % is offset by 
                                                            % 5-min compared 
                                                            % to RRK
%% Data out infomation 
make_plot_flag = 1;
indPlot = 1;

%% If plot, load some aux data
if make_plot_flag
    load('../DATA/aux_data/EAC_coastline.mat');
    load('../DATA/aux_data/EAC_bathy.mat');
    load('../DATA/aux_data/EAC_radarstations.mat');
end

%% Make 2dVar gap-filling surface current maps
U = nan(nx,ny,length(t));
V = nan(nx,ny,length(t));

cDir = pwd;
pDir = [pwd,'/2dVar-COF'];
cd(pDir);
fOut = [pDir,'/out12.dat'];

for indt = 1:length(t)

    nsites = 0;

    xr1=[];yr1=[];x1=[];y1=[];sn1=[];cs1=[];rcomp1=[];

    for nsite=1:length(siteNames)
        % Perform moving average to have the coverage of data time
        % Time index should be like this: NNB: [05 15 25 35 45 55];
        %                                 RRK: [00 10 20 30 40 50];
        
        iind = find( abs( R(nsite).TimeStamp-t(indt) ) < 0.5/24 );

        % Reshape first
        if isempty( R(nsite).SiteOrigin )
            continue
        else
            llOrig = [R(nsite).SiteOrigin(1),R(nsite).SiteOrigin(2)];
        end

        % Then do the average
        lonTmp = mean( R(nsite).Lon(:,:,iind), 3 );
        latTmp = mean( R(nsite).Lat(:,:,iind), 3 );
        bearTmp  = mean( R(nsite).Bear(:,:,iind), 3 );
        rangeTmp = mean( R(nsite).Range(:,:,iind), 3 );
        rcompTmp = mean( R(nsite).RadComp(:,:,iind), 3 );
        %
        xr = reshape( lonTmp, numel(lonTmp), 1 );
        yr = reshape( latTmp, numel(latTmp), 1 );
        br = reshape( bearTmp, numel(bearTmp), 1 );
        rr = reshape( rangeTmp, numel(rangeTmp), 1 );
        rcomp = reshape( rcompTmp, numel(lonTmp), 1 );
        % 
        x = deg2km( llOrig(1)-x0 ) + deg2km( ( xr-llOrig(1) ) );
        y = deg2km( llOrig(2)-y0 ) + deg2km( ( yr-llOrig(2) ) );

        rbr = nan(length(xr),1);

        for ii=1:length(xr)
            uu = deg2km(xr(ii) - llOrig(1));
            vv = deg2km(yr(ii) - llOrig(2));
            rbr(ii) = mod(atan2d(vv,uu),360);
        end

        sn= sind( br );
        cs= cosd( br );

        % Remove NaN
        ind=find(isnan(rcomp));
        xr(ind)  = [];
        yr(ind)  = [];
        rbr(ind)  = [];
        rr(ind)  = [];
        x(ind)  = [];
        y(ind)  = [];
        sn(ind) = [];
        cs(ind) = [];
        rcomp(ind) = [];

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
    
    if  nsites < 2
        continue
    end
    
    % Write to file 
    if exist(fOut,'file')
        delete(fOut);
    end

    dlmwrite( fOut,[cs1 sn1 x1 y1 rcomp1],'delimiter',...
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
    tmp = dlmread([pDir '/u_opt.dat']);
    U(:,:,indt) = tmp';
    tmp = dlmread([pDir '/v_opt.dat']);
    V(:,:,indt) = tmp';
end
%
disp('Finish processing total current velocities!!');
cd(cDir);
%% Set NaN
U(U==0) = NaN;
V(V==0) = NaN;

%% For test plot
if make_plot_flag

    LonLims = [152.9 154.2];
    LatLims = [-31.2 -29.8];
    
    % Flag for plotting radar vectors
    vecScl = 54.5;
    vecHeadangle=60;
    vecShaftwidth=0.055;
    vecHeadwidth=0.25;
    vecHeadlength=0.55;
    dx = 2;
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