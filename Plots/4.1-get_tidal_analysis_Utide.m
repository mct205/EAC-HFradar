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
%
maskNEWC = nan(size(glon1)); % Masking grid with data
indnotnan = find( ~isnan(u1(:,:,40656)) );
maskNEWC(indnotnan) = 1;
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
%
maskCOF = nan(size(glon2)); % Masking grid with data
indnotnan = find( ~isnan(u2(:,:,75865)) );
maskCOF(indnotnan) = 1;

%% Cut data 
tStart = datenum(2020,01,01);
tEnd = datenum(2021,01,01);

% For NEWC
indT1 = find(t1 >= tStart & t1 <= tEnd);
u1 = u1(:,:,indT1).*maskNEWC; v1 = v1(:,:,indT1).*maskNEWC;
t1 = t1(indT1);

% For COF
indT2 = find(t2 >= tStart & t2 <= tEnd);
u2 = u2(:,:,indT2).*maskCOF; v2 = v2(:,:,indT2).*maskCOF;
t2 = t2(indT2);
%% 
[nx1,ny1] = size(glon1);
[nx2,ny2] = size(glon2);
%% Run test first with the cell contains data
uTmp = squeeze(u1(27,28,:));
vTmp = squeeze(v1(27,28,:));

c = ut_solv(t1,uTmp,vTmp,glat1(27,28),'auto','NodsatLinT','GwchLinT','OLS');

tt1.nameu=cell(nx1,ny1);
tt1.fu=nan(nx1,ny1,length(c.aux.frq));
tt1.fmaj = nan(nx1,ny1,length(c.aux.frq));
tt1.emaj = nan(nx1,ny1,length(c.aux.frq));
tt1.fmin = nan(nx1,ny1,length(c.aux.frq));
tt1.emin = nan(nx1,ny1,length(c.aux.frq));
tt1.finc = nan(nx1,ny1,length(c.aux.frq));
tt1.einc = nan(nx1,ny1,length(c.aux.frq));
tt1.fpha = nan(nx1,ny1,length(c.aux.frq));
tt1.epha = nan(nx1,ny1,length(c.aux.frq));
tt1.xout = nan(nx1,ny1,length(t1));
tt1.uvarRatio = nan(nx1,ny1);
tt1.vvarRatio = nan(nx1,ny1);
% 
uTmp = squeeze(u2(34,101,:));
vTmp = squeeze(v2(34,101,:));

c = ut_solv(t2,uTmp,vTmp,glat2(101,34),'auto','NodsatLinT','GwchLinT','OLS');

tt2.nameu=cell(nx2,ny2);
tt2.fu=nan(nx2,ny2,length(c.aux.frq));
tt2.fmaj = nan(nx2,ny2,length(c.aux.frq));
tt2.emaj = nan(nx2,ny2,length(c.aux.frq));
tt2.fmin = nan(nx2,ny2,length(c.aux.frq));
tt2.emin = nan(nx2,ny2,length(c.aux.frq));
tt2.finc = nan(nx2,ny2,length(c.aux.frq));
tt2.einc = nan(nx2,ny2,length(c.aux.frq));
tt2.fpha = nan(nx2,ny2,length(c.aux.frq));
tt2.epha = nan(nx2,ny2,length(c.aux.frq));
tt2.xout = nan(nx2,ny2,length(t2));
tt2.uvarRatio = nan(nx2,ny2);
tt2.vvarRatio = nan(nx2,ny2);
%% Perform tidal analysis
disp('Performing tidal analysis at site NEWC....');
% First radar site 
for xind=1:nx1
    for yind=1:ny1
        nd = sum( isfinite(u1(xind,yind,:)) );
        %
        if nd/length(t1) < 0.6
            continue
        end
        % 
        uTmp = squeeze( u1(xind,yind,:) );
        vTmp = squeeze( v1(xind,yind,:) );

        c = ut_solv(t1,uTmp,vTmp,glat1(xind,yind),'auto','NodsatLinT','GwchLinT','OLS');
        [uout,vout] = ut_reconstr(t1,c);

        % Put data into array
        tt1.nameu{xind,yind} = c.name;
        tt1.fu(xind,yind,:) = c.aux.frq;
        tt1.fmaj(xind,yind,:) = c.Lsmaj;
        tt1.emaj(xind,yind,:) = c.Lsmaj_ci;
        tt1.fmin(xind,yind,:) = c.Lsmin;
        tt1.emin(xind,yind,:) = c.Lsmin_ci;
        tt1.finc(xind,yind,:) = c.theta;
        tt1.einc(xind,yind,:) = c.theta_ci;
        tt1.fpha(xind,yind,:) = c.g;
        tt1.epha(xind,yind,:) = c.g_ci;
        tt1.xout(xind,yind,:) = complex(uout',vout');
        tt1.uvarRatio(xind,yind,:) = nanvar(uout)/nanvar(uTmp)*100;
        tt1.vvarRatio(xind,yind,:) = nanvar(vout)/nanvar(vTmp)*100;
    end
end

% Second radar site 
disp('Performing tidal analysis at site COF....');

for xind=1:nx2
    for yind=1:ny2
        nd = sum( isfinite(v2(xind,yind,:)) );
        %
        if nd/length(t2) < 0.7
            continue
        end
        % 
        uTmp = squeeze( v2(xind,yind,:) );
        vTmp = squeeze( v2(xind,yind,:) );
        c = ut_solv(t2,uTmp,vTmp,glat2(xind,yind),'auto','NodsatLinT','GwchLinT','OLS');
        [uout,vout] = ut_reconstr(t2,c);

        % Put data into array
        tt2.nameu{xind,yind} = c.name;
        tt2.fu(xind,yind,:) = c.aux.frq;
        tt2.fmaj(xind,yind,:) = c.Lsmaj;
        tt2.emaj(xind,yind,:) = c.Lsmaj_ci;
        tt2.fmin(xind,yind,:) = c.Lsmin;
        tt2.emin(xind,yind,:) = c.Lsmin_ci;
        tt2.finc(xind,yind,:) = c.theta;
        tt2.einc(xind,yind,:) = c.theta_ci;
        tt2.fpha(xind,yind,:) = c.g;
        tt2.epha(xind,yind,:) = c.g_ci;
        tt2.xout(xind,yind,:) = complex(uout',vout');
        tt2.uvarRatio(xind,yind,:) = nanvar(uout)/nanvar(uTmp)*100;
        tt2.vvarRatio(xind,yind,:) = nanvar(vout)/nanvar(vTmp)*100;
    end
end

%% Save to file 
disp('Pusing to file....');
save(['DATA/processed-data/Tidal_utide_' datestr(tStart,'yyyymmdd'),'_',datestr(tEnd,'yyyymmdd')],...
    't1','t2','tt1','tt2','u1','v2','v1','v2', ...
    'glon1','glat1','glon2','glat2','-v7.3');