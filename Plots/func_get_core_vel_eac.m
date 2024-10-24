function [lonCoreVel,latCoreVel,lonCoreWidth,latCoreWidth,coreVel,coreAngle,coreWidth,jetxAx,jetAxVel]=func_get_core_vel_eac(glon,glat,U,V,dRx,ddist_from_coastVal,nbJetAx,searchRadius)
%%% This function tries to pop out the core velocity of the eac by shifting 
%%% the geogrpahical coordinate to jet cooordinate following the
%%% method of Archer et al., 2017: doi: 10.1002/2017JC013286

%% Now start identify the core of the jet
lonCoreVel = nan(size(glat,1),1);
latCoreVel = nan(size(glat,1),1);
lonCoreWidth = nan(size(glat,1),2,1);
latCoreWidth = nan(size(glat,1),2,1);

coreVel = nan(size(glat,1),1);
coreAngle = nan(size(glat,1),1);
coreWidth = nan(size(glat,1),1);

jetxAx = [-1*nbJetAx*dRx:dRx:nbJetAx*dRx];
jetxAxInd = round( jetxAx/dRx );
jetAxVel = nan(size(glat,1),length(jetxAx));

%% In search for starting point 
spd = sqrt( U.^2 + V.^2 );
[iFirst,jFirst] = find(spd == nanmax(nanmax(spd)));

while length( find(~isnan(spd(iFirst,:))) ) < 3
    spd( spd == nanmax(nanmax(spd)) ) = NaN;
    [iFirst,jFirst] = find(spd == nanmax(nanmax(spd)));
    break
end

if nanmax(nanmax(spd)) < 0.8 || isnan(nanmax(nanmax(spd)))
    disp('Eac does not pass!! Vel. less than 0.8 m/s');
    return
end

ddist_from_coast = deg2km(sqrt( (glon-glon(:,1)).^2 + (glat-glat(:,1)).^2 ));

latInd = iFirst; % Set first point

% Identify the first core vel
USlice = squeeze(U(latInd,:));
VSlice = squeeze(V(latInd,:));

spdSlice = sqrt( USlice.^2 + VSlice.^2 );
% Interpolate the gaps 
indnan = find(isnan(spdSlice));
indnotnan = find(~isnan(spdSlice));
spdSlice(indnan) = interp1(glon(latInd,indnotnan),spdSlice(indnotnan), ...
    glon(latInd,indnan));
% Running smooth vel 
smthSpdSlice = movmean(spdSlice,5,'omitmissing');
% smthSpdSlice = spdSlice;
% Find all the vel > 0.7 in offshore region 
smthSpdSlice(ddist_from_coast(iFirst,:) < ddist_from_coastVal) = NaN;
indmaxVel = find( smthSpdSlice==max(smthSpdSlice) & VSlice < 0 );

if isempty(indmaxVel)
    disp('Current goes upward!!');
    return
end

% Core location
lonCoreVel(latInd) = glon(latInd,indmaxVel);
latCoreVel(latInd) = glat(latInd,indmaxVel);

% Put into jet coordinate
d = nan(1,length(jetxAxInd));
for ii=1:length(glon(latInd,:))
    d(ii) = m_lldist([glon(latInd,ii),lonCoreVel(latInd)], ...
        [glat(latInd,ii),latCoreVel(latInd)]);
end
dInd = zeros(1,length(d));
% Left hand side
dInd(1:indmaxVel-1) = round(d(1:indmaxVel-1)/dRx)*-1;
% Right hand side
dInd(indmaxVel:end) = round(d(indmaxVel:end)/dRx);
% Then transform geo coord to jet coord
jetAxVel(latInd,ismember(jetxAxInd,dInd)) = spdSlice(ismember(dInd,jetxAxInd));
% Core vel 
coreVel(latInd) = spdSlice(indmaxVel);
% Jet direction
[spd,ang] = uv2spdir(USlice(indmaxVel-1:indmaxVel+1),VSlice(indmaxVel-1:indmaxVel+1));
coreAngle(latInd) = nanmean(ang);
% Jet width: 
% Search for the points that vel reduces by 50% of the core vel 
% Offshore side: 
minJetEdgeBound = coreVel(latInd)*0.5;
minJetEdgeBound(minJetEdgeBound < 0.5) = 0.5;

indSearch = indmaxVel:length(smthSpdSlice);
indW2 = find( smthSpdSlice(indSearch) < minJetEdgeBound );
if ~isempty(indW2)
    indW2 = indW2(1);
    % Position 
    lonCoreWidth(latInd,2) = glon(latInd,indSearch(indW2));
    latCoreWidth(latInd,2) = glat(latInd,indSearch(indW2));
else
    indW2 = [];
    % Position 
    lonCoreWidth(latInd,2) = NaN;
    latCoreWidth(latInd,2) = NaN;
end
% Onshore side: 
indSearch = 1:indmaxVel;
indW1 = find( smthSpdSlice(indSearch) < minJetEdgeBound );
if ~isempty(indW1)
    indW1 = indW1(end);
    % Position 
    lonCoreWidth(latInd,1) = glon(latInd,indSearch(indW1));
    latCoreWidth(latInd,1) = glat(latInd,indSearch(indW1));
else
    indW1 = [];
    % Position 
    lonCoreWidth(latInd,1) = NaN;
    latCoreWidth(latInd,1) = NaN;
end

%% Seach downstream
lonCorePrev = lonCoreVel(iFirst);
latCorePrev = latCoreVel(iFirst);

for latInd = iFirst+1:size(glat,1)
    
    USlice = squeeze(U(latInd,:));
    VSlice = squeeze(V(latInd,:));
    
    spdSlice = sqrt( USlice.^2 + VSlice.^2 );
    % Check nb of data
    if sum(isfinite(spdSlice)) < 5
        continue
    end
    % Interpolate the gaps 
    indnan = find(isnan(spdSlice));
    indnotnan = find(~isnan(spdSlice));
    spdSlice(indnan) = interp1(glon(latInd,indnotnan),spdSlice(indnotnan), ...
        glon(latInd,indnan));
    % Running smooth vel 
    smthSpdSlice = movmean(spdSlice,5,'omitmissing');
    % smthSpdSlice = spdSlice;

    % Identify maximum filtered vel in the eac 
    % Check the core location with the core at the previous location
    d = deg2km(sqrt( (glon(latInd,:)-lonCorePrev).^2 ...
        + (glat(latInd,:)-latCorePrev).^2 ));

    indmaxVel1 = find( smthSpdSlice > 0.6 & d < searchRadius*dRx  & VSlice < 0 );
    if isempty(indmaxVel1)
        continue
    end
    indmaxVel2 = find( smthSpdSlice == max(smthSpdSlice(indmaxVel1)) );
    indmaxVel  = indmaxVel2;
    
    if isempty(indmaxVel) || ddist_from_coast(latInd, indmaxVel) < ddist_from_coastVal
        continue
    end

    % 
    % Core location
    lonCoreVel(latInd) = glon(latInd,indmaxVel);
    latCoreVel(latInd) = glat(latInd,indmaxVel);
    % Put into jet coordinate
    d = nan(1,length(jetxAxInd));
    for ii=1:length(glon(latInd,:))
        d(ii) = m_lldist([glon(latInd,ii),lonCoreVel(latInd)], ...
            [glat(latInd,ii),latCoreVel(latInd)]);
    end
    dInd = zeros(1,length(d));
    % Left hand side
    dInd(1:indmaxVel-1) = round(d(1:indmaxVel-1)/dRx)*-1;
    % Right hand side
    dInd(indmaxVel:end) = round(d(indmaxVel:end)/dRx);
    % Then transform geo coord to jet coord
    jetAxVel(latInd,ismember(jetxAxInd,dInd)) = spdSlice(ismember(dInd,jetxAxInd));
    % Core vel 
    coreVel(latInd) = spdSlice(indmaxVel);
    % Jet direction
    [spd,ang] = uv2spdir(USlice(indmaxVel-1:indmaxVel+1),VSlice(indmaxVel-1:indmaxVel+1));
    coreAngle(latInd) = nanmean(ang);
    % Jet width: 
    minJetEdgeBound = coreVel(latInd)*0.5;
    minJetEdgeBound(minJetEdgeBound < 0.5) = 0.5;
    % Search for the points that vel reduces by 50% of the core vel 
    % Offshore side: 
    indSearch = indmaxVel:length(smthSpdSlice);
    indW2 = find( smthSpdSlice(indSearch) < minJetEdgeBound );
    if ~isempty(indW2)
        indW2 = indW2(1);
        % Position 
        lonCoreWidth(latInd,2) = glon(latInd,indSearch(indW2));
        latCoreWidth(latInd,2) = glat(latInd,indSearch(indW2));
    else
        indW2 = [];
        % Position 
        lonCoreWidth(latInd,2) = NaN;
        latCoreWidth(latInd,2) = NaN;
    end
    % Onshore side: 
    indSearch = 1:indmaxVel;
    indW1 = find( smthSpdSlice(indSearch) < minJetEdgeBound );
    if ~isempty(indW1)
        indW1 = indW1(end);
        % Position 
        lonCoreWidth(latInd,1) = glon(latInd,indSearch(indW1));
        latCoreWidth(latInd,1) = glat(latInd,indSearch(indW1));
    else
        indW1 = [];
        % Position 
        lonCoreWidth(latInd,1) = NaN;
        latCoreWidth(latInd,1) = NaN;
    end

    lonCorePrev = lonCoreVel(latInd);
    latCorePrev = latCoreVel(latInd);

end

%% Search upstream
lonCorePrev = lonCoreVel(iFirst);
latCorePrev = latCoreVel(iFirst);

for latInd = iFirst-1:-1:1
    
    USlice = squeeze(U(latInd,:));
    VSlice = squeeze(V(latInd,:));
    
    spdSlice = sqrt( USlice.^2 + VSlice.^2 );
    % Check nb of data
    if sum(isfinite(spdSlice)) < 5
        continue
    end
    % Interpolate the gaps 
    indnan = find(isnan(spdSlice));
    indnotnan = find(~isnan(spdSlice));
    spdSlice(indnan) = interp1(glon(latInd,indnotnan),spdSlice(indnotnan), ...
        glon(latInd,indnan));
    % Running smooth vel 
    smthSpdSlice = movmean(spdSlice,5,'omitmissing');
    % smthSpdSlice = spdSlice;

    % Identify maximum filtered vel in the eac 
    % Check the core location with the core at the previous location
    d = deg2km(sqrt( (glon(latInd,:)-lonCorePrev).^2 ...
        + (glat(latInd,:)-latCorePrev).^2 ));

    indmaxVel1 = find( smthSpdSlice > 0.6 & d < searchRadius*dRx  & VSlice < 0 );
    if isempty(indmaxVel1)
        continue
    end
    indmaxVel2 = find( smthSpdSlice == max(smthSpdSlice(indmaxVel1)) );
    indmaxVel  = indmaxVel2;
    
    if isempty(indmaxVel) || ddist_from_coast(latInd, indmaxVel) < ddist_from_coastVal
        continue
    end

    % 
    % Core location
    lonCoreVel(latInd) = glon(latInd,indmaxVel);
    latCoreVel(latInd) = glat(latInd,indmaxVel);
    % Put into jet coordinate
    d = nan(1,length(jetxAxInd));
    for ii=1:length(glon(latInd,:))
        d(ii) = m_lldist([glon(latInd,ii),lonCoreVel(latInd)], ...
            [glat(latInd,ii),latCoreVel(latInd)]);
    end
    dInd = zeros(1,length(d));
    % Left hand side
    dInd(1:indmaxVel-1) = round(d(1:indmaxVel-1)/dRx)*-1;
    % Right hand side
    dInd(indmaxVel:end) = round(d(indmaxVel:end)/dRx);
    % Then transform geo coord to jet coord
    jetAxVel(latInd,ismember(jetxAxInd,dInd)) = spdSlice(ismember(dInd,jetxAxInd));
    % Core vel 
    coreVel(latInd) = spdSlice(indmaxVel);
    % Jet direction
    [spd,ang] = uv2spdir(USlice(indmaxVel-1:indmaxVel+1),VSlice(indmaxVel-1:indmaxVel+1));
    coreAngle(latInd) = nanmean(ang);
    % Jet width: 
    minJetEdgeBound = coreVel(latInd)*0.5;
    minJetEdgeBound(minJetEdgeBound < 0.5) = 0.5;
    % Search for the points that vel reduces by 50% of the core vel 
    % Offshore side: 
    indSearch = indmaxVel:length(smthSpdSlice);
    indW2 = find( smthSpdSlice(indSearch) < minJetEdgeBound );
    if ~isempty(indW2)
        indW2 = indW2(1);
        % Position 
        lonCoreWidth(latInd,2) = glon(latInd,indSearch(indW2));
        latCoreWidth(latInd,2) = glat(latInd,indSearch(indW2));
    else
        indW2 = [];
        % Position 
        lonCoreWidth(latInd,2) = NaN;
        latCoreWidth(latInd,2) = NaN;
    end
    % Onshore side: 
    indSearch = 1:indmaxVel;
    indW1 = find( smthSpdSlice(indSearch) < minJetEdgeBound );
    if ~isempty(indW1)
        indW1 = indW1(end);
        % Position 
        lonCoreWidth(latInd,1) = glon(latInd,indSearch(indW1));
        latCoreWidth(latInd,1) = glat(latInd,indSearch(indW1));
    else
        indW1 = [];
        % Position 
        lonCoreWidth(latInd,1) = NaN;
        latCoreWidth(latInd,1) = NaN;
    end

    lonCorePrev = lonCoreVel(latInd);
    latCorePrev = latCoreVel(latInd);

end

return
end