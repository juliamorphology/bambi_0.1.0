%% Bedform analysis of a DEM file. Outputs amplitude, leeside (mean/max) values, wavelength and depth of large and superimposed dunes)
%DEMfile = file location of the .txt file (string)
%flow = flow direction within the channel (double)
%export = the choice to export .asc files (0=no, 1=yes)
%Smalldunes/Largedunes are output tables with X, Y, Amplitude, Leeside
%angle (mean/max),sloping direction, wavelength, and depth as columns for
%each dune.



function [Smalldunes,Largedunes] = BAMBI(DEMfile,flow,flowvar,export)
%% add the matlab path to use keep, arcgiswrite, slopeCalc, Bedform threshold, and findwavelength functions
% change this to wherever you store your most used functions
addpath(genpath('/Users/cisneros/Documents/MATLAB'));
addpath(genpath('/Users/cisneros/Documents/GitHub/bambi_0.1.0'));
%% Upload the spatial matrix of slope and slope direction for pc
% code needs txt files of slope, aspect, and depth in the same folder
% txt files should be created using ArcMap tool "Raster to ASCII"
% txt files should all be for the same area and same resolution

%folder = 'J:\BedformAnalysis\Dune_data\Parana\'; %folder path where txt files are stored


depth = dlmread([DEMfile],' ',6,1); %i had to change this ebcause it looks like its reading the 1st column as "0" change 1 to 0 if it doesnt work
% Change -9999 NoData value to NaN
depth(depth==0) = NaN;
depth(depth==-9999) = NaN;
depth(depth==-999.99) = NaN;

% read in the first 6 lines of the slope txt file to get the lower left
% Easting and Northing coordiantes, and the resolution
[info.category info.value] = textread([DEMfile],'%s %f',6);
res = info.value(5,1); %raster resolution

%% Rotate depth to work along profile lines
depth_rotated = imrotate(depth,flow); %rotated matrix of depth grid
depth_rotated(depth_rotated == 0)= NaN;
depth_rotated(depth_rotated == -9999)= NaN;

depth1 = depth_rotated; %if not smoothing depth1 = depth rotated

clear depth_rotated
clear depth

%% Use the slopeCALC function that I wrote. Makes slope and direction
% raster calculations using same methods as ArcMap
[slope_rotated, direction_rotated] = slopeCALC(depth1,res);

     %% Create Easting and Northing Vectors and location grids
Easting = zeros(1,info.value(1,1)); %create zeros vector of the same size as the data
Northing = zeros(info.value(2,1),1); %create zeros vector of the same size as the data
Easting(1,1) = info.value(3,1); % populate x-left value
Northing(1,1) = info.value(4,1) + info.value(5,1)*(info.value(2,1)-1); %populate y-top value

for j = 2:info.value(1,1);
    Easting(1,j) = Easting(1,j-1)+info.value(5,1); %fill in the rest of the easting matrix
end
for i = 2:info.value(2,1);
    Northing(i,1) = Northing(i-1,1)-info.value(5,1); %fill in the rest of the northing matrix
end
[X,Y]=meshgrid(Easting,Northing); %creates grids of Easting and Northing coordinates


clear Easting Northing i j


 %% rotate the depth, slope, aspect, and location grids.
% % The flow direction will now be pointing upwards towards 0 degrees

X_rotated = imrotate(X,flow); %rotated matrix of x grid points
Y_rotated = imrotate(Y,flow); %rotated matrix of y grid points



%% Create bounds for designating lee side and stoss side
% Define some degree range in the flow direction in which values that fall between that
% range are lee sides and values outside are stoss sides.

a = 360 - flowvar;  %The lower bound for directions of flow
b = 0 + flowvar; %The upper bound for directions of flow

% make sure a and b are within 0 to 360 degrees
if a < 0;
    a = 360 + a;
elseif b > 360;
    b = b - 360;
end

%% Split Slope, Direction, and Depth grids by lee and stoss sides

% create zeros matrix for splitting slope, aspect (dir), and depth into lee
% and stoss sides
slope_lee = zeros(size(slope_rotated));

for i = 1:size(slope_rotated,1);
    for j = 1:size(slope_rotated,2);

        % if the aspect of this grid cell is within the range of flow directions
        % populate that grid cell in the lee side matrices with the slope, depth,
        % and aspect (direction) values. if not, fill that cell with NaN
        if direction_rotated(i,j) >= a || direction_rotated(i,j) <= b
            slope_lee(i,j) = slope_rotated(i,j);
        else
            slope_lee(i,j) = NaN;
        end
    end

end
clear i j

%% Identify transitions from lee to stoss sides
nanID = isnan(slope_lee); %find where the stoss sides are located (i.e. where the lee side matrix is NaN)

for k = 1:size(nanID,1)-1;
    for j = 1:size(nanID,2);
        slopebreaks(k,j) = nanID(k,j)-nanID(k+1,j);  %go down each column of the nanID logical matrix and subtract each cell from the one above it
        % could probably get rid of one of the for loops and replace with
        % the diff function
    end
end
clear k j




%% mean/max lee angles, mean/max stoss angles, wavelengths, amplitudes along flow paths

% since the matrices were rotated, the cell size changed...need
% this to get wavelengths
if flow <=90;
    theta = flow;
elseif flow>90 && flow<=180;
    theta = 180-flow;
elseif flow>180 && flow<=270;
    theta = flow-180;
elseif flow>270 && flow<=360;
    theta =360-flow;
end

if flow == 270 || flow == 0 || flow == 90 || flow == 180;
    cell = res; %vertical size of each cell after rotation;
else
    cell = res/cosd(theta); %vertical size of each cell after rotation;
end


%% calculate mean (uses trimmean and excludes top and bottom 5%) & max lee side angles, widths and amplitudes of each

[sizeX,sizeY] = size(slope_rotated);
cWL = NaN(sizeX,sizeY,'like',slope_rotated);
tWL = NaN(sizeX,sizeY,'like',slope_rotated);
amp = NaN(sizeX,sizeY,'like',slope_rotated);
leemeanval = NaN(sizeX,sizeY,'like',slope_rotated);
leemaxval = NaN(sizeX,sizeY,'like',slope_rotated);
stossmeanval = NaN(sizeX,sizeY,'like',slope_rotated);
stossmaxval = NaN(sizeX,sizeY,'like',slope_rotated);
maxamp = NaN(sizeX,sizeY,'like',slope_rotated);
slopedirection = NaN(sizeX,sizeY,'like',slope_rotated);
%Now all measured values of lee side will be assigned to the trough point INFRONT of
%the dune that the measured values correspond to. This is because we
%rotated the field so flow points upwards and we work down a row. all
%values of stoss side will be at the crest point of the dune we measure.
%Amplitude and wavelength will be at both points and should loosely
%correspond

for j = 1:size(slope_rotated,2)
    tloc = find(slopebreaks(:,j) == 1); %locations of troughs (lee to stoss)
    cloc = find(slopebreaks(:,j) == -1); %locations of crests (stoss to lee)


      %for lee side measurements
    for i = 1:size(tloc)
        if tloc(i) < cloc(i)    %the first location is a trough location
            amp(tloc(i),j) = abs(max(depth1((tloc(i):cloc(i)),j))-min(depth1((tloc(i):cloc(i)),j)));
            topfive = ceil((.05*numel(tloc(i):cloc(i))));
            leemeanval(tloc(i),j) = mean(slope_rotated((tloc(i)+topfive:cloc(i)-topfive),j));
            leemaxval(tloc(i),j) = nanmax(slope_rotated((tloc(i):cloc(i)),j));
            maxind = tloc(i) + find(slope_rotated((tloc(i):cloc(i)),j) == max(slope_rotated((tloc(i):cloc(i)),j))) - 1; %finding
            %the exact location of the max ind within the column. we need to find where the max is within each section,
            %then add the indice to the trough location. but because counting starts at 1 in matlab we have to subtract 1.
            maxamp(tloc(i),j) = depth1(maxind(1,1),j)-depth1(tloc(i),j);
            slopedirection(tloc(i),j) = mean(direction_rotated((tloc(i)+topfive:cloc(i)-topfive),j));
        elseif cloc(i) < tloc(i) && i < size(cloc)
            amp(tloc(i),j) = abs(max(depth1((tloc(i):cloc(i+1)),j))-min(depth1((tloc(i):cloc(i+1)),j)));
            topfive = ceil((.05*numel(tloc(i):cloc(i+1))));
            leemeanval(tloc(i),j) = mean(slope_rotated((tloc(i)+topfive:cloc(i+1)-topfive),j));
            leemaxval(tloc(i),j) = nanmax(slope_rotated((tloc(i):cloc(i+1)),j));
            maxind = tloc(i) + find(slope_rotated((tloc(i):cloc(i+1)),j) == max(slope_rotated((tloc(i):cloc(i+1)),j))) - 1; %finding
            %the exact location of the max ind within the column. we need to find where the max is within each section,
            %then add the indice to the trough location. but because counting starts at 1 in matlab we have to subtract 1.
            maxamp(tloc(i),j) = depth1(maxind(1,1),j)-depth1(tloc(i),j);
            slopedirection(tloc(i),j) = mean(direction_rotated((tloc(i)+topfive:cloc(i+1)-topfive),j));
        else
            amp(tloc(i),j) = NaN;
            leemeanval(tloc(i),j) = NaN;
            leemaxval(tloc(i),j) = NaN;
            maxamp(tloc(i),j) = NaN;
            slopedirection(tloc(i),j) = NaN;
        end

    end


end


maxpercentlee = (maxamp ./ amp)*100;
maxpercentlee(maxpercentlee == 0) = NaN;
maxpercentlee(maxpercentlee == 100) = NaN;
Maxpercentagelee = abs(maxpercentlee);
  %% Now rotate all our data matrices back to their true orientation

X2 = imrotate(X_rotated,-flow);
Y2 = imrotate(Y_rotated,-flow);
meanlee2 = imrotate(leemeanval,-flow);
meanstoss2 = imrotate(stossmeanval,-flow);
amplitude = imrotate(amp,-flow);
maxlee = imrotate(leemaxval,-flow);
maxstoss = imrotate(stossmaxval,-flow);
depth2 = imrotate(depth1,-flow);
Maxpercentagelee = imrotate(Maxpercentagelee,-flow);
slopedirection = imrotate(slopedirection, -flow);


%% change zeros and -9999 to NaN in data matrices
meanlee2(meanlee2==0) = NaN;
meanstoss2(meanstoss2==0) = NaN;
amplitude(amplitude == 0) = NaN;
amplitude(amplitude > 900) = NaN;
maxlee(maxlee == 0) = NaN;
maxstoss(maxstoss == 0) = NaN;
depth2(depth2 == 0) = NaN;
Maxpercentagelee(Maxpercentagelee == 0) = NaN;
slopedirection(slopedirection == 0) = NaN;
Y2(Y2 == 0) = NaN;
X2(X2 == 0) = NaN;


%rename to make things nice
meanleeangle=meanlee2;
meanstossangle=meanstoss2;
depth = depth2;

% only save amplitude for lee sides (1 amp per dune)
leeind = find(meanleeangle);

keep meanleeangle meanstossangle amplitude maxlee maxstoss depth...
    Maxpercentagelee slopedirection leeind flow DEMfile flowvar res export...
    X2 Y2

%% run threshold, plotting, and wavelength scripts & compile/save into two tables.
run('BedformThreshold.m')
run('findWL.m')

%create tables to save small superimposed dunes and large dunes
smallind = find(~isnan(meanLeeSP));
Smalldunes = [X2(smallind),Y2(smallind),ampSP(smallind),meanLeeSP(smallind),maxLeeSP(smallind),slopedirectionSP(smallind),SuperimposedWL(smallind),depthSP(smallind)];

largeind = find(~isnan(meanLeethresh));
Largedunes = [X2(largeind),Y2(largeind),ampthresh(largeind),meanLeethresh(largeind),maxLeethresh(largeind),slopedirection(largeind),FormativeWL(largeind),depththresh(largeind),Maxpercentagelee(largeind)];

%clean the data for large and small dunes
aspect_small = Smalldunes(:,3)./Smalldunes(:,7);
aspect_large = Largedunes(:,3)./Largedunes(:,7);

small_ind = find(~isnan(aspect_small));
large_ind = find(~isnan(aspect_large));

Smalldunes = Smalldunes(small_ind,:);
Largedunes = Largedunes(large_ind,:);


 if export == 1 %to save data files as table in the same location as the DEM
    %move to the folder the DEM is from
    [pathstr] = fileparts(DEMfile);
    cd(pathstr)

    %Make new table in the main folder
    mkdir('DataTables')
    cd ([pwd '/DataTables/'])
    %save tables as txt files
    dlmwrite('Smalldunes.txt',Smalldunes,'precision',10)
    dlmwrite('Largedunes.txt',Largedunes,'precision',10)
 else

    'No export'

end


load gong
%sound(y,Fs)

end
