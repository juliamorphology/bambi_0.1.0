%% gonna try to calculate slope from a random raster using a 3x3 moving window
% testing against ArcMap slope and aspect tool on juaturana dataset: min
% error -1e-4 degree, max error 3.22 degree, mean error 3.2e-6
function [slope,direction] = slopeCALC(DEM,res)

[sizeX,sizeY] = size(DEM);

%if we work within a 3x3 window [a,b,c;d,e,f;g,h,i] where e is the cell we
%populate with a slope value

a = DEM((1:sizeX-2),(1:sizeY-2));
b = DEM((1:sizeX-2),(2:sizeY-1));
c = DEM((1:sizeX-2),(3:sizeY));
d = DEM((2:sizeX-1),(1:sizeY-2));
f = DEM((2:sizeX-1),(3:sizeY));
g = DEM((3:sizeX),(1:sizeY-2));
h = DEM((3:sizeX),(2:sizeY-1));
i = DEM((3:sizeX),(3:sizeY));

%calculate the rise over run through dz/dx and dz/dy
dzdx = NaN(sizeX,sizeY,'like',DEM);
dzdy = NaN(sizeX,sizeY,'like',DEM);
dzdx(2:sizeX-1,2:sizeY-1) = ((c+(2*f)+i)-(a+(2*d)+g))./(8*res);
dzdy(2:sizeX-1,2:sizeY-1) = ((g+(2*h)+i)-(a+(2*b)+c))./(8*res);

% riserun = zeros(sizeX,sizeY);
%fill only e values
riserun = sqrt((dzdx.^2)+(dzdy.^2));

%convert radians to degrees for slope and find direction

slope = atan(riserun)* (180/pi);
direction = atan2(dzdy,dzdx)* (180/pi);

%convert the +/- 0-180 values to azimuth 0-360
ind1 = find(direction(:) >= 90 & direction(:) < 180);
ind2 = find(direction(:) < 90 & direction(:) >= -180);
direction(ind1) = direction(ind1) - 90;
direction(ind2) = direction(ind2) + 270;
end
