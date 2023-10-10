% Determing a threshold for size of dunes to search for based on amplitude
% by fitting a gamma distribution to the depth data and choosing a
% threshold one standard deviation less than the mean or peak.

%% If export == 1 threshold data will be saved as ascii grid files for ArcGIS
% export = 0;

%% Determing a threshold for size of dunes to search for based on amplitude
% by fitting a gamma distribution to the depth data and choosing a
% threshold one standard deviation less than the mean or peak or the mean
% as the threshold

ampvector = amplitude(~isnan(amplitude)); %only use the amplitudes from the lee sides

meanamp = mean(ampvector);
figure(10)
histfit(ampvector,100,'gamma')
pd = fitdist(ampvector,'gamma');
s = std(pd);
m = mean(pd);

threshold = 0.1; %m + s %threshold that defines small and large dunes

ampthresh = zeros(size(amplitude)); %threshold is defined by dune height

for j = 1:size(amplitude,2);
    for i = 1:size(amplitude,1);
        if amplitude(i,j) > threshold
            ampthresh(i,j) = amplitude(i,j);
            meanLeethresh(i,j) = meanleeangle(i,j);
            maxLeethresh(i,j) = maxlee(i,j);
            depththresh(i,j) = depth(i,j);
            Maxpercentagelee(i,j) = Maxpercentagelee(i,j);
            slopedirection(i,j) = slopedirection(i,j);

        else
           ampthresh(i,j) = NaN;
           meanLeethresh(i,j) = NaN;
           maxLeethresh(i,j) = NaN;
           depththresh(i,j) = NaN;
           maxpercentagelee(i,j) = NaN;
           slopedirection(i,j) = NaN;

        end
    end
end


for j = 1:size(amplitude,2);
    for i = 1:size(amplitude,1);
        if amplitude(i,j) < threshold
            ampSP(i,j) = amplitude(i,j);
            meanLeeSP(i,j) = meanleeangle(i,j);
            maxLeeSP(i,j) = maxlee(i,j);
            depthSP(i,j) = depth(i,j);
            slopedirectionSP(i,j) = slopedirection(i,j);

        else
           ampSP(i,j) = NaN;
           meanLeeSP(i,j) = NaN;
           maxLeeSP(i,j) = NaN;
           depthSP(i,j) = NaN;
           slopedirectionSP(i,j) = NaN;

        end
    end
end


keep ampthresh meanLeethresh maxLeethresh depththresh Maxpercentagelee slopedirection...
    ampSP meanLeeSP maxLeeSP depthSP slopedirectionSP flow flowvar export res DEMfile...
    X2 Y2
