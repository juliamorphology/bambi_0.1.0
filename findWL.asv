% okay so here we are gonna use the location points of superimposed and
% formative dunes from the threshold. we will index them something
% different like 1 and 2. then we will do a difference. then we will find
% the distance between -1 a(0-1) for superimposed and 1 (2-1) for
% formative.


LocPt = zeros(size(depththresh));    %empty matrix the size fo the depth matrix

indF = ~isnan(depththresh);         %Location indices of formative dunes
indS = ~isnan(depthSP);             %Location infices of superimposed dunes

LocPt(indF) = 2;                    %Formative dunes are designate by a 2
LocPt(indS) = 1;                    %Superimposed dunes are designate by a 1

LocPt = imrotate(LocPt, flow);      %Rota te the matrix of location so we
                                    %can run down columns
Fwavelength = zeros(size(LocPt));
Swavelength = zeros(size(LocPt));

for j = 1:size(LocPt,2)
    for i = 1:size(LocPt,1)
        Fpt = find(LocPt(:,j) == 2);        %Location indices of formative dunes
        Spt = find(LocPt(:,j) == 1);        %Location indices of superimposed dunes

        Fd = diff(Fpt);    % find space between each Formative trough pt
        Sd = diff(Spt);    % find space between each Superimposed trough pt

        Fwavelength(Fpt(1:size(Fd)),j) = Fd * res;
        Swavelength(Spt(1:size(Sd)),j) = Sd * res;


    end
end

%This makes sure all the values in the wavelenth calculation are either a
%value, or NaN
Fwavelength(find(Fwavelength == 0)) = NaN;
Fwavelength(find(Fwavelength == 1)) = NaN;
Swavelength(find(Swavelength == 0)) = NaN;

%Here we rotate the calculated values back to the original orientation
FormativeWL = imrotate(Fwavelength,-flow);
SuperimposedWL = imrotate(Swavelength,-flow);

%Now because the X and Y value were trimmed earlier in the code, we know we
%trimmed them by even amounts from the end margins of the matrix. So we
%find the difference in size between the calculated matrix and the X and Y
%matrix that corresponds and we divide by two so we can subtract evenly
%from each margin
row = abs(size(FormativeWL,1)-size(X2,1));
row = row/2;
col = abs(size(FormativeWL,2)-size(X2,2));
col = col/2;
row = abs(size(SuperimposedWL,1)-size(X2,1));
row = row/2;
col = abs(size(SuperimposedWL,2)-size(X2,2));
col = col/2;
%Now we fix the indices of which to pull X and Y from by trimming the
%margins by half the difference in size i.e., we trim top, bottom, left,
%and right sides of the matrix.
FormativeWL = FormativeWL(row+1:end-row,col+1:end-col); %add 1 to the beginning because i starts at 1
SuperimposedWL = SuperimposedWL(row+1:end-row,col+1:end-col);



FormativeWL(find(FormativeWL == 0)) = NaN;
SuperimposedWL(find(SuperimposedWL == 0)) = NaN;
