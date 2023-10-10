% to run the bedform analysis from an excel parameter file
% addpath('C:\Users\cisneros\Documents\MATLAB\MATLAB')
parameterfile = '/Users/cisneros/Library/CloudStorage/OneDrive-VirginiaTech/LakeAustinProject/Julia/Split_ascii_for_bambi/data0623/parameterfile.xlsx';
[num,txt,raw] = xlsread(parameterfile);

DEMfile = raw(2:end,1);
flow = raw(2:end,2);
flowvar = raw(2:end,3);
export = raw(2:end,4);

 for i = 1:numel(flow)
       print = strcat('analyzing',' ',num2str(i));
       disp(print)
[Smalldunes,Largedunes] = BAMBI(DEMfile{i},flow{i},flowvar{i},export{i}); 

if export{i} == 1
    print = strcat('saving',' ',num2str(i));
    disp(print)
    %probably need to change somewhere here to save with the real name but
    %in the folder name. that way all data is saved in folders but the
    %large duens files are unique so I can concatenate them later..
    foldername = [DEMfile{i}(1:end-4)];
    if exist(foldername) ~= 7
        mkdir(foldername);
        cd(foldername)
        flowstr = num2str(flow{i});
        flowvarstr = num2str(flowvar{i});
        smallfilename =  ['Smalldunes',flowstr,'_',flowvarstr,'.txt'];
        dlmwrite(smallfilename,Smalldunes,'precision',10,'delimiter',' ');

        largefilename =  ['Largedunes',flowstr,'_',flowvarstr,'.txt'];
        dlmwrite(largefilename,Largedunes,'precision',10,'delimiter',' ');
    else
        cd(foldername)
        flowstr = num2str(flow{i});
        flowvarstr = num2str(flowvar{i});
        smallfilename =  ['Smalldunes',flowstr,'_',flowvarstr,'.txt'];
        dlmwrite(smallfilename,Smalldunes,'precision',10,'delimiter',' ');

        largefilename =  ['Largedunes',flowstr,'_',flowvarstr,'.txt'];
        dlmwrite(largefilename,Largedunes,'precision',10,'delimiter',' ');
    end
end

 keep DEMfile flow flowvar export i

disp(i)


cd ../../
end
