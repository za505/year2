%dataPlots.m 
%Author: Zarina Akbary 
%Date: 30 March 2021 
%Purpose: to gather all the data in one place and plot them to make sound conclusions.

clear, close all
%%%%%%%%%%%%%

%% User Input: 647
%first, let's start with the easier one, 647
maindir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/data_analysis'];
filename={['03172021_Exp1'] ['03252021_Exp2'] ['03292021_Exp2']};
channel=['_647'];

%we'd never included index as a variable in the previous scripts in our
%workflow, so we have to hard code it. for the Mg2+ gradient experiments,
%it is this ['last frame without dye' 'last frame of first perfusion' 'last
%frame of second perfusion' 'last frame of third perfusion' 'last frame of
%fourth perfusion']

%%%%%%%%%%%%%%%%
for f=1:length(filename)
    

    %first, let's open the mat file and get the ratios
    file=char(filename(f));
    cd(maindir)
    load([file '_BTfluoAVG' channel '.mat'], 'ratio');
    mgConc=nan(width(ratio), 1);
    experiment=transpose(repelem(filename(f), width(ratio)));
    dye=transpose(repelem("647", width(ratio)));
    
    load([file '_BTfluoAVG' channel '.mat'], 'mgRange');
    
    if f==1
        index=[18 30 42 54 61];
    elseif f==2
        index=[19 30 40 54 66];
    elseif f==3
        index=[20 29 42 53 67];
    end
    
     for i=1:length(mgRange)
        mgConc(index(i)+1:index(i+1), 1)=mgRange(i);
     end
     
     dataRatio=[array2table(experiment), array2table(dye), array2table(mgConc), array2table(transpose(ratio))];
     writetable(dataRatio,[file channel '_dataRatio.csv']);
     
     %now, the tricky part--concatenating the length values
     load([file '_BTfluoAVG' channel '.mat'], 'basename');
     lcellComp=[];
     
    for b=1:length(basename)
        
        base=char(basename(b))
        
        temp=struct2cell(load([base '_BTphase.mat'], 'lcell'));
        temp=temp{1,1};
        lcellComp=[lcellComp; temp];
        
    end
    
    dataLength=[array2table(experiment), array2table(dye), array2table(mgConc), array2table(transpose(lcellComp))];
    writetable(dataLength,[file channel '_dataLength.csv']);
     
end

  