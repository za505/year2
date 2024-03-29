%dataPlots.m 
%Author: Zarina Akbary 
%Date: 30 March 2021 
%Purpose: to gather all the data in one place and plot them to make sound conclusions.

clear, close all
%%%%%%%%%%%%%

%% User Input: 647
maindir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/data_analysis/FSS'];
filename={['03172021_Exp1'] ['03252021_Exp1'] ['03262021_Exp1'] ['03292021_Exp1']};
channel=['_FSS'];
recrunch=1;
%we'd never included index as a variable in the previous scripts in our
%workflow, so we have to hard code it. for the Mg2+ gradient experiments,
%it is this ['last frame without dye' 'last frame of first perfusion' 'last
%frame of second perfusion' 'last frame of third perfusion' 'last frame of
%fourth perfusion']

%%%%%%%%%%%%%%%%
    
%     %first, 03172021_Exp1
%     file=char(filename(4));
%     cd(maindir)
%     dataLength=readtable([file channel '_dataLength.csv']);
%     lengths=transpose(table2array(dataLength(:,4:end)));
%     load([file '_BTfluoAVG' channel '.mat'], 'time');
%     
%     %avgLength=mean(lengths, 1, 'omitnan');
%     %stdLength=std(lengths, 0, 1, 'omitnan');
%     
%   %let's plot the average length vs Mg2+ 
%    figure 
% clf
% hold on
% for i=1:height(lengths)  
%     plot(time(1:end),lengths(i,1:end)) 
% end
%     xlabel('Time (s)')
%     ylabel('Length (\mum)')
%      title('Cell Length vs. Time')
%     ylim([0 Inf])
%     xline(60, '--k', '*PBS + 5% detergent') %frame 1-18
% xline(110, '--k', '*PBS + FSS') %frame 11-21
% xline(220, '--k', '*PBS + FSS + 6.66 mM Mg^{2+}') %frame 22-33
% xline(340, '--k', '*PBS + FSS + 12.33 mM Mg^{2+}') %frame 34-44 
% xline(450, '--k', '*PBS + FSS + 20 mM Mg^{2+}') %frame 45-59
%     fig2pretty 
%     saveas(gcf, [file '_ltraces' channel '.png'])
    
    
if recrunch==0
    
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
        index=[23 33 45 58 67];
    elseif f==3
        index=[17 28 40 52 64];
    elseif f==4
        index=[10 21 33 44 59];
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

cd(maindir)
 dataRatio=readtable('mgRatio_647.csv');
 dataLength=readtable('mgLength_647.csv');
 mgRange=[0, 6.66, 12, 12.33, 15, 20];
 avgRatio=[];
 avgLength=[];
 
 ratios=table2array(dataRatio(:,4:end));
 lengths=table2array(dataLength(:,4:end));
 
 for i=1:length(mgRange)
    avgRatio(i, :)=mean(ratios(dataRatio.mgConc==mgRange(i),:),1, 'omitnan');
    avgLength(i, :)=mean(lengths(dataLength.mgConc==mgRange(i),:),1, 'omitnan');
 end
 
 avgMgRatio=mean(avgRatio, 2, 'omitnan');
 stdMgRatio=std(avgRatio,0,2, 'omitnan');
 
 avgMgLength=mean(avgLength, 2, 'omitnan');
 stdMgLength=std(avgLength,0,2, 'omitnan');

else
    cd(maindir)
    load(['dataPlots' channel '.mat'])
end

 %let's plot the average ratio vs Mg2+ 
    figure, hold on
    title('Average Intensity/Background vs Mg^{2+} Concentration')
    errorbar(mgRange([1:4,7:8]),avgMgRatio([1:4,7:8]),stdMgRatio([1:4,7:8]), 'both', 'o')
    xlabel('Mg^{2+} concentration (mM)')
    ylabel('Intensity/Background')
    yline(1, '--k')
    xlim([-2 22])
    ylim([0 Inf])
    %xticks(mgRange)
    fig2pretty 
    %saveas(gcf, ['avgMgRatio' channel '.fig'])
% %     
%     %let's plot the average length vs Mg2+ 
%     figure, hold on
%     title('Average Length vs Mg^{2+} Concentration')
%     errorbar(mgRange,avgMgLength,stdMgLength, 'o')
%     xlabel('Mg^{2+} concentration (mM)')
%     ylabel('Length (\mum)')
%     ylim([0 Inf])
%     %xticks(mgRange)
%     fig2pretty 
%     saveas(gcf, ['avgMgLength' channel '.fig'])
    
    %save(['dataPlots' channel '.mat'])