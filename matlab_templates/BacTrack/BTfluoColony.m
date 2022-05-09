%BTfluoColony.m
%Zarina Akbary, updated 03/22/21
%Calculates the pop. average cytoplasmic fluorescence intensity from cell
%tracked with BacTrack.m and plots it vs ion concentration
%based on BTfluo.m

clear, close all

%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save fluorescent image stacks
%directories by themselves. 

%INPUT
%basename: name to save the results to.
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%icell: Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a matrix (ncellxT) with the fluorescent intensities of each
        %cell, where rows are the cells and columns are time points.
%icell_av:  Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a vector containing the population-
        %average of the single-cell fluorescent intensities.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='03172021_Exp1_colony1';
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03172021_analysis/' basename];
channels={[filename '/' basename '_FSS/' basename '_full']};
%channels={[filename '/' basename '_FSS/' basename '_full']};
switch1=20; %frame during switch 1
switch2=31;
switch3=44;
switch4=57;
removeData=0;
recrunch=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==0;

%load data from BTfluo.m
load([filename '/' basename '_FSS/' basename '_BTFSS'])

%average intensity ratio
intensity_ratio(intensity_ratio==0)=NaN;
temp{1}=intensity_ratio;

avg_ratio=cell(length(channels),1);
for i=1:length(channels)
    avg_ratio{i}=nanmean(temp{i});
end

%input Mg+2 concentration
mg = zeros(1, width(intensity_ratio));
mg(1:switch2)=0;
mg(switch2+1:switch3)=6.66;
mg(switch3+1:switch4)=12.33;
mg(switch4+1:end)=20;

%Plot data
figure, hold on, 
plot(mg, avg_ratio{1}, '--r')
xlabel('Mg2+ (mM)')
ylabel('Avg Intensity Ratio (cell intensity/background)')
fig2pretty
xticks([0, 6.66, 12.33, 20])
% xline(60, '--k', '*PBS + 5% detergent')
% xline(114, '--k', '*PBS + FSS + FSS') %frame 19-30
% xline(234, '--k', '*PBS + FSS + FSS + 6.66 mM Mg2+') %frame 31-43
% xline(354, '--k', '*PBS + FSS + FSS + 12.33 mM Mg2+') %frame 44-56
% xline(474, '--k', '*PBS + FSS + FSS + 20 mM Mg2+') %frame 57-69
saveas(gcf, [filename '/' basename '_FSS/' basename,'_FSSratioAVG.png'])

figure, hold on, 
plot(mg ,icell_av{1},'-r')
xlabel('Mg2+ (mM)')
ylabel('Population Avg Intensity (A.U.)')
fig2pretty
xticks([0, 6.66, 12.33, 20])
% xline(60, '--k', '*PBS + 5% detergent')
% xline(114, '--k', '*PBS + FSS + FSS') %frame 19-30
% xline(234, '--k', '*PBS + FSS + FSS + 6.66 mM Mg2+') %frame 31-43
% xline(354, '--k', '*PBS + FSS + FSS + 12.33 mM Mg2+') %frame 44-56
% xline(474, '--k', '*PBS + FSS + FSS + 20 mM Mg2+') %frame 57-69
saveas(gcf, [filename '/' basename '_FSS/' basename,'_FSSAVG.png'])

save([filename '/' basename '_FSS/' basename '_BTFSSavg'])

end


