%BTfluo.m
%Rico Rojas, updated 1/21/19
%Calculates the average cytoplasmic fluorescence intensity from cell
%tracked with BacTrack.m.  

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
basename='051215_1';
channels={['/Users/Rico/Documents/MATLAB/Matlab Ready/' basename '/' basename '_2_a']};
recrunch=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==0;

curdir=cd;
for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

load([basename '_BT'])
load([basename '_BTlab'])

intensities=cell(length(channels),1);


for i=1:length(channels)
    cd(channels{i}); 
    intensities_temp=zeros(size(lcell));
    for t=1:T
        t
        imagename=fluo_directory{i}(t).name;
        im=imread(imagename);
        for j=1:ncells
            intensities_temp(j,t)=mean(im(pixels{j,t}));
            
        end
    end
    intensities_temp(intensities_temp==0)=NaN;
    icell{i}=intensities_temp;
end

icell_av=cell(length(channels),1);
for i=1:length(channels)
    icell_av{i}=nanmean(icell{i});
end

save([basename '_BTfluo'])

elseif recrunch==1
    load ([basename '_BTfluo'])
end

%Plot data
figure, hold on, 
for i=1:ncells
    plot(time,icell{1}(i,:))
end
plot(time,icell_av{1},'-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty

save([basename '_BTfluo'])