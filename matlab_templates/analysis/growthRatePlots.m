%growthRatePlots
%modified from FCTCompile
%Average growth rate and colony size from several runs
%Zarina Akbary and Paola Bardetti
%07/10/2020

%Compiles data from several positions tracked with FluoColonyTrack 
%and overlaps them on the same plot

clear, close all

%%%%%yhdL downshift
%list the files you're compiling data from
folder=['/Users/zarina/Downloads/NYU/Lab_2020_Summer/06162020_nutrientShiftAssay'];
basename='06162020_yhdL_nutrientShift';
filelist={[basename '_colony1'], [basename '_colony2'], [basename '_colony3'], [basename '_colony4']};
smooth=20; 
max=241; %max length of vector

colonylab = {'Colony 1 (S1P1)', 'Colony 2 (S2P1)', 'Colony 3 (S2P2)', 'Colony 4 (S3P1)', 'Average'};

%Extract the variables of interest from the FCT files
for i=1:length(filelist)
    cd([folder '/' filelist{i}])
    load([filelist{i},'_FCT'],'ColonyArea','eAsmooth','Ncells','tmid')
    data{1,i}=ColonyArea;
    data{2,i}=eAsmooth;
    data{3,i}=Ncells;
    %data{4,i}=tmid;
    %data{4,i}=time;
    %data{5,i}=tscale; 
end

%Assure that tmid will match the max vector length
for i=1:max
    tmid(i) = (30+(60*i))/3600;
end

%make sure the eA and Ncells matrices are the same length
for i=1:length(filelist)
    if length(data{2,i}(1,:)) < max
        col= max - length(data{2,i}(1,:));
        data{2,i}(1, max-col+1:max)=NaN(1,col);
        data{3,i}(1, max-col+1:max)=NaN(1,col);
    else
        data{3,i}=data{3,i}(1, 1:max);
    end
end

%Preallocate a matrix for the growth rate and colony size averages
smoothAvg=zeros(1,max);
NcellAvg=zeros(1,max);

%Adjust all eAsmooth values
for i=1:length(filelist)
    for j=1:max
        data{2,i}(1,j)=data{2,i}(1,j)*3600;
    end
end

for j=1:max
nSmooth=[data{2,1}(1,j), data{2,2}(1,j), data{2,3}(1,j), data{2,4}(1,j)];
nCells=[data{3,1}(1,j), data{3,2}(1,j), data{3,3}(1,j), data{3,4}(1,j)];
smoothAvg(j)=nanmean(nSmooth);
NcellAvg(j)=nanmean(nCells);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
colors = [0.85 0.85 0.85;0.75 0.75 0.75; 0.7 0.7 0.7;0.4 0.4 0.4];

%Plot average growth rate
% figure
% plot(tmid, smoothAvg, '--r','LineWidth', 4)
% xlabel('Time (h)')
% ylabel('Average Growth Rate (s^{-1})')
% fig2pretty

%Plot all growth rates 
figure, hold on
for i=1:length(filelist)
    plot(tmid, data{2,i}(1,:),'Color',colors(i, :))
end
plot(tmid, smoothAvg, '--r','LineWidth', 4)
xlabel('Time (h)')
ylabel('Growth Rate (s^{-1})')
legend(colonylab)
fig2pretty
cd(folder)
savefig([basename '_overlay' '.fig'])

%Plot average colony size
figure
plot(tmid, NcellAvg/(NcellAvg(1)), '--r','LineWidth', 4)
xlabel('Time (h)')
ylabel('Average Colony Size (# of cells)')
fig2pretty

%Plot all colony size curves
figure, hold on
for i=1:length(filelist)
    plot(tmid, data{3,i}(1,:)/NcellAvg(1), 'Color',colors(i, :))
end
xlabel('Time (h)')
ylabel('Colony Size (# of cells)')
legend(colonylab)
fig2pretty
cd(folder)
savefig([basename '_overlay' '.fig'])