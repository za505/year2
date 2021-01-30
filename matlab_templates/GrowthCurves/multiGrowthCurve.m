%multiGrowthCurve.m
%Zarina Akbary, 12/22/2020
%Calculates and plots growth growth curves from all wells from Tecan M200.
%Manipulation of TecanGrowthCurve.m, InfiniteOriz.m, and GrowthCurveAnalysis.m
%
%INPUT
%dirname: full path to excel file with plate reader data 
%OUTPUT
%strain:strain rate of cells.

clear, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
filename=('01212021_growthCurve');
dirname=('/Users/zarina/Downloads/NYU/Year2_2021_Spring/growthCurves/01212021/');
filepath=[dirname filename '.xlsx'];

xlRange='B51:EM104'; %what is the range of values for the OD data
timeRange='B49:EM49'; %what is the range of time values
wellRange='A51:A104'; %wells are we pulling data from


OD=xlsread(filepath, xlRange); %OD data
t=xlsread(filepath, timeRange); %time data
[data, wells]=xlsread(filepath, wellRange); %why does this need to be formatted differently?

T=length(t);
tscale=T/600; %why divide by 600? probably bc that's the # second in 10 min
time=([0:T-1]*600);
tmid=(time(2:2:end)+time(1:2:end-1))/2;
tim=tmid./3600;
time2=([0:T-1]*600)/3600;

%Now, let's define conditions
condlab={'FM4-64', '647', 'Fluorescein', 'HADA', 'TADA', 'FDL', 'MV625', 'MV720', 'carboxyfluorescein', 'MV405', 'MVG', 'MTO', 'lysozyme A', 'lysozyme B', 'pronase', 'control'};

%pre-allocate a structure
growth = struct; 

%Let's take an average of the blanks
blank=[49:54];
blank_avg = mean(OD(blank, :));
    
for j=1:length(condlab)
    
    growth(j).name = condlab(j);
    
    index=[1:3];
    index = index + 3*(j-1); 
    index
    
    growth(j).rep1 = OD(index(1), :);
    growth(j).rep2 = OD(index(2), :);
    growth(j).rep3 = OD(index(3), :);
    set = [growth(j).rep1; growth(j).rep2; growth(j).rep3];
    growth(j).ODavg = mean(set);
    growth(j).ODnorm = growth(j).ODavg - blank_avg;
    growth(j).dOD=growth(j).ODnorm(:,2:2:end)-growth(j).ODnorm(:,1:2:end-1);
    growth(j).eOD=growth(j).dOD./(growth(j).ODnorm(:,1:2:end-1))/600;
    growth(j).eODsmooth=movingaverage(growth(j).eOD,10) * 3600; %what is moving average? the derivative?
    growth(j).max = max(growth(j).eODsmooth);
    
end

    
%Now, let's plot it
figure, hold on
for j=1:length(condlab)
plot(tim,growth(j).eODsmooth)
end
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
fig2pretty
legend(condlab)
cd(dirname)
saveas(gcf, [filename '_gr'])
saveas(gcf, [filename '_gr.png'])

figure, hold on
for j=1:length(condlab)  
plot(time2,growth(j).ODnorm)
end
xlabel('Time (h)')
ylabel('OD(AU)')
fig2pretty
legend(condlab)
cd(dirname)
saveas(gcf, [filename '_OD'])
saveas(gcf, [filename '_OD.png'])

