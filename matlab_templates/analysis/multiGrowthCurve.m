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
filename=('02232021_growthCurve');
dirname=('/Users/zarina/Downloads/NYU/Year2_2021_Spring/growthCurves/02232021/');
filepath=[dirname filename '.xlsx'];

xlRange='B51:EM110'; %what is the range of values for the OD data
timeRange='B49:EM49'; %what is the range of time values
wellRange='A51:A110'; %wells are we pulling data from

OD=xlsread(filepath, xlRange); %OD data
t=xlsread(filepath, timeRange); %time data
[data, wells]=xlsread(filepath, wellRange); %why does this need to be formatted differently?

T=length(t);
tscale=T/600; %why divide by 600? probably bc that's the # second in 10 min
time=([0:T-1]*600);
tmid=(time(2:2:end)+time(1:2:end-1))/2;
tim=tmid./3600;
time2=([0:T-1]*600)/3600;

%%%%%First, the LB data
%Now, let's define conditions
condlab={'ER438', 'ER439', 'ER440', 'ER441', 'ER002'};

%pre-allocate a structure
growth = struct; 

%Let's take an average of the blanks
blank1=[9:10:29, 49:10:59]; %LB blank
blank1_avg = mean(OD(blank1, :));

flag = 0;

for j=1:length(condlab)
    
    growth(j).name = condlab(j);
    
    index=[1:10:21];
    index=index + flag; 
    index

    growth(j).rep1 = OD(index(1), :);
    growth(j).rep2 = OD(index(2), :);
    growth(j).rep3 = OD(index(3), :);
    set = [growth(j).rep1; growth(j).rep2; growth(j).rep3];
    growth(j).ODavg = mean(set);
    
    growth(j).ODnorm = growth(j).ODavg - blank1_avg;
    growth(j).dOD=growth(j).ODnorm(:,2:2:end)-growth(j).ODnorm(:,1:2:end-1);
    growth(j).eOD=growth(j).dOD./(growth(j).ODnorm(:,1:2:end-1))/600;
    growth(j).eODsmooth=movingaverage(growth(j).eOD,10) * 3600; %what is moving average? the derivative?
    growth(j).max = max(growth(j).eODsmooth);

    flag = flag + 1;
end

%Now, let's plot it
figure, hold on
for j=1:length(condlab)
plot(tim,growth(j).eODsmooth)
end
xlabel('Time (h)')
ylabel('Growth Rate in LB (h^{-1})')
fig2pretty
legend(condlab)
cd(dirname)
saveas(gcf, [filename '_gr1'])
saveas(gcf, [filename '_gr1.png'])

figure, hold on
for j=1:length(condlab)  
plot(time2,growth(j).ODnorm)
end
xlabel('Time (h)')
ylabel('OD in LB (AU)')
fig2pretty
legend(condlab)
cd(dirname)
saveas(gcf, [filename '_OD1'])
saveas(gcf, [filename '_OD1.png'])

%%%%%Now, the LB + 1 M sorbitol data
%LB + 1 M sorbitol blank
blank2=[10:10:60]; 
blank2_avg = mean(OD(blank2, :));

%pre-allocate a structure
growthSorb = struct; 

flag = 0;

for j=1:length(condlab)
    
    growthSorb(j).name = condlab(j);
    
    index=[31:10:51];
    index=index + flag; 
    index

    growthSorb(j).rep1 = OD(index(1), :);
    growthSorb(j).rep2 = OD(index(2), :);
    growthSorb(j).rep3 = OD(index(3), :);
    set = [growthSorb(j).rep1; growthSorb(j).rep2; growthSorb(j).rep3];
    growthSorb(j).ODavg = mean(set);
    
    growthSorb(j).ODnorm = growthSorb(j).ODavg - blank2_avg;
    growthSorb(j).dOD=growthSorb(j).ODnorm(:,2:2:end)-growthSorb(j).ODnorm(:,1:2:end-1);
    growthSorb(j).eOD=growthSorb(j).dOD./(growthSorb(j).ODnorm(:,1:2:end-1))/600;
    growthSorb(j).eODsmooth=movingaverage(growthSorb(j).eOD,10) * 3600; %what is moving average? the derivative?
    growthSorb(j).max = max(growthSorb(j).eODsmooth);

    flag = flag + 1;
end

%Now, let's plot it
figure, hold on
for j=1:length(condlab)
plot(tim,growthSorb(j).eODsmooth)
end
xlabel('Time (h)')
ylabel('Growth Rate in LB + 1 M sorbitol (h^{-1})')
fig2pretty
legend(condlab)
cd(dirname)
saveas(gcf, [filename '_gr2'])
saveas(gcf, [filename '_gr2.png'])

figure, hold on
for j=1:length(condlab)  
plot(time2,growthSorb(j).ODnorm)
end
xlabel('Time (h)')
ylabel('OD in LB + 1 M sorbitol (AU)')
fig2pretty
legend(condlab)
cd(dirname)
saveas(gcf, [filename '_OD2'])
saveas(gcf, [filename '_OD2.png'])

%%%%%The sytox and FITC-dextran data
condlab2 = {'LB + 0.2 \muM sytox', 'LB + 0.02 \muM sytox', 'LB + 0.002 \muM sytox', 'LB + 100 \mug/mL FITC-dextran (10kDa)', 'LB + 10 \mug/mL FITC-dextran (10kDa)', 'LB + 1 \mug/mL FITC-dextran (10kDa)'};

%pre-allocate a structure
growthDye = struct; 

for j=1:length(condlab2)
    
    growthDye(j).name = condlab2(j);
    
    index=[6:8];
    index=index + 10*(j-1); 
    index

    growthDye(j).rep1 = OD(index(1), :);
    growthDye(j).rep2 = OD(index(2), :);
    growthDye(j).rep3 = OD(index(3), :);
    set = [growthDye(j).rep1; growthDye(j).rep2; growthDye(j).rep3];
    growthDye(j).ODavg = mean(set);
    
    growthDye(j).ODnorm = growthDye(j).ODavg - blank1_avg;
    growthDye(j).dOD=growthDye(j).ODnorm(:,2:2:end)-growthDye(j).ODnorm(:,1:2:end-1);
    growthDye(j).eOD=growthDye(j).dOD./(growthDye(j).ODnorm(:,1:2:end-1))/600;
    growthDye(j).eODsmooth=movingaverage(growthDye(j).eOD,10) * 3600; %what is moving average? the derivative?
    growthDye(j).max = max(growthDye(j).eODsmooth);

    flag = flag + 1;
end

%Now, let's plot it
figure, hold on
for j=1:length(condlab2)
plot(tim,growthDye(j).eODsmooth)
end
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
fig2pretty
legend(condlab2)
cd(dirname)
saveas(gcf, [filename '_gr3'])
saveas(gcf, [filename '_gr3.png'])

figure, hold on
for j=1:length(condlab2)  
plot(time2,growthDye(j).ODnorm)
end
xlabel('Time (h)')
ylabel('OD (AU)')
fig2pretty
legend(condlab2)
cd(dirname)
saveas(gcf, [filename '_OD3'])
saveas(gcf, [filename '_OD3.png'])

save('02232021_growthCurve.mat')

%%%%%%%%%%let's run a final check
for j=1:length(condlab2)
    
    growthDye(j).name = condlab2(j);
    
    index=[6:8];
    index=index + 10*(j-1); 
    index

    growthDye(j).rep1 = OD(index(1), :);
    growthDye(j).rep2 = OD(index(2), :);
    growthDye(j).rep3 = OD(index(3), :);
    set = [growthDye(j).rep1; growthDye(j).rep2; growthDye(j).rep3];
    growthDye(j).ODavg = mean(set);
    
    growthDye(j).ODnorm = growthDye(j).ODavg - blank2_avg;
    growthDye(j).dOD=growthDye(j).ODnorm(:,2:2:end)-growthDye(j).ODnorm(:,1:2:end-1);
    growthDye(j).eOD=growthDye(j).dOD./(growthDye(j).ODnorm(:,1:2:end-1))/600;
    growthDye(j).eODsmooth=movingaverage(growthDye(j).eOD,10) * 3600; %what is moving average? the derivative?
    growthDye(j).max = max(growthDye(j).eODsmooth);

    flag = flag + 1;
end

%Now, let's plot it
figure, hold on
for j=1:length(condlab2)
plot(tim,growthDye(j).eODsmooth)
end
xlabel('Time (h)')
ylabel('Growth Rate (h^{-1})')
fig2pretty
legend(condlab2)
cd(dirname)


figure, hold on
for j=1:length(condlab2)  
plot(time2,growthDye(j).ODnorm)
end
xlabel('Time (h)')
ylabel('OD (AU)')
fig2pretty
legend(condlab2)
cd(dirname)

