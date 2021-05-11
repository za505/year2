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
filename=('05082021_growthCurve');
dirname=('/Users/zarina/Downloads/NYU/Year2_2021_Spring/growthCurves/05082021/');
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

%%%%%First, the 7H9 data
%Now, let's define conditions
condlab={'LB', 'LB + 1 M sorbitol', 'LB + 20 mM Mg^{2+}'};
strainlab={'ER048', 'ER419', 'ER420', 'ER435'};

%pre-allocate a structure
growth = struct; 

%Let's take an average of the blanks
blank1=[37:39, 46:48]; %LB blank
blank1_avg = mean(OD(blank1, :));

%flag = 0;

for j=1:length(strainlab)
    
    growth(j).name = strainlab(j);
    
    index=[1:3];
    index=index + 9 * (j-1); 
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

    %flag = flag + 1;
end

%Now, let's plot it
figure, hold on
for j=1:length(strainlab)
plot(tim,growth(j).eODsmooth)
end
xlabel('Time (h)')
ylabel('Growth Rate in LB (h^{-1})')
fig2pretty
legend(strainlab)
cd(dirname)
saveas(gcf, [filename '_gr1'])
saveas(gcf, [filename '_gr1.png'])

figure, hold on
for j=1:length(strainlab)  
plot(time2,growth(j).ODnorm)
end
xlabel('Time (h)')
ylabel('OD in LB (AU)')
fig2pretty
legend(strainlab)
cd(dirname)
saveas(gcf, [filename '_OD1'])
saveas(gcf, [filename '_OD1.png'])


%%%%%Now, the LB + 1 M sorbitol data
%Let's take an average of the blanks
blank2=[40:42, 49:51]; %LB + 1 M sorbitol blank
blank2_avg = mean(OD(blank2, :));

%flag = 0;

for j=1:length(strainlab)
    
    growth(j+3).name = strainlab(j);
    
    index=[4:6];
    index=index + 9 * (j-1); 
    index

    growth(j+3).rep1 = OD(index(1), :);
    growth(j+3).rep2 = OD(index(2), :);
    growth(j+3).rep3 = OD(index(3), :);
    set = [growth(j+3).rep1; growth(j+3).rep2; growth(j+3).rep3];
    growth(j+3).ODavg = mean(set);
    
    growth(j+3).ODnorm = growth(j+3).ODavg - blank2_avg;
    growth(j+3).dOD=growth(j+3).ODnorm(:,2:2:end)-growth(j+3).ODnorm(:,1:2:end-1);
    growth(j+3).eOD=growth(j+3).dOD./(growth(j+3).ODnorm(:,1:2:end-1))/600;
    growth(j+3).eODsmooth=movingaverage(growth(j+3).eOD,10) * 3600; %what is moving average? the derivative?
    growth(j+3).max = max(growth(j+3).eODsmooth);

    %flag = flag + 1;
end

%Now, let's plot it
figure, hold on
for j=1:length(strainlab)
plot(tim,growth(j+3).eODsmooth)
end
xlabel('Time (h)')
ylabel('Growth Rate in LB + 1 M sorbitol (h^{-1})')
fig2pretty
legend(strainlab)
cd(dirname)
saveas(gcf, [filename '_gr2'])
saveas(gcf, [filename '_gr2.png'])

figure, hold on
for j=1:length(strainlab)  
plot(time2,growth(j+3).ODnorm)
end
xlabel('Time (h)')
ylabel('OD in LB + 1 M sorbitol (AU)')
fig2pretty
legend(strainlab)
cd(dirname)
saveas(gcf, [filename '_OD2'])
saveas(gcf, [filename '_OD2.png'])

%%%%%Now, the LB + 20 mM Mg^{2+} data
%LB + 20 mM Mg^{2+} blank
blank3=[43:45, 52:54]; 
blank3_avg = mean(OD(blank3, :));

%pre-allocate a structure
%growthSorb = struct; 

%flag = 0;

for j=1:length(strainlab)
    
    growth(j+6).name = strainlab(j);
    
    index=[7:9];
    index=index + 9 * (j-1); 
    index

    growth(j+6).rep1 = OD(index(1), :);
    growth(j+6).rep2 = OD(index(2), :);
    growth(j+6).rep3 = OD(index(3), :);
    set = [growth(j+6).rep1; growth(j+6).rep2; growth(j+6).rep3];
    growth(j+6).ODavg = mean(set);
    
    growth(j+6).ODnorm = growth(j+6).ODavg - blank3_avg;
    growth(j+6).dOD=growth(j+6).ODnorm(:,2:2:end)-growth(j+6).ODnorm(:,1:2:end-1);
    growth(j+6).eOD=growth(j+6).dOD./(growth(j+6).ODnorm(:,1:2:end-1))/600;
    growth(j+6).eODsmooth=movingaverage(growth(j+6).eOD,10) * 3600; %what is moving average? the derivative?
    growth(j+6).max = max(growth(j+6).eODsmooth);

    %flag = flag + 1;
end

%Now, let's plot it
figure, hold on
for j=1:length(strainlab)
plot(tim,growth(j+6).eODsmooth)
end
xlabel('Time (h)')
ylabel('Growth Rate in LB + 20 mM Mg^{2+} (h^{-1})')
fig2pretty
legend(strainlab)
cd(dirname)
saveas(gcf, [filename '_gr3'])
saveas(gcf, [filename '_gr3.png'])

figure, hold on
for j=1:length(strainlab)  
plot(time2,growth(j+6).ODnorm)
end
xlabel('Time (h)')
ylabel('OD in LB + 20 mM Mg^{2+} (AU)')
fig2pretty
legend(strainlab)
cd(dirname)
saveas(gcf, [filename '_OD3'])
saveas(gcf, [filename '_OD3.png'])

