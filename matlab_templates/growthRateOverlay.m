%growthRateOverlay.m
%
%Zarina Akbary
%01/20/2020
%
%Plot growth rate data generated by multiple FluoColonyTrack.m runs to one
%figure
%Save figure to the main folder

clear, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder=['C:\Users\zarin\Downloads\NYU\R2_Rojas Lab\growthRate_01232020_yhdLdownshift'];
colonylab = {'Colony 1 (S1P1)', 'Colony 2 (S1P2)', 'Colony 3 (S2P1)', 'Colony 4 (S4P1)'};

%Add path and load FCT data from each figure of interest
addpath(genpath('C:\Users\zarin\Downloads\NYU\R2_Rojas Lab\growthRate_01232020_yhdLdownshift\01232020_Series1_P1'));
fig1= load(['01232020_Series1_P1_FCT.mat'], 'tmid', 'eAsmooth')

addpath(genpath('C:\Users\zarin\Downloads\NYU\R2_Rojas Lab\growthRate_01232020_yhdLdownshift\01232020_Series1_P2'));
fig2=load(['01232020_Series1_P2_FCT.mat'], 'tmid', 'eAsmooth')

addpath(genpath('C:\Users\zarin\Downloads\NYU\R2_Rojas Lab\growthRate_01232020_yhdLdownshift\01232020_Series2_P1'));
fig3=load(['01232020_Series2_P1_FCT.mat'], 'tmid', 'eAsmooth')

addpath(genpath('C:\Users\zarin\Downloads\NYU\R2_Rojas Lab\growthRate_01232020_yhdLdownshift\01232020_Series4_P1'));
fig4=load(['01232020_Series4_P1_FCT.mat'], 'tmid', 'eAsmooth')

%Generate figure
figure
hold on
plot(fig1.tmid/3600,fig1.eAsmooth*3600)
plot(fig2.tmid/3600,fig2.eAsmooth*3600)
plot(fig3.tmid/3600,fig3.eAsmooth*3600)
plot(fig4.tmid/3600,fig4.eAsmooth*3600)
xlabel('Time (h)')
ylabel('Growth Rate (s^{-1})')
legend(colonylab)
fig2pretty
cd(folder)
savefig(['01232020_yhdLdownshift_growthRate.fig'])
