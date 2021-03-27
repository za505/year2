%analysisCheck.m
%Author: Zarina Akbary
%Date: 27 March 2021
%Purpose: go through multiple experiments and quality control the analysis

clear, close all

%Questions
%1Q. How many experiments are we re-analyzing?
%1A. Four. 03262021_Exp1, 03262021_Exp1, 03262021_Exp1, and 03262021_Exp1

%2Q. How many colonies are in each?
%2A. There are three in the first and four in the other three

%3Q. What do you want to do first?
%3A. Verify the accuracy of the cell tracking from BacTrack.m

%%%%Let's start with Exp1 from 03262021_analysis
maindir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03262021_analysis'];
basename=['03262021_Exp1_colony4'];
load([maindir '/' basename '/' basename '_phase/' basename '_figures/' basename '_BTphase'], 'dirname', 'directory', 'T', 'nc', 'boun');
cd([maindir '/' basename '/' basename '_phase/' basename '_erased']);

for t=1:T
    
    %Load image
    imagename=directory(t).name;
    im=imread(imagename);
    
    if t <= 5 | t >=T-5
        t
        figure
        imshow(im)
        hold on
        for k=1:nc(t)
        plot(boun{k,t}(:,1),boun{k,t}(:,2),'-r')
        end

        pause(0.3)
        close all
    end
end 