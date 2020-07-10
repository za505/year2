%PositiveAutoregulationTemplate.m
%
%Zarina Akbary, 5/19/2020
%
%Purpose: Homework 1 for the Ninja Skills Club. To simulate a dynamic
%system that exhibits positive autoregulation

%empty m-file
clear, close all

%We're solving for dm/dt=(a*(k+p))-(b*m) and
% dp/dt=(g*m)-(n*p)
%
%Step 1. Define input parameters
a=10; %rate of transcription in nM/min
b=1; %rate of RNA degradation in 1/min
g=20; %rate of translation in 1/min
n=0.1; %rate of proteolysis in 1/min
k=1000 %protein activation concentration?

%Step 2. Set initial conditions.
m1=0; %initial mRNA concentration
p1=0; %initial protein concentration

%Step 3. Set time increment, time steps, and time vector
dt=0.0167; %time increment (1 second in minutes)
N=60*250; %time increments in 20 minutes
time=[0:N-1]*dt;

%Step 4. Pre-allocate vectors
m=zeros(1, N);
p=zeros(1, N);
m(1)=m1;
p(1)=p1;

%Step 5. Integrate
%f=(a/k+p)-(b*m) AND h=(g*m)-(n*p)
%dm=f*dt AND dp=h*dt
for i=1:N-1
    f=(a*(k+p(i)))-(b*m(i));
    dm=f*dt;
    m(i+1)=m(i)+dm;
    h=(g*m(i))-(n*p(i));
    dp=h*dt;
    p(i+1)=p(i)+dp;
end

%Step 6. Plot
figure
plot(time, m)
hold on
plot(time,p)
xlabel('Time (min)')
ylabel('macromolecule concentration (nM)')
legend('mRNA', 'protein')