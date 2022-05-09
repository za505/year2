%NewtonMethodTemplate.m
%
%Zarina Akbary, 5/14/2020

%
%Purpose: Homework 1 for the Ninja Skills Club. To simulate a dynamic
%system

%empty m-file
clear, close all

%We're solving for dm/dt=a-b*m and
% dp/dt=(g*m)-(w*p)
%
%Step 1. Define input parameters
a=10; %rate of transcription in nM/min
b=1; %rate of RNA degradation in 1/min
g=20; %rate of translation in 1/min
w=0.1; %rate of proteolysis in 1/min

%Step 2. Set initial conditions.
m1=0; %initial mRNA concentration
p1=0.1; %initial protein concentration

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
%f=a-b*m AND h=(g*m)-(w*p)
%dm=f*dt AND dp=h*dt
for i=1:N-1
    f=a-b*m(i);
    dm=f*dt;
    m(i+1)=m(i)+dm;
    h=(g*m(i))-(w*p(i));
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

%% Part II, 07/14/2021
%Step 1. Define input parameters
b=10;
n=1;

%Step 2. Set initial conditions.
Ct1=40000;
Cu1=40000;

%Step 3. Set time increment, time steps, and time vector
dt=1; %time increment (1 second in minutes)
N=60*250; %time increments in 20 minutes
time=[0:N-1]*dt;

%Step 4. Pre-allocate vectors
Ct=zeros(1, N);
Cu=zeros(1, N);
Ct(1)=Ct1;
Cu(1)=Cu1;

%Step 5. Integrate
%dCu/dt = n(Ct-Cu) - b(Cu)
for i=1:N-1
    f=n*(Ct(i)-Cu(i))-b*Cu(i);
    dC=f*dt;
    Cu(i+1)=Cu(i)+dC;
end

%Step 6. Plot
figure
plot(time, Cu)
hold on
plot(time,Ct)
xlabel('Time (s)')
ylabel('Fluor. Concentration')
legend('unbleached', 'total')