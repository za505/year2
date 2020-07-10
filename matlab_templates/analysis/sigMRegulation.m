%sigMRegulation.m
%
%Zarina Akbary, 6/4/2020
%
%Purpose: To simulate a dynamic system that exhibits sigM regulation

%empty m-file
clear, close all

%We're solving for (k1*s/k2+s)-(k3*s)
%
%Step 1. Define input parameters
k1=0.06; %upregulation of transcription, 1/min
k2=1.5; %basal transcription rate, nM/min
k3=0.0167; %degradation rate, 1/min

%Step 2. Set initial conditions.
s1=0.1; %initial sigM concentration, nM

%Step 3. Set time increment, time steps, and time vector
dt=0.0167; %time increment (1 second in minutes)
N=60*800; %time increments in 20 minutes
time=[0:N-1]*dt;

%Step 4. Pre-allocate vectors
s=zeros(1, N);
dS=zeros(1,N);
s(1)=s1;

%Step 5. Integrate
%f=((k1*s)/(k2+s))-(k3*s)
%dS=f*dt 
for i=1:N-1
    f=((k1*s(i))/(k2+s(i)))-(k3*s(i));
    dS(i)=f*dt;
    s(i+1)=s(i)+dS(i);
end


%Step 6. Plot
figure
plot(time, s)
xlabel('Time (min)')
ylabel('SigM concentration (nM)')

figure 
plot(s, dS)
xlabel('SigM concentration (nM)')
ylabel('Change in SigM concentration')


%%%%%%%%
%Let's confirm that we get no feedback if s1=0
%Set s(1)=1
s=zeros(1, N);
dS=zeros(1,N);

%Integrate as before
%f=((k1*s)/(k3+s))-(k3*s)
%dS=f*dt 
for i=1:N-1
    f=((k1*s(i))/(k2+s(i)))-(k3*s(i));
    dS(i)=f*dt;
    s(i+1)=s(i)+dS(i);
end

%{
%Plot as before
figure
plot(time, s)
xlabel('Time (min)')
ylabel('SigM concentration (nM) (S1=0)')

figure 
plot(s, dS)
xlabel('SigM concentration (nM)')
ylabel('Change in SigM concentration (S1=0)')
%%%%%%%
%}

%%%%%%%
%Now, let's plot the RHS as a function of S
%Set s(1)=1
s=[0:0.1:5];

%Calculate the RHS (right-hand side)
RHS=k1*s./(k2+s)-k3*s;

%Plot
figure,plot(s,RHS)
%%%%%%%

%%%%%%%
%Time to add an alpha term
%We're solving for a+(k1*s/k2+s)-(k3*s)
%
%Step 1. Define input parameters
k1=0.06; %upregulation of transcription, 1/min
k2=1.5; 
k3=0.0167; %degradation rate, 1/min
a=k1/10; %basal transcription rate, nM/min, sigA

%Step 2. Set initial conditions.
s1=0.1; %initial sigM concentration, nM

%Step 3. Set time increment, time steps, and time vector
dt=0.0167; %time increment (1 second in minutes)
N=60*800; %time increments in 20 minutes
time=[0:N-1]*dt;

%Step 4. Pre-allocate vectors
s=zeros(1, N);
dS=zeros(1,N);
s(1)=s1;

%Step 5. Integrate
%f=((k1*s)/(k3+s))-(k3*s)
%dS=f*dt 
for i=1:N-1
    f=a+((k1*s(i))/(k2+s(i)))-(k3*s(i));
    dS(i)=f*dt;
    s(i+1)=s(i)+dS(i);
end


%Step 6. Plot
figure
plot(time, s)
xlabel('Time (min)')
ylabel('SigM concentration (nM)')

figure 
plot(s, dS)
xlabel('SigM concentration (nM)')
ylabel('Change in SigM concentration')
%%%%%%%

%%%%%%%
%Let's see what happens if we add a hill co-efficient (which indicates
%cooperative binding)

%We're solving for a+(k1*s^n/k2+s^n)-(k3*s)
%
%Step 1. Define input parameters
k1=0.06; %upregulation of transcription, 1/min
k2=1.5; 
k3=0.0167; %degradation rate, 1/min
a=k1/10; %basal transcription rate, nM/min, sigA

%Step 2. Set initial conditions.
s1=0.1; %initial sigM concentration, nM

%Step 3. Set time increment, time steps, and time vector
dt=0.0167; %time increment (1 second in minutes)
N=60*3000; %time increments in 20 minutes
time=[0:N-1]*dt;

%Step 4. Pre-allocate vectors
n=zeros(1,5);
s=zeros(1, N);
dS=zeros(1,N);

%Step 5. Integrate
%f=((k1*s^n)/(k2+s^n))-(k3*s)
%dS=f*dt
for j=1:5
    n(j)=j;
    s(j,1)=s1;
for i=1:N-1
    f=a+((k1*s(j,i)^n(j))/(k2+s(j,i)^n(j)))-(k3*s(j,i));
    dS(j,i)=f*dt;
    s(j,i+1)=s(j,i)+dS(j,i);
end
end

%Step 6. Plot
figure, hold on
for j=1:5
plot(time, s(j, :))
end
xlabel('Time (min)')
ylabel('SigM concentration (nM)')

%{
figure,hold on
for j=1:5
plot(s, dS(j,:))
end
xlabel('SigM concentration (nM)')
ylabel('Change in SigM concentration')


%Now, let's plot the RHS as a function of n
n=[1:1:5];

%Calculate the RHS (right-hand side)
RHS=k1*s^n./(k2+s^n)-k3*s;

%Plot
figure,plot(n,RHS)
%}
%%%%%%%%%%%%

%%%%%%%
%Now, let's try modeling sigM binding and unbinding

%We're solving for sa=(k1*s*a)-(k2*sa)
%
%Step 1. Define input parameters
k1=0.024; %binding constant
k2=0.06; %unbinding constant

%Step 2. Set initial conditions.
sa1=0.01; %sigma factor-anti sigma factor complex
s1=0.05; %sigma factor concentration
a1=0.05; %anti-sigma factor concentration

%Step 3. Set time increment, time steps, and time vector
dt=0.0167; %time increment (1 second in minutes)
N=60*300; %time increments in 20 minutes
time=[0:N-1]*dt;

%Step 4. Pre-allocate vectors
sa=zeros(1, N);
s=zeros(1, N);
a=zeros(1, N);
dSA=zeros(1,N);
dS=zeros(1,N);
dA=zeros(1,N);
sa(1)=sa1;
s(1)=s1;
a(1)=a1;

%Step 5. Integrate
%f=(k1*s*a)-(k2*sa)
%dSA=f*dt
k3=0.06; %upregulation of transcription, 1/min
k4=1.5; 
k5=0.01;
k6=0.03;

for i=1:N-1
    f=(k1*s(i)*a(i))-(k2*sa(i));
    f2=(-k1*s(i)*a(i))+(k2*sa(i))+(k3*s/k4+s)-(k5*s);
    f3=(-k1*s(i)*a(i))+(k2*sa(i))+(k3*s/k4+s)-(k6*a);
    dSA(i)=f*dt;
    dS(i)=f2*dt;
    dA(i)=f3*dt;
    sa(i+1)=sa(i)+dSA(i);
    s(i+1)=s(i)+dS(i);
    a(i+1)=a(i)+dA(i);
end

%Step 6. Plot
figure
plot(time, sa)
xlabel('Time (min)')
ylabel('SigM Complex concentration (nM)')

figure 
plot(sa, dSA)
xlabel('SigM Complex concentration (nM)')
ylabel('Change in SigM Complex concentration')
%%%%%%%%%%%%

%%%%%%%

%Now, let's adding this new model to the previous one

%We're solving for b+(k1*s/k2+s)-(k3*s)-(k4*a)-(k5*s*a)+(k6*sa)
%
%Step 1. Define input parameters
k1=0.06; %upregulation of transcription, 1/min
k2=1.5; 
k3=0.05; %degradation rate, 1/min
k4=0.01; %degradation rate, 1/min
k5=0.024; %binding constant
k6=0.06; %unbinding constant
b=k1/10; %basal transcription rate, nM/min, sigA

%Step 2. Set initial conditions.
s1=0.01; %initial sigM concentration, nM
a1=0.01
%Step 3. Set time increment, time steps, and time vector
dt=0.0167; %time increment (1 second in minutes)
N=60*3000; %time increments in 20 minutes
time=[0:N-1]*dt;

%Step 4. Pre-allocate vectors
n=zeros(1,5);
s=zeros(1, N);
dS=zeros(1,N);

%Step 5. Integrate
%f=((k1*s^n)/(k2+s^n))-(k3*s)
%dS=f*dt
for j=1:5
    n(j)=j;
    s(j,1)=s1;
for i=1:N-1
    f=a+((k1*s(j,i)^n(j))/(k2+s(j,i)^n(j)))-(k3*s(j,i));
    dS(j,i)=f*dt;
    s(j,i+1)=s(j,i)+dS(j,i);
end
end

%Step 6. Plot
figure, hold on
for j=1:5
plot(time, s(j, :))
end
xlabel('Time (min)')
ylabel('SigM concentration (nM)')

%{
figure,hold on
for j=1:5
plot(s, dS(j,:))
end
xlabel('SigM concentration (nM)')
ylabel('Change in SigM concentration')
%}

%Now, let's plot the RHS as a function of n
n=[1:1:5];

%Calculate the RHS (right-hand side)
RHS=k1*s^n./(k2+s^n)-k3*s;

%Plot
figure,plot(n,RHS)
%%%%%%%%%%%%
%}