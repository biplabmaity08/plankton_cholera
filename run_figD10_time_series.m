clear ALL
clc

%%%%%%%%%%%%%%%%%%%%%%%% Initial condition and time
y0=[220000 1 0 1 0 0.78947368421052631578947368421053 0.2];
time=[0 1000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global lambda  h_b mu w gamma delta xi c
global r_p eta alpha h_p d_z k_p h_z d_b sigma beta_z beta
global h_m 

%%%%%%%%%%%%%%%%%%%%% PARAMETERS phyto-zooplankton
r_p = 0.5;                        %  phytoplankton intrinsic growth rate
k_p = 0.95;                       %  carrying capacity of phytoplankton
d_z = 0.06;                       %  death rate of zooplankton   
h_p = 0.6;                        %  half-saturation constant for phytoplankton                    
alpha = 0.4;                      %  maximum predation rate
eta = 0.6;                        %  conversion coefficient

%%%%%%%%%%%%%%%%%%%%%% B-Z association
c = 5*10^7;                       %  colonization coefficient of bacteria
sigma = 0.03;                     %  rate of bacteria-zooplankton association
h_m = 2*10^6;                     %  half-saturation constant of bacteria-zooplankton association

%%%%%%%%%%%%%%%%%%%%% Bacteria
d_b = 0.33;                       %  removal rate of bacteria

%%%%%%%%%%%%%%%%%%%%%% Human SIR
h_b = 1e9;                        %  half saturation constant of bacterial transmission
h_z = 20;                         %  half saturation constant of zooplankton-mediated transmission
mu = 0;   %3.8*10^(-5);           %  natutal death rate of human        
lambda = 0;    %0.025/365*(y0(1)+y0(2)+y0(3)); % Constant recruitment rate of human population
delta = 0.013;                    %  disease induced mortality rate of humans 
gamma = 1/5;                      %  recovery rate of infected human
w = 0;    %0.00092;               %  rate of immunity loss for recovered individuals
xi = 2000;                        %  bacteria shedding rate of infected human
N0 = y0(1)+y0(2)+y0(3);           %  initial human population size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=0.214;                              %  transmission rate via free-living bacteria 
%beta_z=0.1;                             %  transmission rate via zooplankton 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for beta_z=[0.04 0.06 0.08]
  
%% ode run
options=[];
[t,y]=ode45(@(t,y)(ode_plankton_cholera(t,y)),time,y0,options);

%% plot
hold on
for i=1:5

subplot(3,3,i)
hold on
plot(t,y(:,i),'linewidth',1.5);
hold on
xlabel('Time')
xlim([0 1000])
end

hold on
subplot(3,3,6)
hold on
yyaxis left
plot(t,y(:,6),'linewidth',1.5);
yyaxis right
plot(t,y(:,7),'linewidth',1.5);
hold on
xlabel('Time')
xlim([0 1000])

end