clear ALL
clc

%%%%%%%%%%%%%%%%%%%%%%%% Initial condition and time
y0=[220000 1 0 1 0 0.78947368421052631578947368421053 0.2];
time=[0 3000];

%% parameters
global lambda  h_b mu w gamma delta xi beta c
global r_p eta alpha h_p d_z k_p h_z d_b  %sigma beta_z
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
% sigma = 0.03;                   %  rate of bacteria-zooplankton association
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

%% equilibrium phytoplankton-zooplankton density
P=d_z*h_p/(eta*alpha-d_z);               %  equilibrium phytoplankton density
Z=(r_p/alpha)*(1-P/k_p)*(h_p+P);         %  equilibrium zooplankton density (Z = Z_B + Z_F)
  
i=0;
for sigma = [0 0.01 0.03 0.06]           %  rate of bacteria-zooplankton association
    i=i+1;
    j=0;
    for beta_z = [0.04 0.05 0.08]        %  transmission rate via zooplankton               
        j=j+1;

        R0_B=xi*N0/((d_b+c*sigma*Z/h_m)*(gamma+mu+delta))*(beta/h_b);
        R0_Z=xi*N0/((d_b+c*sigma*Z/h_m)*(gamma+mu+delta))*(beta_z*sigma*Z/(d_z*h_z*h_m));
        R0_BZ(i,j)=R0_B+R0_Z;

        %% ode run
        options=[];
        options = odeset('RelTol',1e-6, 'AbsTol',1e-6);
        [t,y]=ode45(@(t,y) (ode_plankton_cholera_sigma_betaz(t,y,sigma,beta_z)),time,y0,options);

        %%
        [I_max,I_max_time] = max(y(:,2));
        cum_cases(i,j)=y0(1)-y(end,1);
        peak_time(i,j)=t(I_max_time);
        peak(i,j)=I_max;

        dd=t(find(y(:,2)>1));
        if isempty(dd)
            duration(i,j)=0;
        else
            duration(i,j)=dd(end)-dd(1);
        end
    end
end

%% plot
sigma=[0 0.01 0.03 0.06];
beta_z=[0.04 0.05 0.08];

figure;
for i=1:3
    subplot(4,3,1)
    plot(sigma,[R0_BZ(:,i)]','-o','linewidth',1)
    hold on
    ylabel('R_{0_{BZ}}^{out}')
    ylim([1.22 2.65])
    xlim([-0.003 0.065])

    subplot(4,3,4)
    plot(sigma,[peak(:,i)]','-o','linewidth',1)
    hold on
    ylabel('Peak')
    ylim([3200 8000])
    xlim([-0.003 0.065])

    subplot(4,3,7)
    plot(sigma,[peak_time(:,i)]','-o','linewidth',1)
    hold on
    ylabel('Peak timing')
    ylim([150 520])
    xlim([-0.003 0.065])

    subplot(4,3,10)
    plot(sigma,[duration(:,i)]','-o','linewidth',1)
    hold on
    xlabel('\sigma')
    ylabel('Duration')
    ylim([450 1250])
    xlim([-0.003 0.065])

end


%% time series plots
i=0;
for beta_z=[0.08 0.05 0.04]
    i=i+1;
    for sigma=[0 0.01 0.03 0.06]

        options = odeset('RelTol',1e-6, 'AbsTol',1e-6);
        [t,y]=ode45(@(t,y) (ode_plankton_cholera_sigma_betaz(t,y,sigma,beta_z)),time,y0,options);


        subplot(3,3,[(3*i-1)+0.15,(3*i)-0.6])
        hold on
        plot(t,y(:,2),'linewidth',1)
        ylabel('Infected')
        xlabel('Time')
        xlim([0 600])
        ylim([0 8000])
    end
end