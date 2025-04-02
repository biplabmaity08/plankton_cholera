clear ALL
clc

%%%%%%%%%%%%%%%%%%%%%%%% Initial condition and time
y0=[220000 1 0 1 0 0.78947368421052631578947368421053 0.2];
time=[0 5000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global lambda  h_b mu w gamma delta xi c
global r_p eta alpha h_p d_z k_p h_z d_b sigma % beta beta_z
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
sigma = 0;                   %  rate of bacteria-zooplankton association
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

for beta_z=0                      % transmission rate from zooplankton to human

    n=length(0.04:0.01:0.15);
    R0_BZ=zeros(1,n); %[];
    cum_cases_peak=zeros(1,n); %[];          %  cumulative_cases_at_peak
    cum_cases=zeros(1,n); %[];               %  cumulative infection at the end of outbreak
    overshoot=zeros(1,n); %[];               %  overshoot
    overshoot_by_cum_cases=zeros(1,n); %[];  %  ratio of overshoot and cum_cases
    peak_time=zeros(1,n); %[];                   %  peak timing
    peak=zeros(1,n); %[];                        %  peak value
    duration=zeros(1,n); %[];                    %  epidemic duration

    i=0;
    for beta = 0.214:0.02:0.6
        i=i+1;
        %% R0_BZ
        P=d_z*h_p/(eta*alpha-d_z);               %  equilibrium phytoplankton density
        Z=(r_p/alpha)*(1-P/k_p)*(h_p+P);         %  equilibrium zooplankton density (Z = Z_B + Z_F)

        R0_B=xi*N0/((d_b+c*sigma*Z/h_m)*(gamma+mu+delta))*(beta/h_b);
        R0_Z=xi*N0/((d_b+c*sigma*Z/h_m)*(gamma+mu+delta))*(beta_z*sigma*Z/(d_z*h_z*h_m));
        R0_BZ(i)=R0_B+R0_Z;

        if R0_BZ(i)<=1
            cum_cases_peak(i)=NaN;
            cum_cases(i)=0;
            overshoot(i)=NaN;
            overshoot_by_cum_cases(i)=NaN;
            peak_time(i)=0;
            peak(i)=y0(2);
            duration(i)=NaN;

        else
            %%%%%%%%%%%%%%%%%
            options=[];
            options = odeset('RelTol',1e-8, 'AbsTol',1e-8 );
            [t,y]=ode23s(@(t,y) ode_plankton_cholera_overshoot(t,y,beta,beta_z),time,y0,options);

            %%
            [I_max,I_max_time]=max(y(:,2));
            cum_cases_peak(i)=y0(1)-y(I_max_time,1);
            cum_cases(i)=y0(1)-y(end,1);
            overshoot(i)=(cum_cases(i)-cum_cases_peak(i));
            overshoot_by_cum_cases(i)=overshoot(i)/cum_cases(i);
            peak_time(i)=t(I_max_time);
            peak(i)=I_max;

            dd=t(find(y(:,2)>1));
            if isempty(dd)
                duration(i)=0;
            else
                duration(i)=dd(end)-dd(1);
            end

        end
    end
    beta=0.214:0.02:0.6;

    figure;

    hold on
    subplot(2,2,1)
    hold on
    yyaxis left
    plot(beta,R0_BZ,'linewidth',2)
    ylabel('R_{0_{BZ}}^{out}')
    ylim([1 4.2])

    yyaxis right
    plot(beta,cum_cases/N0,'linewidth',2)
    ylabel('Attack rate')
    ylim([0.37 1])
    xlabel('\beta')
    xlim([0.2 0.614])
    axis square

    hold on
    subplot(2,2,2)
    hold on
    yyaxis left
    plot(beta,cum_cases_peak/N0,'linewidth',2)
    ylabel('Cum inf at peak')
    ylim([0.2 0.65])
    yyaxis right
    plot(beta,overshoot/N0,'linewidth',2)
    ylabel('Overshoot')
    ylim([0.2 0.65])
    xlabel('\beta')
    xlim([0.2 0.614])
    axis square

    subplot(2,2,3)
    hold on
    plot(beta,overshoot_by_cum_cases,'linewidth',2)
    ylabel('Overshoot/Attack rate ( \rho_{\rm OA} )')
    xlabel('\beta')
    ylim([0.32 0.5])
    xlim([0.2 0.614])
    axis square

    subplot(2,2,4)
    hold on
    yyaxis left
    plot(beta,peak_time,'linewidth',2)
    ylabel('Peak time')
    ylim([0 600])
    yyaxis right
    plot(beta,duration,'linewidth',2)
    ylabel('Duration')
    xlabel('\beta')
    ylim([0 600])
    xlim([0.2 0.614])
    axis square


end