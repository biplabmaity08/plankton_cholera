clear ALL
clc

%%%%%%%%%%%%%%%%%%%%%%%% Initial condition and time
y0=[220000 1 0 1 0 0.78947368421052631578947368421053 0.2];
time=[0 1500];
time_span=0:time(end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global lambda  h_b mu w gamma delta xi c
global r_p eta alpha h_p d_z k_p h_z d_b  %sigma beta_z beta
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
beta_z=c*beta*d_z*h_z/(h_b*d_b);         %  transmission rate via zooplankton               

%% equilibrium phytoplankton-zooplankton density
P=d_z*h_p/(eta*alpha-d_z);               %  equilibrium phytoplankton density
Z=(r_p/alpha)*(1-P/k_p)*(h_p+P);         %  equilibrium zooplankton density (Z = Z_B + Z_F)

X=xi*N0/(gamma+delta+mu);

R_0=X*beta/(d_b*h_b);                    %  R_0_BZ in absence of B-Z association (sigma = 0)

%% 
figure;
i=0;
for R_0B=[1*R_0 0.8*R_0 0.6*R_0 0.4*R_0]      %  similar as R_0Z=[0*R_0 0.2*R_0 0.4*R_0 0.6*R_0]

    i=i+1;

    R_0Z(i)=R_0-R_0B;
    sigma(i)=d_b*h_m/(c*Z)*(R_0/R_0B-1);      %  rate of B-Z association that captures relative contribution

    % options=[];
    options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
    [t,y]=ode45(@(t,y) ode_plankton_cholera_rel_contb(t,y,beta,beta_z,sigma(i)),time,y0,options);

    %% peak infections, peak timing, cumulative infection, epidemic duration
    [I_max,I_max_time] = max(y(:,2));
    cum_cases_peak(i) = y0(1)-y(I_max_time,1);           %  cumulative_cases_at_peak
    cum_cases(i) = y0(1)-y(end,1);                       %  cumulative infection at the end of outbreak
    peak_time(i) = t(I_max_time);                        %  peak timing
    peak(i) = I_max;                                     %  peak value

    dd = t(find(y(:,2)>1));
    if isempty(dd)
        duration(i) = NaN;
    else
        duration(i)=dd(end)-dd(1);                        %  duration
    end

    %% new cases
    ymodel = new_cases(time_span,y0,beta,beta_z,sigma(i));

    cum_case_B(i) = ymodel(end,4);
    cum_case_ZB(i) = ymodel(end,5);

    %% EGR
    M = [-(gamma+delta) N0*beta/h_b N0*beta_z/h_z;
        xi -(d_b+c*sigma(i)*Z/h_m) 0;
        0 sigma(i)*Z/h_m -d_z];

    Eig_real_M = real(eig(M));
    EGR_M(i) = max(Eig_real_M);

    %% plot

    hold on
    subplot(2,2,1);
    hold on
    plot(t,y(:,2),'linewidth',2);
    ylabel('Infected');
    xlabel('Time');
    xlim([70 700]);
    ylim([0 4700]);


    subplot(2,2,3)
    yyaxis left
    hold on
    plot(time_span,ymodel(:,1),'linewidth',2)
    ylabel('- New Infections (B)');
    ylim([0 1000]);

    yyaxis right
    hold on
    plot(time_span,ymodel(:,2),'linewidth',2)
    ylabel('-- New Infections (Z_B)');
    ylim([0 150]);
    xlabel('Time');
    xlim([70 700]);

end

subplot(2,2,4);
i=length(duration);
yyaxis left
plot(1:i,duration,'-*','linewidth',2)
ylabel('Duration');
ylim([450 900]);
yyaxis right
plot(1:i,EGR_M,'-o','linewidth',2)
ylabel('EGR');
ylim([0.015 0.05]);
xlim([0.8 4.2]);

subplot(2,2,2);
hold on
br=[(cum_case_B/N0)' (cum_case_ZB/N0)' (cum_cases_peak/N0)'];

hold on
yyaxis left
bar(0:50:150,br(:,1:2),'stacked')
ylabel('Cumulative Infections');
ylim([0.35 0.462])

hold on
yyaxis right
bar(20:50:170,br(:,3))
ylabel('Cum Inf at Peak');
ylim([0.22 0.25]);
xlim([-20 188]);
