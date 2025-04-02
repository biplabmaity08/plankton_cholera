clear ALL
clc

%%%%%%%%%%%%%%%%%%%%%%%% Initial condition and time
global time1
y01=[220000 50 0 10^6 0 0.78947368421052631578947368421053 0.2];
time1=[0 250*360];
time2=[time1(end) 500*360];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global lambda  h_b mu w gamma delta xi c
global eta alpha h_p d_z k_p h_z d_b r_p d sigma beta_z beta p
global  h_m

%%%%%%%%%%%%%%%%%%%%% PARAMETERS phyto-zooplankton
r_p=0.5;
d=0.8; 
p=360/2;       
k_p=0.27;             
d_z=0.06;      
h_p=0.6;                        
alpha=0.4;  
eta=0.6;    

%%%%%%%%%%%%%%%%%%%%%% B-Z association
c=5*10^7;
sigma=0.05;
h_m=2*10^6;

%%%%%%%%%%%%%%%%%%%%% Bacteria
d_b=0.33;             

%%%%%%%%%%%%%%%%%%%%%% Human SIR
h_b=1e9;           
h_z=20;   
mu=3.8*10^(-5);        
lambda=15.0719; 
delta=0.013;           
gamma=1/5;               
w=0.00092; 
xi=1000; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_z=0.15; 
beta=0.15;
f_e=0.8;                             % filtration efficacy

% t_i=0;                             % initiated at T_1 (following P abundance)
% t_i=41;                            % initiated at T_2 (best scenario / following Z abundance)
% t_i=70;                            % initiated at T_3 (following I increase)
% t_i=130;                           % initiated at T_4 (worst scenatio) 

figure;
j=0;
for t_i = [130 0 70 41]
    j=j+1;

    %% ode_run

    options=[];
    options = odeset('RelTol',1e-6, 'AbsTol',1e-6 );
    [t,y]=ode45('ode_plankton_cholera_wof',time1,y01,options);

    time_span1=time1(1):time1(end);
    y1=interp1(t,y,time_span1,'spline');

    %%
    y02=[y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7)];

    [t,y]=ode45(@(t,y) ode_plankton_cholera_wf(t,y,t_i,f_e),time2,y02,options);
    time_span2=time2(1):time2(end);
    y2=interp1(t,y,time_span2,'spline');

    time_span3=[time_span1(1:end-1) time_span2];
    y3=[y1(1:end-1,:); y2];

    %% control profile
    C=[];
    for i=1:length(time_span3)
        if i>time1(end)+t_i && rem(floor((i-time1(end)-t_i)/90),2)==0
            C(i)= 1-f_e;
        else
            C(i)= 1;
        end
    end

    %% with and without filtration control function C
    time_wof=time_span1(end-2*p):time_span1(end);
    time_wf=time_span3(end-2*p):time_span3(end);

    C_wof=C(time_wof);
    C_wf=C(time_wf);

    %% new & cumulative cases || with and without filtration
    y0_wof=[y1(time_wof(1),:)];
    y0_wf=[y3(time_wf(1),:)];

    S_wof=y1(time_wof,1); I_wof=y1(time_wof,2); B_wof=y1(time_wof,4); Zb_wof=y1(time_wof,5); Z_wof=y1(time_wof,6); P_wof=y1(time_wof,7);
    S_wf=y3(time_wf,1); I_wf=y3(time_wf,2); B_wf=y3(time_wf,4); Zb_wf=y3(time_wf,5); Z_wf=y3(time_wf,6); P_wf=y3(time_wf,7);

    NC_wof=beta.*S_wof.*B_wof./(h_b+B_wof) + beta_z*S_wof.*C_wof'.*Zb_wof./(h_z+C_wof'.*Zb_wof);
    NC_wf=beta.*S_wf.*B_wf./(h_b+B_wf) + beta_z*S_wf.*C_wf'.*Zb_wf./(h_z+C_wf'.*Zb_wf);

    L=length(time_wof);

    New_cases_wof=zeros(L,1);
    New_cases_wf=zeros(L,1);
    Cum_cases_wof=zeros(L,1);
    Cum_cases_wf=zeros(L,1);
    New_cases_wof(1)=NC_wof(1);
    New_cases_wf(1)=NC_wf(1);
    Cum_cases_wof(1)=New_cases_wof(1);
    Cum_cases_wf(1)=New_cases_wf(1);

    for i=1:(L-1)
        New_cases_wof(i+1)=trapz(NC_wof(i:i+1));
        New_cases_wf(i+1)=trapz(NC_wf(i:i+1));
        Cum_cases_wof(i+1)=Cum_cases_wof(i)+New_cases_wof(i+1);
        Cum_cases_wf(i+1)=Cum_cases_wf(i)+New_cases_wf(i+1);

    end

    redu_cum_cases=(Cum_cases_wof(end)-Cum_cases_wf(end))*100/Cum_cases_wof(end);

    %% plots

    tp=1:length(time_wof);
    hold on
    subplot(3,4,[11-2*(j-1) 12-2*(j-1)])
    yyaxis left
    plot(tp,I_wof,'-.',tp,I_wf,'-','linewidth',2);
    ylabel('Infected')
    ylim([0 900])

    yyaxis right
    plot(tp,Cum_cases_wof,'-.',tp,Cum_cases_wf,'-','linewidth',2);
    ylabel('Cumulative Infections')
    ylim([0 25000])
    xlabel('Time')
    xlim([0 360])

end

hold on
subplot(3,4,[3 4])
yyaxis left
plot(tp,P_wf,'-',tp,Z_wf+Zb_wf,'--','linewidth',2);
ylabel('P, Z')
ylim([0 0.55])

yyaxis right
plot(tp,I_wof,'linewidth',2);
ylabel('Infected')
ylim([0 900])
xlabel('Time')
xlim([0 360])