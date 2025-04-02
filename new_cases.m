function y=new_cases(time,y0,beta,beta_z,sigma)
global lambda  h_b mu w gamma delta xi c
global r_p eta alpha h_p d_z k_p h_z d_b  %sigma beta_z beta
global h_m 

T=length(time);

[t,y2] = ode45(@(t,y)ode_plankton_cholera_rel_contb(t,y,beta,beta_z,sigma),0:0.1:T-1,y0,[]);
y1=interp1(t,y2,time,'spline');

S=y1(:,1); I=y1(:,2); R=y1(:,3); B=y1(:,4); Z_B=y1(:,5); Z_F=y1(:,6); P=y1(:,7);

NC_B = beta*S.*B./(h_b+B);
NC_ZB = beta_z*S.*Z_B./(h_z+Z_B);
NC = beta*S.*B./(h_b+B)+beta_z*S.*Z_B./(h_z+Z_B);

L=length(time);
New_case_B = zeros(L,1);
New_case_ZB = zeros(L,1);
New_cases = zeros(L,1);
Cumulative_cases_B = zeros(L,1);
Cumulative_cases_ZB = zeros(L,1);
Cumulative_cases = zeros(L,1);

New_case_B(1) = I(1);
New_case_ZB(1) = 0;
New_cases(1) = I(1); 
Cumulative_cases_B(1) = New_case_B(1);
Cumulative_cases_ZB(1) = New_case_ZB(1);
Cumulative_cases(1) = New_cases(1);

for i=1:(L-1)
    New_case_B(i+1) = trapz(NC_B(i:i+1));
    New_case_ZB(i+1) = trapz(NC_ZB(i:i+1));
    New_cases(i+1) = trapz(NC(i:i+1));
    Cumulative_cases_B(i+1) = Cumulative_cases_B(i) + New_case_B(i+1);
    Cumulative_cases_ZB(i+1) = Cumulative_cases_ZB(i) + New_case_ZB(i+1);
    Cumulative_cases(i+1) = Cumulative_cases(i) + New_cases(i+1);
    
 end
 
 y=[New_case_B New_case_ZB New_cases Cumulative_cases_B Cumulative_cases_ZB Cumulative_cases];
