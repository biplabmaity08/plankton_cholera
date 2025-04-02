function ydot =ode_plankton_cholera_wf(t,y,t_i,f_e)

global lambda beta h_b mu w gamma delta xi c

global r_p eta alpha h_p d_z k_p h_z d_b beta_z sigma p

global h_m d

global time1

S=y(1);
I=y(2);
R=y(3);
B=y(4);
Z_B=y(5);
Z_F=y(6);
P=y(7);

%periodic
k=@(t) k_p*(1+d*sin(2*pi*t/p));

% filtration on Zb
if t>time1(end)+t_i && rem(floor((t-time1(end)-t_i)/90),2)==0
    C=@(t) 1-f_e;
else
    C=@(t) 1;
end

Sdot=lambda-beta*S*B/(h_b+B)-beta_z*S*C(t)*Z_B/(h_z+C(t)*Z_B)-mu*S+w*R;

Idot=beta*S*B/(h_b+B)+beta_z*S*C(t)*Z_B/(h_z+C(t)*Z_B)-(gamma+mu+delta)*I;

Rdot=gamma*I-(mu+w)*R;

Bdot=xi*I-d_b*B-c*sigma*B*Z_F/(h_m+B);
%bacteria

Z_Bdot=sigma*B*Z_F/(h_m+B)-d_z*Z_B;
%zooplankton associated with bacteria

Z_Fdot=eta*alpha*P*(Z_F+Z_B)/(h_p+P)-d_z*Z_F-sigma*B*Z_F/(h_m+B);
% zooplankton

Pdot=r_p*P*(1-P/k(t))-alpha*P*(Z_F+Z_B)/(h_p+P);
% phytoplankton
 

ydot=[Sdot;Idot;Rdot;Bdot;Z_Bdot;Z_Fdot;Pdot];
end