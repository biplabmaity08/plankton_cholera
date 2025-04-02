function ydot =ode_plankton_cholera_wof(t,y)

global lambda beta h_b mu w gamma delta xi c

global r_p eta alpha h_p d_z k_p h_z d_b beta_z sigma p

global h_m d

S=y(1);
I=y(2);
R=y(3);
B=y(4);
Z_B=y(5);
Z_F=y(6);
P=y(7);

%periodic
% if t>0*360
k=@(t) k_p*(1+d*sin(2*pi*t/p));
% else
%  k=@(t) k_p;
% end

Sdot=lambda-beta*S*B/(h_b+B)-beta_z*S*Z_B/(h_z+Z_B)-mu*S+w*R;

Idot=beta*S*B/(h_b+B)+beta_z*S*Z_B/(h_z+Z_B)-(gamma+mu+delta)*I;

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