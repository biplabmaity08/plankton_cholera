function ydot =ode_plankton_cholera_rel_contb(t,y,beta,beta_z,sigma)

global lambda h_b mu w gamma delta xi c

global r_p eta alpha h_p d_z k_p h_z d_b %sigma beta_z beta

global h_m

S=y(1);
I=y(2);
R=y(3);
B=y(4);
Z_B=y(5);
Z_F=y(6);
P=y(7);

Sdot = lambda-beta*S*B/(h_b+B)-beta_z*S*Z_B/(h_z+Z_B)-mu*S+w*R;

Idot = beta*S*B/(h_b+B)+beta_z*S*Z_B/(h_z+Z_B)-(gamma+mu+delta)*I;

Rdot = gamma*I-(mu+w)*R;

Bdot = xi*I-d_b*B-c*sigma*B*Z_F/(h_m+B);
% free-living bacteria

Z_Bdot=sigma*B*Z_F/(h_m+B)-d_z*Z_B;
% zooplankton associated with bacteria

Z_Fdot=eta*alpha*P*(Z_B+Z_F)/(h_p+P)-d_z*Z_F-sigma*B*Z_F/(h_m+B);
% uncolonized zooplankton

Pdot=r_p*P*(1-P/k_p)-alpha*P*(Z_B+Z_F)/(h_p+P);
% phytoplankton
 

ydot=[Sdot;Idot;Rdot;Bdot;Z_Bdot;Z_Fdot;Pdot];
end