r_c=20;
cell_dim=13;
perc_lead=0.3;
% alpha=1;% damping factor
alpha = 0.7;

% alpha = 1; % gelma;
beta=0.1; % più è basso più arriva velocemente ad una velocità uniforme

a0=3*r_c; % repulsive pot
as=r_c; % repulsive from surface
a1=r_c*3; % attractive pot
al=r_c*4; %repulsive over surface from leader

r_lead=3*r_c; % è leader è vicina alle cellule di superfice
t0=4;
r_surface=2*r_c; % superfice vicine alla leader

% U0=10; % repulsivo lontano a0
U0=1; % repulsivo lontano a0
% U1=0.1; % attrattivo lontano a1
U1=10; % attrattivo lontano a1
Us=5; % repulsivo su normali
Ul=10; % repulsivo da leader su superfice

D_lead=1; % modulo drift leader
D0=2; % modulo drift cellule
D0_x = 20;%10
D0_y = 20;%10
vmax=0.5;

sigma_lead=1;

R_int_c=2*r_c; % generica definizione di vicini per i potenziali attrattivi e repulsivi

sigmax_0=3; % volatilità si può anche variare secondo la regola della densità rimanente rho
sigmax_1=5;
sigmay_0=3;
sigmay_1=5;
rho0=1;
n_cancer=50;% numero cellule

r_hbmscs = 30;
n_hbmscs = 8;



R_int_hbmscs = 7*r_hbmscs;
a0_hbmscs= r_hbmscs + r_c; % repulsive pot
% a0_hbmscs= 6*R_int_hbmscs; % repulsive pot

a1_hbmscs= 6*R_int_hbmscs; 
% a1_hbmscs= r_hbmscs; 





U0_hbmscs = 0.1;
% U1_hbmscs = 50;%50;%10
U1_hbmscs = 60;
% U1_hbmscs = 

R_int_c_far = 600;
U0_far = 9;%3
a0_far = R_int_c_far;


save('settings_file_B2B.mat')