

  %BOUBACAR KANDE Mémoire Modélisation DSGE  M2 APP DIJON


//clear all;
var   y_e mc_e ppi_e r_e a_e a ygap ppi_h rnat r
y ybar ppi s q e p_h  v v_e  c nx mc_he b n be ne;


varexo eps_a_e eps_a eps_v eps_v_e eps_n eps_b eps_ne eps_be eps_ppi;

parameters  alpha rho sigma upsilon varphi  beta lambda mu varphi_e varphi_ea omega_a
kappa_a rho_a theta rho_pi rho_ygap rho_e Ggamma_a eta tau gTheta_a

rho_ae rho_v rho_ve rho_n rho_ne rho_b rho_be ;

mu=1.2;
beta = 0.99;
sigma = 1;
alpha = 0.5;
eta = 1;
varphi = 3;
theta = 0.75;
rho=-log(beta);
tau=mu+log(1-alpha);
upsilon=-log(1-tau);
lambda = (1-(beta*theta))*(1-theta)/theta;
omega_a = 1 + alpha*(2-alpha)*(sigma*eta - 1);
kappa_a = lambda*(varphi + (sigma/omega_a));
gTheta_a = sigma*(1-omega_a)/(sigma + varphi*omega_a);

varphi_ea = -(sigma*(1+varphi)*(1-rho_ae)) / (varphi+sigma);
 varphi_e = 1.01;
rho_a = 0.75;
rho_ae = 0.75;
rho_v = 0.5;
rho_v_e = 0.5;
rho_r=0.75;
rho_pi= 1.5;
rho_ygap=0.1;
rho_e=0.2;
rho_n=0.65;
rho_ne=0.65;
rho_b=0.25;
rho_be=0.25;

model(linear);
y_e = be(+1)-be + y_e(+1)-(r_e-ppi_e(+1)-rho)/sigma;
mc_e = -upsilon+ne+(sigma+varphi)*y_e-(1+varphi)*a_e;

ppi_e = beta*ppi_e(+1)+lambda*mc_he;
mc_he=mc_e+mu;

r_e=varphi_e*ppi_e(+1)+varphi_ea*a_e+v_e;
ygap=be(+1)-be+ygap(+1)-(omega_a/sigma)*(r-ppi_h-rnat);

ppi_h=beta*ppi_h(+1)+kappa_a*ygap;

rnat=rho-(sigma*(1+varphi)*(1-rho_a))/(sigma+varphi*omega_a)*a-gTheta_a*(y_e(+1)-y_e);
r=rho_pi*ppi+rho_ygap*ygap+rho_e*(e-e(-1))+v;

ybar=Ggamma_a+(sigma*(1-omega_a)/(sigma+varphi*omega_a))*y_e;
y=ybar+ygap;

y=y_e+(omega_a/sigma)*s;
q=(1-alpha)*s;

ppi=ppi_h+alpha*(s-s(-1))+eps_ppi;
s=e+ppi_e-ppi_h;


p_h=p_h(-1)+ppi_h;

c=(1-alpha)*y+alpha*y_e;
nx=alpha*((2-alpha)*(sigma*eta-1)+(1-sigma))/omega_a*(y-y_e);

a=rho_a*a(-1)+eps_a;
a_e=rho_ae*a_e(-1)+eps_a_e;
v=rho_v*v(-1)+eps_v;
v_e=rho_ve*v_e(-1)+eps_v_e;
n=rho_n*n(-1)+eps_n;
ne=rho_ne*ne(-1)+eps_ne;
b=rho_b*b(-1)+eps_b;
be=rho_be*be(-1)+eps_be;

end;

  

steady;
check;

shocks;
var eps_v; stderr 0.075;
var eps_v_e; stderr 0.05;
var eps_a; stderr 0.1;
var eps_a_e; stderr 0.1;
var eps_n; stderr 0.1;
var eps_ne; stderr 0.1;
var eps_b ;stderr 0.1;
var eps_be ;stderr 0.1;
var eps_ppi ;stderr 0.1;
end;

stoch_simul(hp_filter = 1600, order = 1, irf = 90);


