function ft=WENO3rec_x(f,fu_var,dx,alpha)

imax = size(f,2);
u(:,3:imax+2)=f;
% periodic ghost nodes
u(:,1)=f(:,imax-2);
u(:,2)=f(:,imax-1);
u(:,imax+3)=f(:,2);
u(:,imax+4)=f(:,3);

fu(:,3:imax+2)=fu_var;
% periodic ghost nodes
fu(:,1)=fu_var(:,imax-2);
fu(:,2)=fu_var(:,imax-1);
fu(:,imax+3)=fu_var(:,2);
fu(:,imax+4)=fu_var(:,3);

% Wert für \epsilon
e=1e-6;

% Tabelle in Shu WENO notes
C00 = 1/2; C01 = 1/2;
C10 = -1/2; C11 = 3/2;

Ct00 = 3/2; Ct01 = -1/2;
Ct10 = C00; Ct11 = C01;


% Polynome für u_i+1/2^(r)
vip0(:,1:imax) = C00*u(:,3:imax+2)+C01*u(:,4:imax+3);
vip1(:,1:imax) = C10*u(:,2:imax+1)+C11*u(:,3:imax+2);

% Polynome für u_i-1/2^(r)
vim0(:,1:imax) = Ct00*u(:,3:imax+2)+Ct01*u(:,4:imax+3);
vim1(:,1:imax) = Ct10*u(:,2:imax+1)+Ct11*u(:,3:imax+2);

% Polynome für u_i+1/2^(r)
vip0_fu(:,1:imax) = C00*fu(:,3:imax+2)+C01*fu(:,4:imax+3);
vip1_fu(:,1:imax) = C10*fu(:,2:imax+1)+C11*fu(:,3:imax+2);

% Polynome für u_i-1/2^(r)
vim0_fu(:,1:imax) = Ct00*fu(:,3:imax+2)+Ct01*fu(:,4:imax+3);
vim1_fu(:,1:imax) = Ct10*fu(:,2:imax+1)+Ct11*fu(:,3:imax+2);


d0 = 1/2;%2/3; 
d1 = 1/2;%1/3;
dt0 = d1; dt1 = d0;


beta0 = (u(:,4:imax+3) - u(:,3:imax+2)).^2; 
beta1 = (u(:,3:imax+2) - u(:,2:imax+1)).^2;

beta0_fu = (fu(:,4:imax+3) - fu(:,3:imax+2)).^2; 
beta1_fu = (fu(:,3:imax+2) - fu(:,2:imax+1)).^2;

alpha0 = d0./(e+beta0).^2;
alpha1 = d1./(e+beta1).^2;

alphat0=dt0./(e+beta0).^2;
alphat1=dt1./(e+beta1).^2;

alpha0_fu = d0./(e+beta0_fu).^2;
alpha1_fu = d1./(e+beta1_fu).^2;

alphat0_fu=dt0./(e+beta0_fu).^2;
alphat1_fu=dt1./(e+beta1_fu).^2;

suma=alpha0+alpha1;
omega0=alpha0./suma;
omega1=alpha1./suma;

sumat=alphat0+alphat1;
omegat0=alphat0./sumat;
omegat1=alphat1./sumat;

suma_fu=alpha0_fu+alpha1_fu;
omega0_fu=alpha0_fu./suma_fu;
omega1_fu=alpha1_fu./suma_fu;

sumat_fu=alphat0_fu+alphat1_fu;
omegat0_fu=alphat0_fu./sumat_fu;
omegat1_fu=alphat1_fu./sumat_fu;

vipm=omega0.*vip0+omega1.*vip1; % u_i+1/2^(-)
vimp=omegat0.*vim0+omegat1.*vim1; % u_i-1/2^(+)

vipm_fu=omega0_fu.*vip0_fu+omega1_fu.*vip1_fu; % u_i+1/2^(-)
vimp_fu=omegat0_fu.*vim0_fu+omegat1_fu.*vim1_fu; % u_i-1/2^(+)

% Berechnung von  u_i+1/2^(+)
vim0(:,1:imax) = Ct00*u(:,4:imax+3)+Ct01*u(:,5:imax+4);
vim1(:,1:imax) = Ct10*u(:,3:imax+2)+Ct11*u(:,4:imax+3);

vim0_fu(:,1:imax) = Ct00*fu(:,4:imax+3)+Ct01*fu(:,5:imax+4);
vim1_fu(:,1:imax) = Ct10*fu(:,3:imax+2)+Ct11*fu(:,4:imax+3);


beta0 = (u(:,5:imax+4) - u(:,4:imax+3)).^2; 
beta1 = (u(:,4:imax+3) - u(:,3:imax+2)).^2;

beta0_fu = (fu(:,5:imax+4) - fu(:,4:imax+3)).^2; 
beta1_fu = (fu(:,4:imax+3) - fu(:,3:imax+2)).^2;

alphat0=dt0./(e+beta0).^2;
alphat1=dt1./(e+beta1).^2;

sumat=alphat0+alphat1;
omegat0=alphat0./sumat;
omegat1=alphat1./sumat;

alphat0_fu=dt0./(e+beta0_fu).^2;
alphat1_fu=dt1./(e+beta1_fu).^2;

sumat_fu=alphat0_fu+alphat1_fu;
omegat0_fu=alphat0_fu./sumat_fu;
omegat1_fu=alphat1_fu./sumat_fu;


vipp=omegat0.*vim0+omegat1.*vim1; % u_i+1/2^(+)
vipp_fu=omegat0_fu.*vim0_fu+omegat1_fu.*vim1_fu; % u_i+1/2^(+)

% Berechnung von  u_i-1/2^(-)
vip0(:,1:imax) = C00*u(:,2:imax+1)+C01*u(:,3:imax+2);
vip1(:,1:imax) = C10*u(:,1:imax)+C11*u(:,2:imax+1);

vip0_fu(:,1:imax) = C00*fu(:,2:imax+1)+C01*fu(:,3:imax+2);
vip1_fu(:,1:imax) = C10*fu(:,1:imax)+C11*fu(:,2:imax+1);


beta0 = (u(:,3:imax+2) - u(:,2:imax+1)).^2; 
beta1 = (u(:,2:imax+1) - u(:,1:imax)).^2;

beta0_fu = (fu(:,3:imax+2) - fu(:,2:imax+1)).^2; 
beta1_fu = (fu(:,2:imax+1) - fu(:,1:imax)).^2;

alpha0=d0./(e+beta0).^2;
alpha1=d1./(e+beta1).^2;

suma=alpha0+alpha1;
omega0=alpha0./suma;
omega1=alpha1./suma;

alpha0_fu=d0./(e+beta0_fu).^2;
alpha1_fu=d1./(e+beta1_fu).^2;

suma_fu=alpha0_fu+alpha1_fu;
omega0_fu=alpha0_fu./suma_fu;
omega1_fu=alpha1_fu./suma_fu;

vimm=omega0.*vip0+omega1.*vip1; % u_i-1/2^(-)
vimm_fu=omega0_fu.*vip0_fu+omega1_fu.*vip1_fu; % u_i-1/2^(-)

% Flux

ftp=1/2*(vipm+vipp-alpha*(vipp_fu-vipm_fu));%...

ftm=1/2*(vimp+vimm-alpha*(vimp_fu-vimm_fu));...


ft=(ftp-ftm)/dx;

