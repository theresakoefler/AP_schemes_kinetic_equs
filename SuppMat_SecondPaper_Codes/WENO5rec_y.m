function ft=WENO5rec_y(f,fu_var,dy,alpha)
imax = size(f,1);
u(4:imax+3,:)=f;
u(3,:)=f(imax-1,:);
u(2,:)=f(imax-2,:);
u(1,:)=f(imax-3,:);
u(imax+4,:)=f(2,:);
u(imax+5,:)=f(3,:);
u(imax+6,:)=f(4,:);

fu(4:imax+3,:)=fu_var;
fu(3,:)=fu_var(imax-1,:);
fu(2,:)=fu_var(imax-2,:);
fu(1,:)=fu_var(imax-3,:);
fu(imax+4,:)=fu_var(2,:);
fu(imax+5,:)=fu_var(3,:);
fu(imax+6,:)=fu_var(4,:);

C00=1/3;
C01=5/6;
C02=-1/6;
C10=-1/6;
C11=5/6;
C12=1/3;
C20=1/3;
C21=-7/6;
C22=11/6;

Ct00=11/6;
Ct01=-7/6;
Ct02=1/3;
Ct10=C00;
Ct11=C01;
Ct12=C02;
Ct20=C10;
Ct21=C11;
Ct22=C12;


vip0(1:imax,:)=C00*u(4:imax+3,:)+C01*u(5:imax+4,:)+C02*u(6:imax+5,:);
vip1(1:imax,:)=C10*u(3:imax+2,:)+C11*u(4:imax+3,:)+C12*u(5:imax+4,:);
vip2(1:imax,:)=C20*u(2:imax+1,:)+C21*u(3:imax+2,:)+C22*u(4:imax+3,:);

vim0(1:imax,:)=Ct00*u(4:imax+3,:)+Ct01*u(5:imax+4,:)+Ct02*u(6:imax+5,:);
vim1(1:imax,:)=Ct10*u(3:imax+2,:)+Ct11*u(4:imax+3,:)+Ct12*u(5:imax+4,:);
vim2(1:imax,:)=Ct20*u(2:imax+1,:)+Ct21*u(3:imax+2,:)+Ct22*u(4:imax+3,:);

vip0_fu(1:imax,:)=C00*fu(4:imax+3,:)+C01*fu(5:imax+4,:)+C02*fu(6:imax+5,:);
vip1_fu(1:imax,:)=C10*fu(3:imax+2,:)+C11*fu(4:imax+3,:)+C12*fu(5:imax+4,:);
vip2_fu(1:imax,:)=C20*fu(2:imax+1,:)+C21*fu(3:imax+2,:)+C22*fu(4:imax+3,:);

vim0_fu(1:imax,:)=Ct00*fu(4:imax+3,:)+Ct01*fu(5:imax+4,:)+Ct02*fu(6:imax+5,:);
vim1_fu(1:imax,:)=Ct10*fu(3:imax+2,:)+Ct11*fu(4:imax+3,:)+Ct12*fu(5:imax+4,:);
vim2_fu(1:imax,:)=Ct20*fu(2:imax+1,:)+Ct21*fu(3:imax+2,:)+Ct22*fu(4:imax+3,:);

d0=3/10;
d1=3/5;
d2=1/10;

dt0=d2;
dt1=d1;
dt2=d0;

e=10^-6;

beta0=13/12*(u(4:imax+3,:)-2*u(5:imax+4,:)+u(6:imax+5,:)).^2+1/4*(3*u(4:imax+3,:)-4*u(5:imax+4,:)+u(6:imax+5,:)).^2;
beta1=13/12*(u(3:imax+2,:)-2*u(4:imax+3,:)+u(5:imax+4,:)).^2+1/4*(u(3:imax+2,:)-u(5:imax+4,:)).^2;
beta2=13/12*(u(2:imax+1,:)-2*u(3:imax+2,:)+u(4:imax+3,:)).^2+1/4*(u(2:imax+1,:)-4*u(3:imax+2,:)+3*u(4:imax+3,:)).^2;

beta0_fu=13/12*(fu(4:imax+3,:)-2*fu(5:imax+4,:)+fu(6:imax+5,:)).^2+1/4*(3*fu(4:imax+3,:)-4*fu(5:imax+4,:)+fu(6:imax+5,:)).^2;
beta1_fu=13/12*(fu(3:imax+2,:)-2*fu(4:imax+3,:)+fu(5:imax+4,:)).^2+1/4*(fu(3:imax+2,:)-fu(5:imax+4,:)).^2;
beta2_fu=13/12*(fu(2:imax+1,:)-2*fu(3:imax+2,:)+fu(4:imax+3,:)).^2+1/4*(fu(2:imax+1,:)-4*fu(3:imax+2,:)+3*fu(4:imax+3,:)).^2;

alpha0=d0./(e+beta0).^2;
alpha1=d1./(e+beta1).^2;
alpha2=d2./(e+beta2).^2;

alphat0=dt0./(e+beta0).^2;
alphat1=dt1./(e+beta1).^2;
alphat2=dt2./(e+beta2).^2;

alpha0_fu=d0./(e+beta0_fu).^2;
alpha1_fu=d1./(e+beta1_fu).^2;
alpha2_fu=d2./(e+beta2_fu).^2;

alphat0_fu=dt0./(e+beta0_fu).^2;
alphat1_fu=dt1./(e+beta1_fu).^2;
alphat2_fu=dt2./(e+beta2_fu).^2;

suma=alpha0+alpha1+alpha2;
omega0=alpha0./suma;
omega1=alpha1./suma;
omega2=alpha2./suma;

sumat=alphat0+alphat1+alphat2;
omegat0=alphat0./sumat;
omegat1=alphat1./sumat;
omegat2=alphat2./sumat;

suma_fu=alpha0_fu+alpha1_fu+alpha2_fu;
omega0_fu=alpha0_fu./suma_fu;
omega1_fu=alpha1_fu./suma_fu;
omega2_fu=alpha2_fu./suma_fu;

sumat_fu=alphat0_fu+alphat1_fu+alphat2_fu;
omegat0_fu=alphat0_fu./sumat_fu;
omegat1_fu=alphat1_fu./sumat_fu;
omegat2_fu=alphat2_fu./sumat_fu;


vipm=omega0.*vip0+omega1.*vip1+omega2.*vip2;
vimp=omegat0.*vim0+omegat1.*vim1+omegat2.*vim2;

vipm_fu=omega0_fu.*vip0_fu+omega1_fu.*vip1_fu+omega2_fu.*vip2_fu;
vimp_fu=omegat0_fu.*vim0_fu+omegat1_fu.*vim1_fu+omegat2_fu.*vim2_fu;

% Berechnung von u_i+1/2^(+)

vim0(1:imax,:)=Ct00*u(5:imax+4,:)+Ct01*u(6:imax+5,:)+Ct02*u(7:imax+6,:);
vim1(1:imax,:)=Ct10*u(4:imax+3,:)+Ct11*u(5:imax+4,:)+Ct12*u(6:imax+5,:);
vim2(1:imax,:)=Ct20*u(3:imax+2,:)+Ct21*u(4:imax+3,:)+Ct22*u(5:imax+4,:);

vim0_fu(1:imax,:)=Ct00*fu(5:imax+4,:)+Ct01*fu(6:imax+5,:)+Ct02*fu(7:imax+6,:);
vim1_fu(1:imax,:)=Ct10*fu(4:imax+3,:)+Ct11*fu(5:imax+4,:)+Ct12*fu(6:imax+5,:);
vim2_fu(1:imax,:)=Ct20*fu(3:imax+2,:)+Ct21*fu(4:imax+3,:)+Ct22*fu(5:imax+4,:);

dt0=d2;
dt1=d1;
dt2=d0;

beta0=13/12*(u(5:imax+4,:)-2*u(6:imax+5,:)+u(7:imax+6,:)).^2+1/4*(3*u(5:imax+4,:)-4*u(6:imax+5,:)+u(7:imax+6,:)).^2;
beta1=13/12*(u(4:imax+3,:)-2*u(5:imax+4,:)+u(6:imax+5,:)).^2+1/4*(u(4:imax+3,:)-u(6:imax+5,:)).^2;
beta2=13/12*(u(3:imax+2,:)-2*u(4:imax+3,:)+u(5:imax+4,:)).^2+1/4*(u(3:imax+2,:)-4*u(4:imax+3,:)+3*u(5:imax+4,:)).^2;

beta0_fu=13/12*(fu(5:imax+4,:)-2*fu(6:imax+5,:)+fu(7:imax+6,:)).^2+1/4*(3*fu(5:imax+4,:)-4*fu(6:imax+5,:)+fu(7:imax+6,:)).^2;
beta1_fu=13/12*(fu(4:imax+3,:)-2*fu(5:imax+4,:)+fu(6:imax+5,:)).^2+1/4*(fu(4:imax+3,:)-fu(6:imax+5,:)).^2;
beta2_fu=13/12*(fu(3:imax+2,:)-2*fu(4:imax+3,:)+fu(5:imax+4,:)).^2+1/4*(fu(3:imax+2,:)-4*fu(4:imax+3,:)+3*fu(5:imax+4,:)).^2;

alphat0=dt0./(e+beta0).^2;
alphat1=dt1./(e+beta1).^2;
alphat2=dt2./(e+beta2).^2;

sumat=alphat0+alphat1+alphat2;
omegat0=alphat0./sumat;
omegat1=alphat1./sumat;
omegat2=alphat2./sumat;

alphat0_fu=dt0./(e+beta0_fu).^2;
alphat1_fu=dt1./(e+beta1_fu).^2;
alphat2_fu=dt2./(e+beta2_fu).^2;

sumat_fu=alphat0_fu+alphat1_fu+alphat2_fu;
omegat0_fu=alphat0_fu./sumat_fu;
omegat1_fu=alphat1_fu./sumat_fu;
omegat2_fu=alphat2_fu./sumat_fu;

vipp=omegat0.*vim0+omegat1.*vim1+omegat2.*vim2;

vipp_fu=omegat0_fu.*vim0_fu+omegat1_fu.*vim1_fu+omegat2_fu.*vim2_fu;

% Berechnung u_i-1/2^(-)

vip0(1:imax,:)=C00*u(3:imax+2,:)+C01*u(4:imax+3,:)+C02*u(5:imax+4,:);
vip1(1:imax,:)=C10*u(2:imax+1,:)+C11*u(3:imax+2,:)+C12*u(4:imax+3,:);
vip2(1:imax,:)=C20*u(1:imax,:)+C21*u(2:imax+1,:)+C22*u(3:imax+2,:);

vip0_fu(1:imax,:)=C00*fu(3:imax+2,:)+C01*fu(4:imax+3,:)+C02*fu(5:imax+4,:);
vip1_fu(1:imax,:)=C10*fu(2:imax+1,:)+C11*fu(3:imax+2,:)+C12*fu(4:imax+3,:);
vip2_fu(1:imax,:)=C20*fu(1:imax,:)+C21*fu(2:imax+1,:)+C22*fu(3:imax+2,:);

d0=3/10;
d1=3/5;
d2=1/10;

e=10^-6;

beta0=13/12*(u(3:imax+2,:)-2*u(4:imax+3,:)+u(5:imax+4,:)).^2+1/4*(3*u(3:imax+2,:)-4*u(4:imax+3,:)+u(5:imax+4,:)).^2;
beta1=13/12*(u(2:imax+1,:)-2*u(3:imax+2,:)+u(4:imax+3,:)).^2+1/4*(u(2:imax+1,:)-u(4:imax+3,:)).^2;
beta2=13/12*(u(1:imax,:)-2*u(2:imax+1,:)+u(3:imax+2,:)).^2+1/4*(u(1:imax,:)-4*u(2:imax+1,:)+3*u(3:imax+2,:)).^2;

beta0_fu=13/12*(fu(3:imax+2,:)-2*fu(4:imax+3,:)+fu(5:imax+4,:)).^2+1/4*(3*fu(3:imax+2,:)-4*fu(4:imax+3,:)+fu(5:imax+4,:)).^2;
beta1_fu=13/12*(fu(2:imax+1,:)-2*fu(3:imax+2,:)+fu(4:imax+3,:)).^2+1/4*(fu(2:imax+1,:)-fu(4:imax+3,:)).^2;
beta2_fu=13/12*(fu(1:imax,:)-2*fu(2:imax+1,:)+fu(3:imax+2,:)).^2+1/4*(fu(1:imax,:)-4*fu(2:imax+1,:)+3*fu(3:imax+2,:)).^2;

alpha0=d0./(e+beta0).^2;
alpha1=d1./(e+beta1).^2;
alpha2=d2./(e+beta2).^2;

suma=alpha0+alpha1+alpha2;
omega0=alpha0./suma;
omega1=alpha1./suma;
omega2=alpha2./suma;

alpha0_fu=d0./(e+beta0_fu).^2;
alpha1_fu=d1./(e+beta1_fu).^2;
alpha2_fu=d2./(e+beta2_fu).^2;

suma_fu=alpha0_fu+alpha1_fu+alpha2_fu;
omega0_fu=alpha0_fu./suma_fu;
omega1_fu=alpha1_fu./suma_fu;
omega2_fu=alpha2_fu./suma_fu;

vimm=omega0.*vip0+omega1.*vip1+omega2.*vip2;
vimm_fu=omega0_fu.*vip0_fu+omega1_fu.*vip1_fu+omega2_fu.*vip2_fu;

% Reconstruction in x direction

ftp=1/2*(vipm+vipp-alpha*(vipp_fu-vipm_fu));
ftm=1/2*(vimp+vimm-alpha*(vimp_fu-vimm_fu));

ft=(ftp-ftm)/dy;