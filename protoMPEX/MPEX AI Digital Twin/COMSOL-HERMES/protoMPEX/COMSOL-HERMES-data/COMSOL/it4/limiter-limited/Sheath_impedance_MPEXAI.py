from math import sqrt, pi, atan2, tanh, log

def get_phi0avg(LayerData, EdgeData, init=0, unmagnetized=0):
 phi0avgList=[]
 for lRow,eRow in zip(LayerData,EdgeData):
  d=0.002; #[m] assuming data is at middle of layer
  x=lRow[0].real
  y=lRow[1].real
  z=lRow[2].real
  phi = atan2(y,x);
  E_norm = (lRow[3]*x+lRow[4]*y)/sqrt(x**2+y**2)
  Vlayer = abs(E_norm)*d;
  Density = eRow[3].real;
  Br = eRow[4];
  B0 = eRow[5];
  bn = abs(Br/B0);
  # Variable Creation
  mu = 24.17; # for Deuterium; upgrade to arbitrary single ion species is in progress
  e=1.602E-19; #J to eV, also elementary charge
  epsilon_0org=8.8419E-12;  #F/m or C^2/J or s^2*C^2/m^3*kg %Permittivity of free space
  epsilon_0=epsilon_0org*e; #C^2/eV %Permittivity of free space
  m_D_kg=3.34E-27; #kg %mass of D ion
  w=13.56e6*2*pi; #rad/s
  #omega=0; %0 for unmagnetized case; omega_ci for magnetized case
  #bx=1*sind(90);
  #xi=14;
  #bn=1; %1 for unmagnetized limit
  LayerThickness=d; #m
  Te=8.0; #eV
  Ti=0; #0; %eV
  Z=1; #D
  A=2; #D
  thick_debye=sqrt((epsilon_0*Te)/(Density*e**2)); #m %Chen2015
  omegapi=1.32E3*Z*sqrt(Density*1e-6/A); #cgs units %rad/s
  omega_hat=w/omegapi;
  w=omega_hat;
  omegaci=7.63e6*B0*6.28; #rad/s for D, Dolan2013 eq.5.33
  omega_ci_hat=omegaci/omegapi;
  omega=omega_ci_hat;
  bx=bn;
  xi=1.0*Vlayer/Te;
  if unmagnetized == 1:
   omega = 0; 
   bx = 1; 
  j=0;
  upar0=1.1;
  # phi0avg -> Myra says this is DC voltage
  a1 = 3.70285;
  a2 = 3.81991;
  b1 = 1.13352;
  b2 = 1.24171;
  a3 = 2.0*b2/pi;
  c0=0.966463;
  c1=0.141639;
  gg=c0+c1*tanh(w);
  xi1=gg*xi;
  ff=((log(mu)+xi1*a1+xi1**2.*a2+xi1**3*a3)/(1+xi1*b1+xi1**2*b2))-log(1-(j/upar0))+log(mu/24.17);
  phi0avg=ff;
  phi0avgList.append(phi0avg)
 return phi0avgList;


def Sheath_impedance_MPEXAI(filename, filename2, init=0, unmagnetized=0):
 computedEpsilon=[]
 computedSigma=[]
 LayerData=[]
 with open(filename) as f:
  for line in f.readlines():
   line=line.strip()
   if not line.startswith('%'):
    line=line.split(',')
    LayerData.append([])
    for v in line:
     if 'i' in v or '+' in v:
      #print(v)
      LayerData[-1].append(complex(v.replace('i','j')))
     else:
      LayerData[-1].append(float(v))
 EdgeData=[]
 with open(filename2) as f:
  for line in f.readlines():
   line=line.strip()
   if not line.startswith('%'):
    line=line.split(',')
    EdgeData.append([])
    for v in line:
     EdgeData[-1].append(float(v))
 for lRow,eRow in zip(LayerData,EdgeData):
  d=0.002; #[m] assuming data is at middle of layer
  x=lRow[0]
  y=lRow[1]
  z=lRow[2]
  phi = atan2(y,x);
  E_norm = (lRow[3]*x+lRow[4]*y)/sqrt(x**2+y**2)
  Vlayer = abs(E_norm)*d;
  Density = eRow[3];
  Br = eRow[4];
  B0 = eRow[5];
  bn = abs(Br/B0);
  # Variable Creation
  mu = 24.17; # for Deuterium; upgrade to arbitrary single ion species is in progress
  e=1.602E-19; #J to eV, also elementary charge
  epsilon_0org=8.8419E-12;  #F/m or C^2/J or s^2*C^2/m^3*kg %Permittivity of free space
  epsilon_0=epsilon_0org*e; #C^2/eV %Permittivity of free space
  m_D_kg=3.34E-27; #kg %mass of D ion
  w=13.56e6*2*pi; #rad/s
  #omega=0; %0 for unmagnetized case; omega_ci for magnetized case
  #bx=1*sind(90);
  #xi=14;
  #bn=1; %1 for unmagnetized limit
  LayerThickness=d; #m
  Te=8.0; #eV
  Ti=0; #0; %eV
  Z=1; #D
  A=2; #D
  thick_debye=sqrt((epsilon_0*Te)/(Density*e**2)); #m %Chen2015
  omegapi=1.32E3*Z*sqrt(Density*1e-6/A); #cgs units %rad/s
  omega_hat=w/omegapi;
  w=omega_hat;
  omegaci=7.63e6*B0*6.28; #rad/s for D, Dolan2013 eq.5.33
  omega_ci_hat=omegaci/omegapi;
  omega=omega_ci_hat;
  bx=bn;
  xi=1.0*Vlayer/Te;
  if unmagnetized == 1:
   omega = 0; 
   bx = 1; 
  j=0;
  upar0=1.1;
  # phi0avg -> Myra says this is DC voltage
  a1 = 3.70285;
  a2 = 3.81991;
  b1 = 1.13352;
  b2 = 1.24171;
  a3 = 2.0*b2/pi;
  c0=0.966463;
  c1=0.141639;
  gg=c0+c1*tanh(w);
  xi1=gg*xi;
  ff=((log(mu)+xi1*a1+xi1**2.*a2+xi1**3*a3)/(1+xi1*b1+xi1**2*b2))-log(1-(j/upar0))+log(mu/24.17);
  phi0avg=ff;
  # ni1avg
  # niw = ion density at the wall for a static sheath
  # vv is the total dc potential drop 
  #niw(omega,bx, xi)
  #niww(w,omega,bx,xi)
  # w=0.2;
  # omega=0.3;
  # bx=0.4;
  # xi=13;
  k0=3.7616962640756197;
  k1=0.2220204461728174;
  phiavg=phi0avg;
  philowomega=k0+k1*(xi-k0)-log(1-(j/upar0));
  phimod=philowomega+(phiavg-philowomega)*tanh(w);
  d0=0.7944430930529499;
  d1=0.803531266389172;
  d2=0.18237897510951012;
  d3=0.9957212047604492;
  nu1=1.4555923231100891;
  arg=sqrt(((mu**2*bx**2+1)/(mu**2+1)));
  fff=-log(arg)/(1+d3*omega**2);
  Phip=phimod-fff;
  Phip1=0.0;
  Phip1=Phip;
  if Phip<0:
   Phip1=0.0;
  omegaPhi=omega*Phip1**(1.0/4.0);
  d4=(d2**2)/((mu**2*d0**2)-d2**2);
  niw=(d0/(d2+sqrt(Phip1)))*sqrt((bx**2+d4+d1**2*omegaPhi**(2.0*nu1))/(1+d4+d1**2*omegaPhi**(2.0*nu1)));
  # niwOmega is the ion density at the wall for an rf sheath
  # Xi is the 0-peak rf voltage
  niwomega=niw.real;
  # yd
  s0=1.1241547327789232;
  phi0a=phi0avg;
  niwomegaa=niwomega;
  Delta=sqrt(phi0a/niwomegaa);
  yd=-1j*s0*w/Delta;
  # ye
  # bx=0.4;
  # xi=3.6;
  h1=0.607405123251634;
  h2=0.3254965671158986;
  g1=0.6243920388599393;
  g2=0.5005946718280853;
  g3=(pi/4.0)*h2;
  he=(1+xi*h1+xi**2*h2)/(1+xi*g1+xi**2*g2+xi**3*g3);
  h0=1.05704235;
  ye=h0*abs(bx)*he*(1-(j/upar0));
  # yi
  parp0=1.0555369617763768;
  parp1=0.7976591020008023;
  parp2=1.47404874815277;
  parp3=0.8096145628336325;
  wcup=parp3*w/sqrt(niwomegaa);
  ycup=abs(bx)/(niwomegaa*sqrt(phi0a));
  epsilon=0.0001;
  gsmall=(w**2-bx**2*omega**2+1j*epsilon)/(w**2-omega**2+1j*epsilon);
  yi0=niwomegaa/sqrt(phi0a);
  yi=parp0*yi0*((1j*wcup))/((wcup**2/gsmall)-parp1+1j*parp2*ycup*wcup);
  # ytot and ztot
  ytot=yi+ye+yd;
  ztot=1./ytot;
  # Epsilon and Sigma Calculations
  computedEpsilon.append([phi,z, -((ytot.imag)/omega_hat)*(LayerThickness/thick_debye)]); 
  computedSigma.append([phi,z,epsilon_0org*omegapi*(LayerThickness/thick_debye)*(ytot.real)]);
 with open('epsilon_py.csv', 'w') as f:
  for line in computedEpsilon:
   f.write(str(line[0])+','+str(line[1])+','+str(line[2])+'\n');
 with open('sigma_py.csv', 'w') as f:
  for line in computedSigma:
   f.write(str(line[0])+','+str(line[1])+','+str(line[2])+'\n');
