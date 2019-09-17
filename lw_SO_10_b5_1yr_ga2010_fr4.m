%=====================================================================================
% Motivation | To calculate the generation rate of lee waves in the Southern Ocean.
%=====================================================================================
% Author     | Luwei Yang
%=====================================================================================
% Log        | 2017-05-19 Mask shallow areas; 
%            |            Add eddy contribution to mean (Cdisp,eperct);
%            |            Consider finite topography effects.
%            |            GA2010.
%            |            Steepness parameter = 0.4 & 4.0
%=====================================================================================

clear

disp('Running')

tic

dir2 = '/short/v45/lxy581/pyfiles/decomp/decomp_KE_disp_SO_10_160804.nc';
hu  = ncread([dir2],'hu');

dir3 = '/short/v45/lxy581/pyfiles/decomp/decomp_KE_bot_500_SO_10_160824_1yr.nc';

ub  = ncread([dir3],'u5');
vb  = ncread([dir3],'v5');
um  = ncread([dir3],'u5m');
vm  = ncread([dir3],'v5m');
up  = ncread([dir3],'u5p');
vp  = ncread([dir3],'v5p');

lon = ncread([dir3],'lon');lon=double(lon);
lat = ncread([dir3],'lat');lat=double(lat);
f   = sw_f(lat); f=double(f); 

%N  = load('/short/v45/lxy581/data/n2_mat/n2_SO_10_171117_b5_1yr.mat','n2s');
%ns = real(sqrt(N.n2s));
N  = load('/short/v45/lxy581/data/n2_mat/n2_SO_10_171204_b5_1yr.mat','nsm_b5');
ns = real(N.nsm_b5);

nx = size(lon,1);
ny = size(lat,1);
nt = size(ub,3);

iwx = zeros(nx,ny);
iwy = zeros(nx,ny);
PW  = zeros(nx,ny);
Edisp = zeros(nx,ny);
Cdisp = zeros(nx,ny);

hu(isnan(hu))=0;
% mask - added on 10 Nov 2016
ub(hu<1000)=nan;
vb(hu<1000)=nan;
um(hu<1000)=nan;
vm(hu<1000)=nan;
up(hu<1000)=nan;
vp(hu<1000)=nan;
ns(hu<1000)=nan;


% Topography - GA2010 (less coverage one)
topo = load('goff.mat');
hrms = topo.hrms;
ks   = topo.ks; 
kn   = topo.kn;
nu   = topo.nu;
theta = topo.theta;
tlon = topo.lon(1,:); 
tlat = topo.lat(:,1); 

%!!!!!!!
ks = ks*2*pi;
kn = kn*2*pi;

% Topo for SO (180W~180E,40S~65S)
ii = find(tlon>=-180 & tlon<=180);
jj = find(tlat>=-65 & tlat<=-40);

shrms = hrms(jj,ii);
sks = ks(jj,ii);
skn = kn(jj,ii);
stheta = theta(jj,ii);
snu = nu(jj,ii);
stlon = tlon(ii);
stlat = tlat(jj);

clear topo hrms ks kn nu theta

disp('Interpolation')
% Interpolation
[lonn,latt] = meshgrid(stlon,stlat);
slon = -180:0.1:180-0.1; slat = lat;
[SLON,SLAT] = meshgrid(slon,slat);
ihrms  = interp2(lonn,latt,shrms,SLON,SLAT);
iks    = interp2(lonn,latt,sks,SLON,SLAT);
ikn    = interp2(lonn,latt,skn,SLON,SLAT);
itheta = interp2(lonn,latt,stheta,SLON,SLAT);
inu    = 10.^interp2(lonn,latt,log10(snu),SLON,SLAT);  % log

clear shrms sks skn stheta snu lonn latt tlon tlat

% Change the topo data from 180W~180E to 80E~180E~0~80E
shrms(:,1:1000) = ihrms(:,2601:3600);
shrms(:,1001:3600) = ihrms(:,1:2600);
sks(:,1:1000) = iks(:,2601:3600);
sks(:,1001:3600) = iks(:,1:2600);
skn(:,1:1000) = ikn(:,2601:3600);
skn(:,1001:3600) = ikn(:,1:2600);
stheta(:,1:1000) = itheta(:,2601:3600);
stheta(:,1001:3600) = itheta(:,1:2600);
snu(:,1:1000) = inu(:,2601:3600);
snu(:,1001:3600) = inu(:,1:2600);

clear ihrms iks ikn itheta inu

% Wavenumber space
kmax = 2*pi/576;
k = linspace(-kmax,kmax,25);
l = k;

dk = k(2) - k(1);
dl = l(2) - l(1);
intk = dk*dl;

[k,l] = meshgrid(k,l);

%Cartesian -> Polar coordinate
kappa = sqrt(k.*k + l.*l);
phi = pi/2 - atan2(l,k);

%Topography spectrum
P0 = 4*pi.*snu.*shrms.^2./sks./skn;
num = 0;
cdi = 4*pi*pi; % common divide

%Correction for finite topography effects
iFrc = 0.4; % critical inverse Froude number
%iFrc = 4.0; % critical inverse Froude number

% 2016-11-10 add the Cdisp - only use the mean velocity in lee wave energy flux calculation
% Cdisp = tao(ubar) x ubar
for i = 1:nx;
  for j = 1:ny;
    sigma = um(i,j)*k + vm(i,j)*l;
    sf = (sigma.*sigma - f(j).*f(j));
    mu = get_mu(k,l,sigma,ns(i,j),f(j),0,2*hu(i,j),false,true);
    sq1 = (kappa/sks(j,i).*cos(phi-stheta(j,i)));
    sq2 = (kappa/skn(j,i).*sin(phi-stheta(j,i)));
    P = P0(j,i) ./ ((1 + sq1.*sq1 + sq2.*sq2).^(snu(j,i)+1));
    P(mu==0 | isnan(mu) | abs(mu)<(pi/hu(i,j))) = NaN;
    hrms(j,i) = nansum(P(:)*intk);
    hrms(j,i) = sqrt(hrms(j,i))/2/pi;
    P = P./(k.*k + l.*l).*mu.*sf;
    iFr = hrms(j,i)*ns(i,j)/sqrt(um(i,j).^2+vm(i,j).^2);
    if iFr <= iFrc;
      P = P;
    else
      P = P*(iFrc/iFr).^2;
    end  
    pw = P.*sigma;
    Cdisp_nl(i,j) = 1e+3*nansum(pw(:)*intk)/cdi;  %1e+3 stands for density 
  end
end
toc

disp('Start computing drag and disp for every time steps')

tic
for t = 1:nt;
  disp('==========')
  disp(t)
  num = num + 1;
  for i = 1:nx;
    for j = 1:ny;
      sigma = ub(i,j,t)*k + vb(i,j,t)*l;
      sf = sigma.*sigma - f(j).*f(j);
      mu = get_mu(k,l,sigma,ns(i,j),f(j),0,2*hu(i,j),false,true);
      sq1 = (kappa/sks(j,i).*cos(phi-stheta(j,i)));
      sq2 = (kappa/skn(j,i).*sin(phi-stheta(j,i)));
      P = P0(j,i) ./ (1 + sq1.*sq1 + sq2.*sq2).^(snu(j,i)+1);
      P(mu==0 | isnan(mu) | abs(mu)<(pi/hu(i,j))) = NaN;
      hrms(j,i) = nansum(P(:)*intk);
      hrms(j,i) = sqrt(hrms(j,i))/2/pi;
      P = P./(k.*k + l.*l).*mu.*sf;
      iFr = hrms(j,i)*ns(i,j)/sqrt(ub(i,j,t).^2+vb(i,j,t).^2);
      if iFr <= iFrc;
        P = P;
      else
        P = P*(iFrc/iFr).^2;
      end  
      pw = P.*sigma;
      pw_int = 1e+3*nansum(pw(:)*intk)/cdi;
      PW(i,j) = PW(i,j) + pw_int;
      tlwx = P.*k; 
      tlwy = P.*l;
      iwx(i,j) = iwx(i,j) + nansum(tlwx(:)*intk);
      iwy(i,j) = iwy(i,j) + nansum(tlwy(:)*intk);  
    end
  end
end
toc
disp('Loop1 Finished')

miwx = 1e+3.*iwx/num./cdi;
miwy = 1e+3.*iwy/num./cdi;

%total - mdisp
mdisp_nl = PW/num;

%mean - Mdisp
Mdisp_nl = miwx.*(um) + miwy.*(vm);

tic
for t = 1:nt;
  disp('-----')
  disp(t)
  for i = 1:nx;
    for j = 1:ny;
      sigma = ub(i,j,t)*k + vb(i,j,t)*l;
      sf = sigma.*sigma - f(j).*f(j);
      mu = get_mu(k,l,sigma,ns(i,j),f(j),0,2*hu(i,j),false,true);
      sq1 = (kappa/sks(j,i).*cos(phi-stheta(j,i)));
      sq2 = (kappa/skn(j,i).*sin(phi-stheta(j,i)));
      P = P0(j,i) ./ ((1 + sq1.*sq1 + sq2.*sq2).^(snu(j,i)+1));
      P(mu==0 | isnan(mu) | abs(mu)<(pi/hu(i,j))) = NaN;
      hrms(j,i) = nansum(P(:)*intk);
      hrms(j,i) = sqrt(hrms(j,i))/2/pi;
      P = P./(k.*k + l.*l).*mu.*sf;      
      iFr = hrms(j,i)*ns(i,j)/sqrt(ub(i,j,t).^2+vb(i,j,t).^2);
      if iFr <= iFrc;
        P = P;
      else
        P = P*(iFrc/iFr).^2;
      end      
      tlwx = P.*k; 
      tlwy = P.*l;
      eiwx = 1e+3*nansum(tlwx(:)*intk)/cdi - miwx(i,j);
      eiwy = 1e+3*nansum(tlwy(:)*intk)/cdi - miwy(i,j);
      edisp = eiwx.*up(i,j,t) + eiwy.*vp(i,j,t);
      Edisp(i,j) = Edisp(i,j) + edisp;         
    end
  end
end
toc

mEdisp_nl = Edisp/num;

% More accurate results (mask+Cdisp+finite_topo)

% iFrc = 0.4
% save /short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_ga2010_170519_fr4.mat lon lat mdisp_nl Mdisp_nl mEdisp_nl Cdisp_nl 
% save /short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_ga2010_170527_fr4.mat lon lat mdisp_nl Mdisp_nl mEdisp_nl Cdisp_nl 
% save /short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_ga2010_170531_fr4.mat lon lat mdisp_nl Mdisp_nl mEdisp_nl Cdisp_nl 
% Use bottom 500m-averaged last year-averaged stratification
%save /short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_ga2010_171118_fr4.mat lon lat mdisp_nl Mdisp_nl mEdisp_nl Cdisp_nl
% Correct the unit for stratification (s-1 to rad s-1)
%save /short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_ga2010_180123_fr4.mat lon lat mdisp_nl Mdisp_nl mEdisp_nl Cdisp_nl
% remove *2pi for N
save /short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_ga2010_180129_fr4.mat lon lat mdisp_nl Mdisp_nl mEdisp_nl Cdisp_nl miwx miwy

% iFrc = 4.0
% save /short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_ga2010_170528_fr4.mat lon lat mdisp_nl Mdisp_nl mEdisp_nl Cdisp_nl 

 
disp('Finished!')
