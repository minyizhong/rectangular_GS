%% generate the grids and GS_wigner
close all 
clear all

path='D:\ILLUMINATION PROJECT\Figures & Reports\step phase mask\delta_z=2um, lambda=2um\';
% path='C:\Users\mzhong\Desktop\temp';

wo=2;%0.08;        % beam waist 
z=2;         % lateral distance from the beam waist
lambda=2e-3;%0.7e-3; % wavelength in millimeters
lc=1;%0.01;        % lateral coherent length in millimeters

% nx=512*2*2;  max_x=4;%4;
nx=512;     %length(Fi);  % sum of sampling points in x direction
max_x=4;    % maximum value on x axis
xi=gnr_array(nx,max_x);
dx=abs(xi(1)-xi(2));

x_prime=xi*2;
dx_prime=x_prime(2)-x_prime(1);

nu=nx;      % sum of sampling points in u direction
max_u=1/dx_prime*lambda/2; %0.0025;%  unit: radian, 1/dx_prime: half of a fourier transform, scale: lambda/2: half of a period
ui=gnr_array(nu,max_u);

[x,xx] = meshgrid(xi, x_prime);%(xi,xxi);
% % x changes in columns, xx changes in lines
[xx2,u]=meshgrid(x_prime,ui);%(xxi,ui);
% %xx2 is actually xx, but changes in columns. u is the frequencies changing in
% %lines.

x1=x+xx/2;
x2=x-xx/2;

%% GS beam

[GS_Corr,Wig_GS]=func_gnr_GS_wig(nx,nu,max_x,max_u,wo,z,lambda,lc,path);

%% Rectangular Gaussian Schell Model beam

sigma=6e-4; % x1-x2 direction; inverse propotional
% proportional to angle expansion in wigner

delta=5e-4; % x1+x2 direction; inverse propotional
% inverse proportional to x expansion in wigner

theta=0.1;

M=50; 

r=1e6;

[Corr_rect_GS,Wig_rect]=func_gnr_rectGS_wig(nx,max_x,lambda,M,r,sigma,delta,theta,path);

%% Step Phase Mask (SPM)

SPM=zeros(1,nx);

fac=1;

SPM(xi<0)=lambda*fac;

% SPM(xi==0)=lambda*fac*0.5;
% 
% SPM(find(xi==0)+1)=lambda*fac*0.25;
% 
% SPM(find(xi==0)-1)=lambda*fac*0.75;
 
Fi_true=SPM;

plot_1d(xi,'x (mm)',Fi_true,'z (mm)','r-o','Step Phase Mask',52,path)

%% wigner calculation

F1=interp1(xi,Fi_true,x1);
F2=interp1(xi,Fi_true,x2);

F1(x1>max(xi)| x1<min(xi))=0;
F2(x2>max(xi)| x2<min(xi))=0;

F1(isnan(F1))=0;
F2(isnan(F2))=0;

n=2;

F_phase1=exp(-1i*2*pi*(n-1)/lambda*F1);
F_phase2=exp(-1i*2*pi*(n-1)/lambda*F2);

% Corr=Corr_rect_GS.*F_phase1.*conj(F_phase2);

Corr=GS_Corr.*F_phase1.*conj(F_phase2);

Phase=exp(-2i*pi/lambda*u.*xx2);

Wig=real(Phase*Corr*dx);

% plot_2d(xi,'x (mm)',ui,'u (radian)',Wig,'Wigner after phase mask',16); %
% without normalization

I_x_step=sum(Wig)/max(sum(Wig));

I_u_step=sum(Wig,2)/max(sum(Wig,2));

Wig_step=Wig/max(max(Wig));

plot_1d(xi,'x(mm)',I_x_step,'I','-','I(x) after phase mask',15,path);

plot_1d(ui,'u(radian)',I_u_step,'I','-','I(u) after phase mask',18,path);

plot_2d(xi,'x (mm)',ui,'u (radian)',Wig_step,'Wigner after phase mask',14,path); % with normalization



%%
z=100;

step_z=100;

Wig=Wig/max(max(Wig));

wigner_z_propagation(nx,xi,ui,z,step_z,Wig,path);
