close all
clear all

lambda=2e-3;

nx=512;     %length(Fi);  % sum of sampling points in x direction
max_x=4;    % maximum value on x axis
nu=nx;      % sum of sampling points in u direction
dx = max_x*2 / (nx-2);
max_u=1/dx/2*lambda/2; %0.0318  unit: radian, 1/dx/2: half of a fourier transform, scale: lambda/2: half of a period

[xi,ui]=gnr_array(nx,nx,max_x,max_u);

x_prime=xi*2;

[x,xx] = meshgrid(xi, x_prime);%(xi,xxi);
% % x changes in columns, xx changes in lines
[xx2,u]=meshgrid(x_prime,ui);%(xxi,ui);
% %xx2 is actually xx, but changes in columns. u is the frequencies changing in
% %lines.

x1=x+xx/2;
x2=x-xx/2;

% equation (12) 'cross spectral density in far field'

sigma=3e-3; 
% x1-x2 direction; inverse propotional
% also proptionalto the u expansion in wigner


delta=5e-4; % x1+x2 direction; inverse propotional

% wavelength 2.25 micro, sigma=delta=5e-4

M=40; 

k=2*pi/lambda;

r=1e6;

theta1=0.1; theta2=0.1;

C=0;  

SUM=0;

for m=1:M
    
C_m=factorial(M)/factorial(m)/factorial(M-m)*(-1)^(m-1)/sqrt(m);

C=C+C_m;

a_mx=sigma*sqrt((2*m*delta^2+4*sigma^2)/(m*delta^2+4*sigma^2));

b_mx=sqrt(2/m/delta^2+1/sigma^2);

c_mx=k^2*sigma^2*m*delta^2/(m*delta^2+4*sigma^2);

d_mx=2*k^2*sigma^4/(m*delta^2+4*sigma^2);


SUM_m=C_m*a_mx/b_mx*exp(-c_mx*(x1.^2+x2.^2)-d_mx*(x1-x2).^2);

% SUM_m=C_m*a_mx/b_mx*exp(-c_mx*(x1+x2).^2-d_mx*(x1-x2).^2);

SUM=SUM+SUM_m;

end

Corr_rectang_GS=2*k^2*cos(theta1)*cos(theta2)/C^2/r^2.*SUM;

plot_2d(xi,'x=(x1+x2)/2',x_prime,'x prime=x1-x2',real(Corr_rectang_GS),'Cross spectral density (real)',18);

figure(22)
surf(xi,x_prime,abs(Corr_rectang_GS));shading interp


Phase=exp(-2i*pi/lambda*u.*xx2);
% corresponding-pixel multiplication 
% phase shift factor, horizontal axis is xx (rr), vertical axis is q.

Wig=real(Phase*Corr_rectang_GS*dx);

plot_2d(xi,'x (mm)',ui,'u (radian)',Wig,'Wigner of Rectangular GS',16);


%% equation (2) 'degree of coherence'

sigma=1000; delta=0.8; M=40; k=2*pi/lambda;

r=1e3;

theta1=0.1; theta2=0.1;

C=0;  

SUM=0;
%
for m=1:M
    
C_m=factorial(M)/factorial(m)/factorial(M-m)*(-1)^(m-1)/sqrt(m);

C=C+C_m;

% a_mx=sigma*sqrt((2*m*delta^2+4*sigma^2)/(m*delta^2+4*sigma^2));
% 
% b_mx=sqrt(2/m/delta^2+1/sigma^2);
% 
% c_mx=k^2*sigma^2*m*delta^2/(m*delta^2+4*sigma^2);
% 
% d_mx=2*k^2*sigma^4/(m*delta^2+4*sigma^2);

SUM_m=C_m*exp(-(x1-x2).^2/2/m/delta^2);

SUM=SUM+SUM_m;

end

I=1/C^2*SUM;

% plot(xi,I);

% figure(31)
plot_2d(xi,'x',x_prime,'x prime',real(I),'degree of coherence',18);


%%

sigma=10;

delta=1;

M=10;

lambda=0.633e-3;

k=2*pi/lambda;

z=100e3;

B=4*pi^4*z/3/lambda^2*B2;

B2=A_a*C_n^2/2/(alfa-2)*(kapa);

C=0; 

D=0;

for m=1:M
    
C_m=factorial(M)/factorial(m)/factorial(M-m)*(-1)^(m-1)/sqrt(m);

A_m=1/8/sigma^2+1/2/m/delta^2;

D_m=(-1)^(m-1)/sqrt(2*m*(sigma^2*k^2/2/z^2+A_m+B))*factorial(M)./factorial(m)./factorial(M-m)...
    .*exp(-(A_m*B*(x1-x2).^2+k^2/4/z^2*(x1+x2).^2+B*1i*k/z*(x2.^2-x1.^2))/4/(A_m+B))...
    .*exp(-(sigma^2*k^2*((2*A_m+B)*(x2-x1)-1i*k/2/z*(x1+x2))^2)./(8*z^2*(A_m+B)*(sigma^2*k/2/z^2+A_m+B)));

C=C+C_m;   

D=D+D_m;

end

Corr=(sigma*k/C/z)^2*exp(-1i*k/2/z*(x1.^2-x2.^2)).*exp(-3/4*B*(x1-x2).^2).*D;

plot_2d(xi,'x',x_prime,'x prime',real(Corr),'Correlation rectangular GS',18);

