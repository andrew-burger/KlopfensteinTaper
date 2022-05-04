format shortG
clear
close all
c = 3e8;
%%
f = 2e9; %lowest freq
er = 3.2;
ZS = 50;
ZL = 100;
Z0 = 50;
MaxRL = -20;
numSections = 16;
%%
GammaMax = 10^(MaxRL/20);
rho0 = log(ZL/ZS)/2;
A = acosh(rho0/GammaMax)
lambda = c/f;
lambdaeff = lambda/sqrt(er);
Beta = 2*pi/lambdaeff;
L = A/Beta; %meter
l  = L/numSections;
z = linspace(0,L,numSections);
syms y
phi = NaN(1,length(z));
for x = 1:length(z)
    phi(x) = vpaintegral(besseli(1,A*sqrt(1-y.^2))./(A*sqrt(1-y.^2)),y,0,2*z(x)/L-1);
end
Z = exp(log(ZL*ZS)/2+rho0*A^2*phi/cosh(A));
Z = Z'
freq = 0:10e6:50e9;
Beta = 2*pi*freq/c*sqrt(er);
Gamma = rho0*exp(-1j*Beta*L).*cos(sqrt((Beta*L).^2-A^2))/cosh(A);
figure
plot(freq,abs(Gamma))
for i = 1:length(freq)
    for j = 1:numSections
        ABCDCube(:,:,i,j) = [cos(Beta(i)*l) 1j*Z(j)*sin(Beta(i)*l); 1j*1/Z(j)*sin(Beta(i)*l) cos(Beta(i)*l)];
    end
end
for i = 1:numSections-1
    ABCDCube(:,:,:,i+1) = pagemtimes(ABCDCube(:,:,:,i),ABCDCube(:,:,:,i+1));
end
ABCDCube = squeeze(ABCDCube(:,:,:,numSections));
SParam  = [(ABCDCube(1,1,:)+ABCDCube(1,2,:)/Z0-ABCDCube(2,1,:)*Z0-ABCDCube(2,2,:))./(ABCDCube(1,1,:)+ABCDCube(1,2,:)/Z0+ABCDCube(2,1,:)*Z0+ABCDCube(2,2,:)) ...
     2*(ABCDCube(1,1,:).*ABCDCube(2,2,:)-ABCDCube(1,2,:).*ABCDCube(2,1,:))./(ABCDCube(1,1,:)+ABCDCube(1,2,:)/Z0+ABCDCube(2,1,:)*Z0+ABCDCube(2,2,:));
    2./(ABCDCube(1,1,:)+ABCDCube(1,2,:)/Z0+ABCDCube(2,1,:)*Z0+ABCDCube(2,2,:)) ...
     (-ABCDCube(1,1,:)+ABCDCube(1,2,:)/Z0-ABCDCube(2,1,:)*Z0+ABCDCube(2,2,:))./(ABCDCube(1,1,:)+ABCDCube(1,2,:)/Z0+ABCDCube(2,1,:)*Z0+ABCDCube(2,2,:))];
GammaL = (100-50)./(100+50);
GammaIn = squeeze(SParam(1,1,:)+SParam(1,2,:).*(SParam(2,1,:)*GammaL./(1-SParam(2,2,:).*GammaL)));
figure
plot(freq,abs(GammaIn))