using Plots
using QuadGK
using SpecialFunctions
using Printf

c = 3e8;
f = 2e9;
ϵᵣ = 3.2;
Zₛ = 50;
Zₗ = 100;
Z₀ = 50;
MaxRL = -20;
numSections = 16;
Γₘ = 10^(MaxRL/20);
ρ₀ = log(Zₗ/Zₛ)/2;
A = acosh(ρ₀/Γₘ);
λ = c/f; 
λeff = λ/sqrt(ϵᵣ);
β = 2*π/λeff;
L = A/β; #meter
@printf("L = %g cm",L*100)
@printf("\n")
@printf("L = %g mils",L*1e6/25.4)
l  = L/numSections;
z = range(0,L,numSections);
I(y)=besseli(1,A*sqrt(1-y.^2))./(A*sqrt(1-y.^2));
ϕ = quadgk(I,0,2*z/L-1,rtol=1e-8) 
#Z = exp(log(ZL*ZS)/2+ρ₀*A^2*ϕ/cosh(A));