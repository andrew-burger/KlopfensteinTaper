using Plots
using QuadGK
using SpecialFunctions
using Printf

#initilization
C = 3e8;
f = 2e9; #lowest frequency for taper to operate at
startFreq = 0; #starting frequency for analysis/plot
freqStep = 10e6; #frequency step for analysis/plot
stopFreq = 50e9; #highest frequency for analysis/plot
ϵᵣ = 3.2;
Zₛ = 50;
Zₗ = 100;
Z₀ = 50;
MaxRL = -20;
numSections = 16;

Γₘ = 10^(MaxRL/20);
ρ₀ = log(Zₗ/Zₛ)/2;
A = acosh(ρ₀/Γₘ);
λ = C/f; 
λeff = λ/sqrt(ϵᵣ);
β = 2*π/λeff;
L = A/β; #meter
@printf("L = %g cm",L*100)
@printf("\n")
@printf("L = %g mils",L*1e6/25.4)
@printf("\n")
l  = L/numSections;
x = range(0,L,numSections);
function I(y)
    besseli(1,A*sqrt(1-y.^2))./(A*sqrt(1-y.^2));
end;
function ϕ(xsym)
    results, error = quadgk(I,0,2 .*xsym /L .- 1,rtol=1e-8);
    return results
end;
Z = exp.(log(Zₗ*Zₛ)/2 .+ ρ₀*A^2*ϕ.(x)/cosh(A));
freq = startFreq:freqStep:stopFreq;
β = 2π*freq/C*sqrt(ϵᵣ);
Γ = ρ₀*exp.(-1im*β*L).*cos.(sqrt.(complex((β*L).^2 .- A^2)))/cosh(A);
plot(freq,abs(Γ));