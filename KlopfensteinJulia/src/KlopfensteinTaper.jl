using Plots
using QuadGK
using SpecialFunctions
using Printf

default(show=true)

#initilization
C = 299792458;
fₗ = 13e9; #lowest frequency for taper to operate at
fₕ = 40e9; # highest frequency taper needs to work at
fₑ = fₕ+fₗ; # sampling frequency margin due to 'lobe'
startFreq = 0; #starting frequency for analysis/plot
freqStep = 10e6; #frequency step for analysis/plot
stopFreq = 50e9; #highest frequency for analysis/plot
ϵᵣ = 3.2;
Zₛ = 50; #source impedance
Zₗ = 100; #load impedance, must be larger than Zₛ
Z₀ = Zₛ;
MaxRL = -40; #maximum allowable return loss

Γₘ = 10^(MaxRL/20);
ρ₀ = log(Zₗ/Zₛ)/2;
A = acosh(ρ₀/Γₘ);
λ = C/fₗ; 
λeff = λ/sqrt(ϵᵣ);
β = 2*π/λeff;
L = A/β; #meter
numSections = ceil(Int,L/(C/fₑ/sqrt(ϵᵣ)/2));
@printf("# of sections: %g\n",numSections)
@printf("L = %g cm\n",L*100)
@printf("L = %g mils\n",L*1e6/25.4)
l  = L/numSections;
x = range(0+l/2,L-l/2,numSections);
function I(y)
    besseli(1,A*sqrt(1-y.^2))./(A*sqrt(1-y.^2));
end;
function ϕ(xsym)
    results, error = quadgk(I,0,2 *xsym /L - 1,rtol=1e-16);
    return results
end;
Z = similar(x);
@. Z = exp(log(Zₗ*Zₛ)/2 + ρ₀*A^2*ϕ(x)/cosh(A));
freq = startFreq:freqStep:stopFreq;
numFreq = length(freq);
β = similar(freq);
@. β = 2π/(C/freq/sqrt(ϵᵣ));
Γ = Vector{ComplexF64}(undef,numFreq);
@. Γ = ρ₀*exp(-1im*β*L)*cos(sqrt(complex((β*L)^2 - A^2)))/cosh(A);
Γmag = abs.(Γ);
display(plot(freq,Γmag));
println("Impedance of each section (Ω):")
println.(Z);

A = complex(ones(numFreq,numSections));
@. A = cos(β*l)*Z'/Z';
A = reshape(A,1,numFreq,numSections);
B = complex(ones(numFreq,numSections));
@. B = 1im*sin(β*l)*Z';
B = reshape(B,1,numFreq,numSections);
C = complex(ones(numFreq,numSections));
@. C = 1im*sin(β*l)./Z';
C = reshape(C,1,numFreq,numSections);
#D = A; %this is given so D will not be calculated.
temp = [A;C;B;A]
ABCD = reshape(temp,2,2,numFreq,numSections)
for i in 1:numFreq
    for j in 1:numSections-1
        ABCD[:,:,i,j+1] .= ABCD[:,:,i,j]*ABCD[:,:,i,j+1];
    end
end
ABCD = ABCD[:,:,:,numSections]
A = ABCD[1,1,:];
B = ABCD[1,2,:];
C = ABCD[2,1,:];
D = ABCD[2,2,:];
Δ = complex(ones(numFreq));
@. Δ = A+B/Zₛ+C*Zₛ+D;
S11 = complex(ones(numFreq));
S21 = complex(ones(numFreq));
S12 = complex(ones(numFreq));
S22 = complex(ones(numFreq));
@. S11 = (A+B/Zₛ-C*Zₛ-D)/Δ;
@. S12 = 2*(A*D-B*C)/Δ;
@. S21 = 2/Δ;
@. S22 = (-A+B/Zₛ-C*Zₛ+D)/Δ;
Γₗ = (Zₗ-Zₛ)/(Zₗ+Zₛ);
Γᵢ = complex(ones(numFreq));
@. Γᵢ = S11+S12*S21*Γₗ/(1-S22*Γₗ);
Γᵢmag = abs.(Γᵢ);
display(plot(freq,Γᵢmag));
dB = 20 .* log10.(Γᵢmag);
display(plot(freq,dB));
ylims!(-40, 0)