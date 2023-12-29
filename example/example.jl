# Implementation of the Peng-Robinson (PR) equation of state for methane
using DelimitedFiles, Roots
using PyCall, PyPlot
using EntropyScaling

# Read experiemntal data
header = split(readline("example/methane_vis.csv"),',')
dat = readdlm("example/methane_vis.csv",',',skipstart=1)
T_exp = Float64.(dat[:,findfirst(header.=="T")])                # K
p_exp = Float64.(dat[:,findfirst(header.=="p")])                # MPa
η_exp = Float64.(dat[:,findfirst(header.=="eta")])              # Pa s
state = dat[:,findfirst(header.=="state")]                      # -
ref = dat[:,findfirst(header.=="shortref")]                     # -

# Parameters for methane
Tc = 190.4                                                      # K
pc = 4.6e6                                                      # Pa
ω = 0.011                                                       # -
M = 0.016043                                                    # kg mol⁻¹

# Load PR equation of state functions
include("PR.jl")
(p_PR, s_PR, B_PR) = PR(Tc, pc, ω, M)

# Implicit calculation of density
ϱ_exp = Float64[]
for (i,state_i) in enumerate(state)
    rts = find_zeros(x -> p_PR(T_exp[i],x,1.0) - p_exp[i]*1e6,0,M/(0.07780*EntropyScaling.R*Tc/pc ))
    if state_i in ["G","SCR"]
        push!(ϱ_exp,rts[1])
    else
        push!(ϱ_exp,rts[end])
    end
end

# Fit entropy scaling
α_CH₄, ηˢ, s = fit_entropy_scaling(T_exp, ϱ_exp, η_exp, "vis"; sfun=s_PR, Bfun=B_PR, Tc=Tc, pc=pc, M=M)

# Plot results
figure(figsize=(6,4))
xs = [0.:0.01:1.1*maximum(s);]
plot(xs,EntropyScaling.fun_es(xs,α_CH₄;prop="vis"),"k-",zorder=0)
scat=scatter(s,log.(ηˢ),c=T_exp,marker="o",s=15,zorder=10)
xlabel(L"\tilde{s}_{\rm conf}",fontsize=10)
ylabel(L"\ln\left(\widehat{\eta}^{+}\right)",fontsize=10)
xlim(minimum(xs),maximum(xs))
cb = colorbar(scat)
cb.ax.set_title(L"T\,/\,{\rm K}",fontsize=10)
tight_layout()

# Calculation of the viscosity along isotherms
what_ref = ref .== "Hellemans et al. (1970) [A]" .&& T_exp .< 140.0
T_iso = unique(T_exp[what_ref])

# Plot experimental data
figure(figsize=(6,4))
scatter(p_exp[what_ref],η_exp[what_ref].*1e3,marker="o",c=T_exp[what_ref],facecolor="white",s=15.0)
cb = colorbar()
cb.ax.set_title(L"T\,/\,{\rm K}",fontsize=10)

for Ti in T_iso
    what_i = what_ref .&& T_exp .== Ti
    pi = LinRange(0.01,11.0,20)
    coli = get_cmap()((Ti-minimum(T_iso))/(maximum(T_iso)-minimum(T_iso)))
    
    # Implicit calculation of density
    ϱi = Float64[]
    for pj in pi push!(ϱi,find_zeros(x -> p_PR(Ti,x,1.0) - pj*1e6,0,M/(0.07780*EntropyScaling.R*Tc/pc ))[end]) end
    
    # Calculate and plot isotherm
    η_i = call_entropy_scaling(repeat([Ti], length(pi)), ϱi, [α_CH₄], "vis"; sfun=s_PR, Bfun=[B_PR], Tc=[Tc], pc=[pc], M=[M])
    plot(pi,η_i.*1e3,"-",c=coli,linewidth=1.0)
end
plot(NaN,NaN,"-k",label="Entropy scaling")
scatter(NaN,NaN,marker="o",s=15,c="k",label="Hellemans et al. (1970)")
legend(loc="upper left",fontsize=10,frameon=false)
xlim(0,11)
ylim(0.06,0.22)
xlabel(L"p\;/\;{\rm MPa}",fontsize=10)
ylabel(L"\eta\;/\;{\rm mPa\,s}",fontsize=10)
tight_layout()