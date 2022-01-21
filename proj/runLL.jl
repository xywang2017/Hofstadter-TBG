using PyPlot
using Printf
using DelimitedFiles
using JLD
fpath = pwd()
include(joinpath(fpath,"libs/HofstadterLL_modv3.jl"))

##
w1=96.056
w0=0.0
params = Params(dθ=1.05π/180,w0=w0,w1=w1,ϵ=0.0)  #chiral limit
qs = collect(2:18)
# qs = [4]
p = 1
data = Dict()
σz = []
for iq in eachindex(qs)
    println("q=$(qs[iq])")
    q = qs[iq]
    lk = 10 
    if (q>=10)
        lk = 4 
    end
    if (q>20)
        lk = 1
    end
    hof = constructHofstadterLL(params,q=q,p=p,nLL=25q,lk=lk)
    data["$iq"] =  hof.spectrum[:]
    push!(σz,hof.σz)
end
save("BMresults/theta1.05_LL_w00.0.jld","data",data,"σz",σz)

##
fig = figure(figsize=(3,4))
qs = collect(2:18)
# d=load("Angle1.05/LL_results.jld")
# data = d["data"]
# σz = d["σz"]
for iq in eachindex(qs)
    p=1
    q = qs[iq]
    lk = 10 
    if (q>=10)
        lk = 4 
    end
    if (q>20)
        lk = 1
    end
    l1,l2 = lk, lk
    # println(length(data["$iq"]/(2q*l1*l2)))
    ϵ = sort(data["$iq"],by=abs)
    ϵflat = ϵ[1:(2q*l1*l2)]
    ϵout = ϵ[(2q*l1*l2+1):(4q*l1*l2)]
    ϕflat = ones(size(ϵflat)) * p/q
    ϕout = ones(size(ϵout)) * p/q
    plot(ϕflat,ϵflat,".",c="b",ms=1.5,markeredgecolor="none")
    plot(ϕout,ϵout,".",c="gray",ms=2,markeredgecolor="none")
end

# ylabel(L"$\frac{ϵ}{ϵ_0}$",rotation=0)
xlabel(L"ϕ/ϕ_0")
xticks([0,0.1,0.2,0.3,0.4,0.5])
ylim([-0.55,0.55])
yticks([-0.4,0,0.4])
tight_layout()
display(fig)
savefig("LL_spectrum.png",dpi=500,transparent=true)
close(fig)

##
# here I show the sublattice polarization as a function of ϕ/ϕ0
fig = figure(figsize=(3,3))

d=load("BMresults/theta1.05_LL_w00.0.jld")
data = d["data"]
σz = d["σz"]
ϕ = p./qs
nσz = [sum(σz[i])/length(σz[i]) for i in eachindex(σz)]
plot(ϕ,nσz,"ro",markerfacecolor="none",ms=3,label=L"w_0/w_1=0.0")

d=load("BMresults/theta1.05_LL_w00.7.jld")
data = d["data"]
σz = d["σz"]
ϕ = p./qs
nσz = [sum(σz[i])/length(σz[i]) for i in eachindex(σz)]
plot(ϕ,nσz,"bo",markerfacecolor="none",ms=3,label=L"w_0/w_1=0.7")

xlabel(L"ϕ/ϕ_0")
ylabel(L"Tr⟨σ_z⟩")
legend(loc="right")
ylim([0,2.2])
tight_layout()
savefig("projected_sigmaz.pdf",transparent=true)
display(fig)
close(fig)

##
# params = Params(dθ=1.38π/180);
# hof = constructHofstadterLL(params,q=3,p=1,nLL=40,lk=128);

# kmesh = ( reshape(hof.k1,:,1)*params.g1 .+ reshape(hof.k2,1,:)*params.g2 ) / abs(params.g2);
# kx = real(kmesh);
# ky = imag(kmesh);

# ##
# fig,ax = subplots(2,2,figsize=(8,4))
# cnt = -1
# for r in 1:2 
#     for c in 1:2
#         pl = ax[r,c].contourf(kx,ky,hof.spectrum[hof.nH+cnt,:,:],cmap="Blues_r")
#         colorbar(pl,ax=ax[r,c])
#         ax[r,c].axis("equal")
#         cnt = cnt + 1
#     end
# end
# tight_layout()
# savefig("Angle1.38/lowest_LLsq3.pdf")
# display(fig)
# close(fig)

# ## 
# # cut 
# fig = figure(figsize=(4,3)) 
# maxLL = 2
# ϵcut = hof.spectrum[(hof.nH-maxLL):(hof.nH+maxLL+1),1:(hof.lk÷2+1),24]
# ϵcut = [ϵcut  hof.spectrum[(hof.nH-maxLL):(hof.nH+maxLL+1),hof.lk÷2+1,25:end]]

# for n in 1:size(ϵcut,1)
#     plot(ϵcut[n,:],"-",lw=0.5)
# end
# tight_layout()
# savefig("Angle1.38/lowest_LLs_cutq3.pdf")
# display(fig)
# close(fig)
