using ITensors
using ITensorTDVP
using PyCall
np = pyimport("numpy")

# ==============================================================================
# TDVP FOR ISING
# ==============================================================================

function H_ising(sites, h)
    N = length(sites)    
    os = OpSum()
    for j in 1:(N - 1)
        os += -1, "X", j, "X", j + 1
        os += -h, "Z", j
    end
    os += -1, "X", 1, "X", N
    os += -h, "Z", N
    return MPO(os, sites)
end

function Sx(n) #
    os = OpSum()
    for j in 1:n
        os += +1,"Sx", j
    end
    return MPO(os, sites)
end

N = 50
sites = siteinds("S=1/2", N)
psi0 = productMPS(sites, n -> "X+")


dt = 0.1
Nsteps = 1500
χ = 50
h = np.linspace(0,4,30)
Sz_t = zeros(Nsteps) #Stores the total magnetization at each time step
Sz_av = zeros(length(h)) #Stores the averaged magnetization values

tdvp_kwargs = (time_step=-1im*dt, normalize=false, maxdim=χ, cutoff=1.0e-10, solver_backend="exponentiate")


for i in 1:length(h)
    
    psi_previous = psi0
    H = H_ising(sites, h[i])
    
    for t in 1:Nsteps
    
        psi_current = tdvp(H, psi_previous, -1im*dt; tdvp_kwargs)
        
        Sz_t[t] = real(inner(psi_current', Sx(N), psi_current))
        
        psi_previous = psi_current
    end
    
    Sz_av[i] = (1/(Nsteps*dt))*np.trapz(Sz_t, np.linspace(0, Nsteps*dt, Nsteps))
        
end

np.savetxt("Magnetization_Ising_N100.txt", Sz_av)



