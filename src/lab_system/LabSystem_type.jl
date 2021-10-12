################################
#  LAB SYSTEM TYPE DEFINITION  #
################################

"""
    mutable struct LabSystem

Object describing the current lab implementation.

# Fields

- `hamiltonian :: AbstractOperator`, hamiltonian
- `eigensys :: Dict{Symbol,Any}`, its eigensystem
- `dipole_hor :: AbstractOperator`,  complete dipole operator with horizontal outgoing polarization
- `dipole_ver :: AbstractOperator` complete dipole operator, vertical outgoing polarization
- `dipoles_hor :: Vector{DipoleOperator}`,list of horizontal dipole operators of the individual sites
- `dipoles_ver :: Vector{DipoleOperator}`, list of vertical dipole operators of the individual sites
- `sites :: Vector{CoordinateFrame}`, the individual sites as coordinate frames
- `sample :: CoordinateFrame`, the sample orientation as a coordinate frame
- `edge :: Int64`, the global edge
- `spin_quantization :: Symbol`, spin quantization frame (description)
- `q_in :: Vector{Float64}`, global incoming wavevector
- `epsilon_in :: Vector{Float64}`
- `q_out :: Vector{Float64}`, global outgoing wavevector
- `epsilon_out_hor :: Vector{Float64}`, horizontal polarization (relative to q_out)
- `epsilon_out_ver :: Vector{Float64}`, vertical polarization   (relative to q_out)
"""
mutable struct LabSystem

    # Hamiltonian
    hamiltonian :: AbstractOperator
    # eigensystem
    eigensys :: Dict{Symbol,Any}

    # complete dipole operator, horizontal outgoing polarization
    dipole_hor :: AbstractOperator
    # complete dipole operator, vertical outgoing polarization
    dipole_ver :: AbstractOperator

    # liste of horizontal dipole operators of the individual sites
    dipoles_hor :: Vector{DipoleOperator}
    # liste of vertical dipole operators of the individual sites
    dipoles_ver :: Vector{DipoleOperator}

    # the individual sites as coordinate frames
    sites :: Vector{CoordinateFrame}

    # the sample orientation as a coordinate frame
    sample :: CoordinateFrame

    # the global edge
    edge :: Int64
    # spin quantization frame (description)
    spin_quantization :: Symbol

    # global incoming wavevector
    q_in            :: Vector{Float64}
    epsilon_in      :: Vector{Float64}

    # global outgoing wavevector
    q_out           :: Vector{Float64}
    epsilon_out_hor :: Vector{Float64}  # horizontal polarization (relative to q_out)
    epsilon_out_ver :: Vector{Float64}  # vertical polarization   (relative to q_out)

end


function LabSystem(hamiltonian :: AbstractOperator)
    ls = LabSystem(
        hamiltonian,
        Dict{Symbol,Any}(),
        ZeroOperator(getT2GBasisXYZ()),
        ZeroOperator(getT2GBasisXYZ()),
        Vector{DipoleOperator}(undef, 1),
        Vector{DipoleOperator}(undef, 1),
        CoordinateFrame[CoordinateFrame() for s in 1:length(get_sites(basis(hamiltonian)))],
        CoordinateFrame(),
        3,
        :local,
        [1,0,1],
        [0,1,0],
        [1,0,-1],
        [0,1,0],
        [1,0,1]
    )
    recalculate_dipole_operators!(ls, new_objects=true)
    recalculate_hamiltonian!(ls, true, true)
    return ls
end




# custom show function
import Base.show
function show(io::IO, lab::LabSystem)
    print(io, "LabSystem object\n")
    # beam properties
    print(io, "Incoming beam parameters:\n")
    print(io, " -> q vector     = $(round.(lab.q_in, digits=3))\n")
    print(io, " -> polarization = $(round.(lab.epsilon_in, digits=3))\n")
    print(io, "Outgoing beam parameters:\n")
    print(io, " -> q vector       = $(round.(lab.q_out, digits=3))\n")
    print(io, " -> polarization h = $(round.(lab.epsilon_out_hor, digits=3))\n")
    print(io, " -> polarization v = $(round.(lab.epsilon_out_ver, digits=3))\n")
    # print the sites
    print(io, "\nSample alignment: (X=$(round.(lab.sample.X, digits=2)), Y=$(round.(lab.sample.Y, digits=2)), Z=$(round.(lab.sample.Z, digits=2)))\n")
    print(io, "\nFollowing sites are known:\n")
    for s in 1:length(lab.sites)
        print(io, " -> site $(s) @ r_$(s) = $(round.(lab.sites[s].position, digits=3)) within sample\n")
        print(io, "    coordinates w.r.t. sample: (X=$(round.(lab.sites[s].X, digits=2)), Y=$(round.(lab.sites[s].Y, digits=2)), Z=$(round.(lab.sites[s].Z, digits=2)))\n")
        print(io, "    Local beam parameters: $(round.(lab.dipoles_hor[s].q_in, digits=3)) -> $(round.(lab.dipoles_hor[s].q_out, digits=3))\n")
        print(io, "      -- > eps horizontal: eps_in || $(round.(lab.dipoles_hor[s].eps_in, digits=3)) -> eps_out || $(round.(lab.dipoles_hor[s].eps_out, digits=3))\n")
        print(io, "      -- > eps vertical:   eps_in || $(round.(lab.dipoles_ver[s].eps_in, digits=3)) -> eps_out || $(round.(lab.dipoles_ver[s].eps_out, digits=3))\n")
    end
    # print the hamiltonian
    print(io, "\nCurrent hamiltonian is given as ")
    show(io, lab.hamiltonian)
end