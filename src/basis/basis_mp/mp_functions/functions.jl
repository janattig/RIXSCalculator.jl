function printMPState(
            state :: Vector{<:Number},
            b     :: MPBasis{N,SPBS}
            ;
            cutoff :: Real = 0.001,
            digits :: Integer = 3,
            max_contributions :: Integer = 10,
            energy :: Union{Symbol, Real} = :None,
            es_index :: Union{Symbol, Integer} = :None
        ) where {N, SPBS<:AbstractSPBasisState}

    # print the index
    if es_index != :None
        print("Multi-Particle Eigenstate $(es_index)")
    else
        print("Multi-Particle state")
    end
    # print the energy
    if energy != :None
        println(" with energy $(round(energy, digits=digits)):")
    else
        println(":")
    end
    # calculate basis state contributions
    max_weight = sum([abs(v)^2 for v in state])
    contrib = [(i, state[i], abs(state[i])^2/max_weight) for i in 1:length(state) if abs(state[i])^2/max_weight>cutoff]
    contrib = sort(contrib, by=c->-c[3])
    # print the largest contributions
    println("State has contributions from $(length(contrib)) of $(length(b)) basis state(s) with amplitude(s) larger then $(round(cutoff*100,digits=2))% of total weight:")
    # print all contributions
    ci = 0
    for c in contrib
        bs = "(+)"
        for sps in b[c[1]].occupation
             bs *= " d_" * summary(b.single_particle_basis[sps])
        end
        bs *= " |vac⟩"
        println(bs, " * (", round.(c[2], digits=3), ")\t(", round(c[3]*100, digits=2), "%)")
        ci += 1
        if ci == max_contributions && length(contrib) != max_contributions
            println("... (+ $(length(contrib)-ci) more)")
            break
        end
    end
end
export printMPState








# function to get a multiparticle basis
# for N particles and a given single particle basis
"""
    getMultiParticleBasis(sp_basis :: SPBasis{SPBS}, N :: Integer  :: MPBasis{N<0 ? length(sp_basis)+N : N,SPBS} where {SPBS<:AbstractSPBasisState}

This function returns the multiparticle basis given a single particle basis `sp_basis` and a number `N` of particles.
"""
function getMultiParticleBasis(
            sp_basis :: SPBasis{SPBS},
            N :: Integer
        ) :: MPBasis{N<0 ? length(sp_basis)+N : N,SPBS} where {SPBS<:AbstractSPBasisState}

    # maybe invert N
    if N<0
        N = length(sp_basis) + N
    end

    # get all permutations of length N of numbers between 1 and the size of the sp basis
    perms = getPositivePermutations(N, length(sp_basis))

    # create multiparticle states out of the permutations
    mp_states = map(p->MPBasisState{N}(p), perms)

    # define a new multiparticle basis
    # since lookup calculation is included in this, can return it directly
    return MPBasis{N,SPBS}(
        mp_states,
        sp_basis
    )
end
export getMultiParticleBasis
