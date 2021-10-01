################################################################################
#
#   Spectrum generating functions
#
################################################################################


# Calculate the dipole amplitude when transitioning from -> D -> to
# which corresponds to <to|D|from>
function get_amplitude(
        dipole_operator :: AbstractOperator{B},
        state_from :: Vector{<:Number},
        state_to   :: Vector{<:Number}
    ) ::Complex{Float64} where {B <: AbstractBasis}

    # build amplitude like <to|D|from>
    return dot(state_to,  matrix_representation(dipole_operator) * state_from)
end

# Calculate the dipole amplitude when transitioning i_from -> D -> i_to
# which corresponds to <i_to|D|i_from>
function get_amplitude(
        eigensys        :: Dict{Symbol, Any},
        dipole_operator :: AbstractOperator{B},
        i_from :: Int64,
        i_to   :: Int64
    ) ::Complex{Float64} where {B <: AbstractBasis}

    # build amplitude like <to|D|from>
    return get_amplitude(dipole_operator, eigensys[:vectors][i_from], eigensys[:vectors][i_to])
end


# Calculate a spectrum based on a RIXSSample object
# polarizations are in global coordinates (NOT in sample coordinates)
# basically a wrapper around getSingleSiteTransitionIntensities(sample,...)
function get_spectrum(
        eigensys        :: Dict{Symbol, Any},
        dipole_operator :: AbstractOperator{B}
        ;
        linewidth :: Real = 25.0
    ) :: Spectrum where {B <: AbstractBasis}

    # allocate the transitions
    transitions = Vector{Transition}(undef, length(eigensys[:values]))

    # fill all transitionss
    for i in 1:length(transitions)
        transitions[i] = Transition(
            abs( get_amplitude(eigensys, dipole_operator, 1, i) )^2,
            linewidth,
            eigensys[:values][i] - eigensys[:values][1]
        )
    end

    # construct spectrum
    return Spectrum(transitions)
end

# get spectrum from hamiltonian with ground states in array states
function get_spectrum(
        eigensys        :: Dict{Symbol, Any},
        dipole_operator :: AbstractOperator{B},
        states          :: Vector{<:Integer}
        ;
        linewidth :: Real = 25.0
    ) :: Spectrum where {B <: AbstractBasis}

    return sum([
        begin
            # allocate the transitions
            transitions = Vector{Transition}(undef, length(eigensys[:values]))

            # fill all transitionss
            for i in 1:length(transitions)
                transitions[i] = Transition(
                    (1.0 ./ length(states)) * abs( dot(eigensys[:vectors][i],  matrix_representation(dipole_operator) * eigensys[:vectors][j]) )^2,
                    linewidth,
                    eigensys[:values][i] - eigensys[:values][j]
                )
            end

            # construct spectrum
            Spectrum(transitions)
        end

        for j in states
    ])
end

function get_spectrum(
        eigensys        :: Dict{Symbol, Any},
        dipole_operator :: AbstractOperator{B},
        GS_energy       :: Real,
        args...
        ;
        precision :: Real = 1e-8,
        kwargs...
    ) :: Spectrum where {B <: AbstractBasis}

    # find out the states of that particular energy
    states = Int64[i for i in 1:length(eigensys[:values]) if isapprox(eigensys[:values][i], GS_energy, atol=precision)]

    # construct spectrum
    return get_spectrum(
            eigensys,
            dipole_operator,
            states,
            args...
            ;
            kwargs...
        )
end


# get spectrum from hamiltonian instead of eigensystem
function get_spectrum(
        hamiltonian     :: AbstractOperator{B},
        args...
        ;
        kwargs...
    ) :: Spectrum where {B <: AbstractBasis}

    # build eigensystem and pass on
    return get_spectrum(
        eigensystem(hamiltonian),
        args...
        ;
        kwargs...
    )
end

# export the function
export getSpectrum



################################################################################
#
#   Spectrum information functions
#
################################################################################

# compare the heights of peaks in a certain energy range
function comparePeakHeights(
        spectrum   :: S,
        energy_min :: Real,
        energy_max :: Real,
        discretization :: Integer = 10000
    ) where {S<:AbstractSpectrum}

    # construct the energy range to search
    energy_range = range(energy_min, stop=energy_max, length=discretization)
    intensities  = Float64[intensity(spectrum, energy) for energy in energy_range]

    # search the whole energy range for for points who are higher than their neighbors and construct
    # list of peak positions (list integers)
    peakpositions_index = Int64[]
    for i in 2:discretization-1
        if intensities[i] > intensities[i+1] && intensities[i] > intensities[i-1]
            push!(peakpositions_index, i)
        end
    end

    # list of actual positions
    peakpositions = Float64[]

    # filter again
    for i in peakpositions_index
        energy_range_sub = range(energy_range[i-1], stop=energy_range[i+1], length=discretization)
        intensities_sub  = Float64[intensity(spectrum, energy) for energy in energy_range_sub]
        for i in 2:discretization-1
            if intensities_sub[i] > intensities_sub[i+1] && intensities_sub[i] > intensities_sub[i-1]
                push!(peakpositions, energy_range_sub[i])
                break
            end
        end
    end

    # print all peak data
    println("Found ", length(peakpositions), " peaks:")
    for (i,p) in enumerate(peakpositions)
        println("Peak ", i, " @ ", p, " with intensity ", intensity(spectrum, p))
    end
    for i1 in 1:length(peakpositions)
    for i2 in 1:length(peakpositions)
        if i1==i2
            continue
        end
        println("I",i1,"/I",i2," = ", intensity(spectrum, peakpositions[i1])/intensity(spectrum, peakpositions[i2]))
    end
    end

end

# export the function
export comparePeakHeights





function identify_peaks(
        spectrum,
        E_min,
        E_max,
        N_grid_points,
        refinements = 3
    )

    # make a grid of energies
    E_vals = range(E_min, stop=E_max, length=N_grid_points)
    dE = E_vals[2] - E_vals[1]
    I_vals = Float64[intensity(spectrum, E) for E in E_vals]

    # make al ist of peaks
    peak_list = Float64[]

    # search for peaks
    for i in 2:N_grid_points-1
        if I_vals[i] > I_vals[i+1] && I_vals[i] > I_vals[i-1]
            push!(peak_list, E_vals[i])
        end
    end

    # maybe refine
    if refinements > 0
        peak_list_unrefined = peak_list
        peak_list = Float64[]
        for p in peak_list_unrefined
            append!(peak_list, identify_peaks(
                spectrum,
                p-dE,p+dE,
                N_grid_points,
                refinements-1
            ))
        end
    end

    # return the peak list
    return sort(unique(peak_list))

end