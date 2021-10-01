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