"""
    set_spin_quantization_to_sample_z!(lab :: LabSystem)

The function sets the spin quantization axis of the object `lab` to the z-axis of the sample.
"""
function set_spin_quantization_to_sample_z!(lab :: LabSystem)
    lab.spin_quantization = :sample
    for s in 1:length(get_sites(basis(lab)))
        set_parameter!(lab, :spin_quantization, lab.sites[s], site=s, print_result=false)
    end
    recalculate!(lab, true, true)
    return nothing
end
"""
    set_spin_quantization_to_local_z!(lab :: LabSystem)

The function sets the spin quantization axis of the object `lab` to the local z-axis.
"""
function set_spin_quantization_to_local_z!(lab :: LabSystem)
    lab.spin_quantization = :local
    for s in 1:length(get_sites(basis(lab)))
        set_parameter!(lab, :spin_quantization, CoordinateFrame(lab.sites[s].position), site=s, print_result=false)
    end
    recalculate!(lab, true, true)
    return nothing
end
export set_spin_quantization_to_sample_z!, set_spin_quantization_to_local_z!