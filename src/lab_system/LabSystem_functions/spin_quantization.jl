function set_spin_quantization_to_sample_z!(lab :: LabSystem)
    lab.spin_quantization = :sample
    for s in 1:length(get_sites(basis(lab)))
        set_parameter!(lab, :spin_quantization, lab.sites[s], site=s, print_result=false)
    end
    recalculate!(lab, true, true)
    return nothing
end
function set_spin_quantization_to_local_z!(lab :: LabSystem)
    lab.spin_quantization = :local
    for s in 1:length(get_sites(basis(lab)))
        set_parameter!(lab, :spin_quantization, CoordinateFrame(lab.sites[s].position), site=s, print_result=false)
    end
    recalculate!(lab, true, true)
    return nothing
end