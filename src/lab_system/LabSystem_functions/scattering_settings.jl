# reset the incoming polarization
# eps_in is perpendicular to q_in and to the axis given (by default the Y axis)
function reset_polarization_in!(lab::LabSystem, axis_perpendicular::Vector{<:Real} = [0,1,0])
    # build the cross product of the given axis and q_in
    eps_direction  = cross(lab.q_in, axis_perpendicular)
    # normalize
    eps_direction  = eps_direction / norm(eps_direction)
    # set as polarization
    lab.epsilon_in = eps_direction
    # return nothing
    return nothing
end
export reset_polarization_in!

# reset the outgoing polarizations
# eps_out_horizontal is perpendicular to q_out and to the axis given (by default the Y axis)
# eps_out_vertical is perpendicular to q_out and to eps_out_horizontal
function reset_polarization_out!(lab::LabSystem, axis_perpendicular::Vector{<:Real} = [0,1,0])
    # build the cross product of the given axis and q_out to determine horizontal polarization
    eps_direction       = cross(lab.q_out, axis_perpendicular)
    # normalize
    eps_direction       = eps_direction / norm(eps_direction)
    # set as horizontal outgoing polarization
    lab.epsilon_out_hor = eps_direction
    # build the cross product of the given axis and q_out to determine vertical polarization
    eps_direction       = cross(lab.q_out, lab.epsilon_out_hor)
    # normalize
    eps_direction       = eps_direction / norm(eps_direction)
    # set as horizontal outgoing polarization
    lab.epsilon_out_ver = eps_direction
    # return nothing
    return nothing
end
export reset_polarization_out!



# set the incoming q vector
function set_q_in!(lab::LabSystem, q_in::Vector{<:Real}, reset_polarization::Bool=true)
    # set the q vector
    lab.q_in = q_in
    # maybe reset polarization
    if reset_polarization
        reset_polarization_in!(lab)
    end
    # return nothing
    return nothing
end
# set the incoming q vector as well as the polarization perp to the axis given
function set_q_in!(lab::LabSystem, q_in::Vector{<:Real}, axis_perpendicular::Vector{<:Real})
    # set the q vector
    lab.q_in = q_in
    # reset polarization
    reset_polarization_in!(lab, axis_perpendicular)
    # return nothing
    return nothing
end

# set the outgoing q vector
function set_q_out!(lab::LabSystem, q_out::Vector{<:Real}, reset_polarization::Bool=true)
    # set the q vector
    lab.q_out = q_out
    # maybe reset polarization
    if reset_polarization
        reset_polarization_out!(lab)
    end
    # return nothing
    return nothing
end
# set the outgoing q vector as well as the polarization perp to the axis given
function set_q_out!(lab::LabSystem, q_out::Vector{<:Real}, axis_perpendicular::Vector{<:Real})
    # set the q vector
    lab.q_out = q_out
    # reset polarization
    reset_polarization_out!(lab, axis_perpendicular)
    # return nothing
    return nothing
end

# export them
export set_q_in!, set_q_out!


# set the outgoing q vector as well as the polarization perp to the axis given
function setup_dQ!(lab::LabSystem, dq::Vector{<:Real}, eps_in::Vector{<:Real}, q_beam :: Real)
    # setup sample
    ls.sample.X = [1,0,0]
    ls.sample.Y = [0,1,0]
    ls.sample.Z = [0,0,1]
    # find out q direction
    dq_dir = dq ./ norm(dq)
    # find out perp to scattering plane
    n = cross(dq_dir, eps_in)
    n = n ./ norm(n)
    # find out q perp direction
    q_perp = cross(dq_dir, n)
    q_perp = q_perp ./ norm(q_perp)
    # find out contributions in the directions
    dq_perp = sqrt(q_beam^2 - (norm(dq)/2)^2)
    dq_par  = norm(dq)/2
    # build up q in and q out
    q_in  =  dq_par.*dq_dir .+ dq_perp.*q_perp
    q_out = -dq_par.*dq_dir .+ dq_perp.*q_perp
    # set the q vectors in the lab frame
    lab.q_in  = q_in
    lab.q_out = q_out
    # reset polarization
    reset_polarization_in!(lab, n)
    reset_polarization_out!(lab, n)
    # return nothing
    return nothing
end
export setup_dQ!







# set the beams depending on THETA and TWO THETA angles in geometry
function set_scattering_angles!(
            lab             :: LabSystem,
            angle_theta     :: Real,
            angle_two_theta :: Real,
            dQ              :: Real = 1.0,
            axis_perpendicular::Vector{<:Real} = [0,1,0]
        )
    # calculate q's
    q_in  = [-cos(angle_theta)                  , 0,  sin(angle_theta)                  ] .* dQ/(sin(angle_theta) + sin(angle_two_theta - angle_theta))
    q_out = [-cos(angle_two_theta - angle_theta), 0, -sin(angle_two_theta - angle_theta)] .* dQ/(sin(angle_theta) + sin(angle_two_theta - angle_theta))
    # pass to lab system
    set_q_in!(lab, q_in, axis_perpendicular)
    set_q_out!(lab, q_out, axis_perpendicular)
    # return nothing
    return nothing
end
# set the beams depending on THETA and TWO THETA angles in geometry
function set_scattering_angles_deg!(
            lab             :: LabSystem,
            angle_theta     :: Real,
            angle_two_theta :: Real,
            dQ              :: Real = 1.0,
            axis_perpendicular::Vector{<:Real} = [0,1,0]
        )
    # pass to the normal function
    set_scattering_angles!(lab, angle_theta*pi/180, angle_two_theta*pi/180, dQ, axis_perpendicular)
end
export set_scattering_angles!, set_scattering_angles_deg!


# set transferred momentum dq
function set_dQ!(
            lab             :: LabSystem,
            dQ              :: Real,
            q_beam          :: Real = 10.0,
            axis_perpendicular::Vector{<:Real} = [0,1,0]
        )
    # calculate q's
    q_hor = sqrt((q_beam)^2 - (dQ/2)^2)
    q_in  = [q_hor, 0,  dQ/2]
    q_out = [q_hor, 0, -dQ/2]
    # pass to lab system
    set_q_in!(lab, q_in, axis_perpendicular)
    set_q_out!(lab, q_out, axis_perpendicular)
    # return nothing
    return nothing
end
export set_dQ!