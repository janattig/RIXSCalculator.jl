##################
#  LAB GEOMETRY  #
##################


# set the site positions
"""
    set_site_position!(
        lab      :: LabSystem,
        site     :: Integer,
        position :: Vector{<:Real}
    )

The function sets the position of the site to the input vector `position`.
"""
function set_site_position!(
            lab      :: LabSystem,
            site     :: Integer,
            position :: Vector{<:Real}
        )

    # set the position
    lab.sites[site].position = position
end
export set_site_position!

"""
    rotate_site_in_sample!(
       lab   :: LabSystem,
       site  :: Integer,
       axis  :: Vector{<:Real},
       angle :: Real
    )

The function rotates the site of an `angle` around the given `axis`.
"""
function rotate_site_in_sample!(
            lab   :: LabSystem,
            site  :: Integer,
            axis  :: Vector{<:Real},
            angle :: Real
        )

    # get the rotation matrix
    R = get_rotation_matrix(axis, angle)
    # apply rotation matrix to coordinate system
    set_coordinates!(lab.sites[site], R*lab.sites[site].X, R*lab.sites[site].Y, R*lab.sites[site].Z)
end

"""
    rotate_site_in_sample_deg!(
        lab   :: LabSystem,
        site  :: Integer,
        axis  :: Vector{<:Real},
        angle :: Real
    )

The function rotates the site of an `angle` around the given `axis`. The angle must be given in arc degrees.
"""
function rotate_site_in_sample_deg!(
            lab   :: LabSystem,
            site  :: Integer,
            axis  :: Vector{<:Real},
            angle :: Real
        )

    # pass to the rotation function
    rotate_site_in_sample!(lab, site, axis, angle * pi / 180)
end
export export rotate_site_in_sample!, rotate_site_in_sample_deg!


function rotate_site_z_axis_to_sample_axis!(
            lab   :: LabSystem,
            site  :: Integer,
            axis  :: Vector{<:Real}
        )

    # construct rotation axis
    rotation_axis = cross(axis, lab.sites[site].Z)
    # construct rotation angle
    angle = -acos(dot(axis, lab.sites[site].Z)/(norm(axis)*norm(lab.sites[site].Z)))
    # rotate
    rotate_site_in_sample!(lab, site, rotation_axis, angle)
end
export rotate_site_z_axis_to_sample_axis!

function rotate_site_axis_to_sample_axis!(
            lab   :: LabSystem,
            site  :: Integer,
            site_axis  :: Vector{<:Real},
            sample_axis  :: Vector{<:Real}
        )

    # construct global representation of site axis
    site_axis = get_in_global_coordinates(lab.sites[site], site_axis)
    # construct rotation axis
    rotation_axis = cross(site_axis, sample_axis)
    if norm(rotation_axis) < 1e-15
        if dot(site_axis, sample_axis) > 0
            # print
            println("try to rotate $(site_axis) to $(sample_axis), but already aligned")
            # the axis are already aligned
            return
        else
            # build random rotation axis
            rotation_axis = cross(cross([rand(), rand(), rand()], sample_axis), sample_axis)
            # print
            println("try to rotate $(site_axis) to $(sample_axis), but anti-aligned")
            println("--> using rotation axis $(rotation_axis)")
            # rotate
            rotate_site_in_sample!(lab, site, rotation_axis, pi)
        end
    else
        # construct rotation angle
        angle = acos(dot(sample_axis, site_axis)/(norm(sample_axis)*norm(site_axis)))
        # rotate
        if isfinite(angle)
            rotate_site_in_sample!(lab, site, rotation_axis, angle)
        else
            println("angle not defined")
        end
    end
end
function rotate_sample_axis_to_global_axis!(
            lab   :: LabSystem,
            sample_axis  :: Vector{<:Real},
            global_axis  :: Vector{<:Real}
        )

    # construct global representation of sample axis
    sample_axis = get_in_global_coordinates(lab.sample, sample_axis)
    # construct rotation axis
    rotation_axis = cross(sample_axis, global_axis)
    if norm(rotation_axis) < 1e-15
        if dot(sample_axis, global_axis) > 0
            # print
            println("try to rotate $(sample_axis) to $(global_axis), but already aligned")
            # the axis are already aligned
            return
        else
            # build random rotation axis
            rotation_axis = cross(cross([rand(), rand(), rand()], global_axis), global_axis)
            # print
            println("try to rotate $(sample_axis) to $(global_axis), but anti-aligned")
            println("--> using rotation axis $(rotation_axis)")
            # rotate
            rotate_sample!(lab, global_axis, pi)
        end
    else
        # construct rotation angle
        angle = acos(dot(global_axis, sample_axis)/(norm(global_axis)*norm(sample_axis)))
        # rotate
        if isfinite(angle)
            rotate_sample!(lab, rotation_axis, angle)
        else
            println("angle not defined")
        end
    end
end
export rotate_site_axis_to_sample_axis!, rotate_sample_axis_to_global_axis!

"""
    rotate_sample!(
       lab   :: LabSystem,
       axis  :: Vector{<:Real},
       angle :: Real
   )

The function rotates the sample of an `angle` around `axis`. The angle must be in radians.
"""
function rotate_sample!(
            lab   :: LabSystem,
            axis  :: Vector{<:Real},
            angle :: Real
        )

    # get the rotation matrix
    R = get_rotation_matrix(axis, angle)
    # apply rotation matrix to coordinate system
    set_coordinates!(lab.sample, R*lab.sample.X, R*lab.sample.Y, R*lab.sample.Z)
end
"""
    rotate_sample_deg!(
       lab   :: LabSystem,
       axis  :: Vector{<:Real},
       angle :: Real
   )

The function rotates the sample of an `angle` around `axis`. The angle must be in arc degrees.
"""
function rotate_sample_deg!(
            lab   :: LabSystem,
            axis  :: Vector{<:Real},
            angle :: Real
        )

    # pass to the rotation function
    rotate_sample!(lab, axis, angle * pi / 180)
end
export rotate_sample!, rotate_sample_deg!

"""
    set_site_in_sample_facesharing_towards!(
       ls :: LabSystem,
       site_to_orient :: Integer,
       site_reference :: Integer,
       face_of_reference :: Vector{<:Real} = [1,1,1],
       face_to_orient :: Vector{<:Real} = [1,1,1]
   )

The function sets a site `site_to_orient` so that it faces `site_referece` in a face-sharing configuration.
"""
function set_site_in_sample_facesharing_towards!(
            ls :: LabSystem,
            site_to_orient :: Integer,
            site_reference :: Integer,
            face_of_reference :: Vector{<:Real} = [1,1,1],
            face_to_orient :: Vector{<:Real} = [1,1,1]
        )

    # normalize directions
    face_of_reference = face_of_reference ./ norm(face_of_reference)
    face_to_orient    = face_to_orient ./ norm(face_to_orient)

    # obtain the position within the sample
    position_reference = ls.sites[site_reference].position
    # obtain the direction within the lab frame
    face_direction_reference = get_in_global_coordinates(ls.sites[site_reference], face_of_reference)

    # set the position of the other site within the sample
    set_site_position!(ls, site_to_orient, position_reference .+ face_direction_reference)
    # rotate the site in 5 steps
    # 1. align the coordinate system
    set_coordinates!(ls.sites[site_to_orient], ls.sites[site_reference].X, ls.sites[site_reference].Y, ls.sites[site_reference].Z)
    # 2. rotate z axis around small amount
    rotate_site_in_sample!(ls, site_to_orient, get_in_global_coordinates(ls.sites[site_reference], [0,0,1]), rand()*pi)
    # 3. rotate face to face
    rotate_site_axis_to_sample_axis!(ls, site_to_orient, face_to_orient, face_direction_reference)
    # 4. rotate around face normal by random amount
    rotate_site_in_sample!(ls, site_to_orient, face_direction_reference, rand()*pi)
    # 5. align z projection to face planes respecting the projection of z onto face normal
    if dot([0,0,1], face_to_orient) < 0
        p_orient = cross(cross([0,0,-1], face_to_orient),    face_to_orient)
    else
        p_orient = cross(cross([0,0,+1], face_to_orient),    face_to_orient)
    end
    if dot([0,0,1], face_of_reference) < 0
        p_ref    = cross(cross([0,0,+1], face_of_reference), face_of_reference)
    else
        p_ref    = cross(cross([0,0,-1], face_of_reference), face_of_reference)
    end
    rotate_site_axis_to_sample_axis!(ls, site_to_orient, p_orient, get_in_global_coordinates(ls.sites[site_reference], p_ref))
end
export set_site_in_sample_facesharing_towards!

"""
    center_sites!(ls :: LabSystem)

The function returns the chosen configuration of sites but centered in the mean position.
"""
function center_sites!(
            ls :: LabSystem
        )

    # find mean position
    pos_mean = [
        sum([s.position[i] for s in ls.sites]) for i in 1:3
    ] ./ length(ls.sites)

    # move all sites by mean position
    for s in ls.sites
        s.position .-= pos_mean
    end
end
export center_sites!




# convenience functions
"""
    setup_sites_facesharing!(ls :: LabSystem)

The function returns a face-sharing configuration of the given lab system.
"""
function setup_sites_facesharing!(
            ls :: LabSystem
        )

    # set the first site to be aligned with sample
    ls.sites[1].position = [0,0,0]
    ls.sites[1].X        = [1,0,0]
    ls.sites[1].Y        = [0,1,0]
    ls.sites[1].Z        = [0,0,1]
    # orient towards z axis
    rotate_site_axis_to_sample_axis!(ls, 1, [1,1,0], [0,1,0])
    rotate_site_axis_to_sample_axis!(ls, 1, [1,1,1], [0,0,1])
    rotate_site_in_sample_deg!(ls, 1, [0,0,1], 90)
    # set the rest of the sites to be facesharing chain
    for s in 2:length(ls.sites)
        set_site_in_sample_facesharing_towards!(ls,s,s-1)
    end
    # center the sites
    center_sites!(ls)
    # return nothing
    return nothing
end
export setup_sites_facesharing!