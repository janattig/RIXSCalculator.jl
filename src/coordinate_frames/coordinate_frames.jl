"""
    CoordinateFrame

This object defines a coordinate frame located at a `position::Vector{Float64}`, whose orientation is determined by three `Vector{Float64}`: `X`, `Y`, and `Z`.
"""
mutable struct CoordinateFrame

    # the position of the coordinate frame
    position :: Vector{Float64}

    # orientation as coordinate system
    X :: Vector{Float64}
    Y :: Vector{Float64}
    Z :: Vector{Float64}

end

"""
    CoordinateFrame(position::Vector{<:Real})

This function returns an object of type `CoordinateFrame` located at `position`, with standard orientation given by X=[1,0,0],Y=[0,1,0],Z=[0,0,1].

The default position is set to be the origin [0,0,0].
"""
function CoordinateFrame(
        position ::Vector{<:Real} = [0,0,0]
    )
    return CoordinateFrame(position, [1,0,0],[0,1,0],[0,0,1])
end

# export
export CoordinateFrame


function get_in_inner_coordinates(s :: CoordinateFrame, v :: Vector{<:Real})
    return Float64[dot(v, s.X), dot(v, s.Y), dot(v, s.Z)]
end
function get_in_global_coordinates(s :: CoordinateFrame, v :: Vector{<:Real})
    return (s.X .* v[1]) .+ (s.Y .* v[2]) .+ (s.Z .* v[3])
end

function invert(s :: CoordinateFrame) :: CoordinateFrame
    return CoordinateFrame(
        s.position .* -1,
        get_in_inner_coordinates(s, [1,0,0]),
        get_in_inner_coordinates(s, [0,1,0]),
        get_in_inner_coordinates(s, [0,0,1])
    )
end


function set_coordinates!(s :: CoordinateFrame, X::Vector{<:Real}, Y::Vector{<:Real}, Z::Vector{<:Real})
    s.X = X ./ norm(X)
    s.Y = Y ./ norm(X)
    s.Z = Z ./ norm(X)
    return nothing
end








"""
    get_rotation_matrix(Vec :: Vector{<:Real}, Angle :: Real) :: Matrix{Float64}
Given a vector `Vec` found at an angle `Angle`, the function returns the corresponding 3x3 rotation matrix.
"""
function get_rotation_matrix(Vec :: Vector{<:Real}, Angle :: Real) :: Matrix{Float64}
    # normalize the axis
    Axis = Vec ./ norm(Vec)
    R = Matrix{Float64}(undef, 3, 3)
    R[1, 1] = Axis[1]^2 + (1.0 - Axis[1]^2) * cos(Angle)
    R[1, 2] = Axis[1] * Axis[2] * (1.0 - cos(Angle)) - Axis[3] * sin(Angle)
    R[1, 3] = Axis[1] * Axis[3] * (1.0 - cos(Angle)) + Axis[2] * sin(Angle)
    R[2, 1] = Axis[1] * Axis[2] * (1.0 - cos(Angle)) + Axis[3] * sin(Angle)
    R[2, 2] = Axis[2]^2 + (1.0 - Axis[2]^2) * cos(Angle)
    R[2, 3] = Axis[2] * Axis[3] * (1.0 - cos(Angle)) - Axis[1] * sin(Angle)
    R[3, 1] = Axis[1] * Axis[3] * (1.0 - cos(Angle)) - Axis[2] * sin(Angle)
    R[3, 2] = Axis[2] * Axis[3] * (1.0 - cos(Angle)) + Axis[1] * sin(Angle)
    R[3, 3] = Axis[3]^2 + (1.0 - Axis[3]^2) * cos(Angle)
    return R
end
export get_rotation_matrix
