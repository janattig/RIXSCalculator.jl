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





















################################################################################
#
#   INDIVIDUAL OPERATORS FOR L AND S
#   PURELY BASED ON THEIR NATURAL QUANTUM NUMBERS
#
################################################################################

# delta function
"""
    delta(i,j)

Standard Kroenecker delta ``\\delta_{i,j}``.
"""
function delta(i, j)
    if i == j
        return 1
    else
        return 0
    end
end

# dagger of something
"""
    dagger(x)

The function computes the complex conjugate of a number, ``x^\\dagger``.
"""
function dagger(x)
    return x'
end


# export delta and dagger
export delta, dagger





# operator functions for L
# corresponding to matrix elements
# < l,mlp | operator | l,ml >
"""
    operatorLz(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L_z \\right| l,m_l \\right>.
"""
function operatorLz(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return delta(mlp, ml)   * ml
end

"""
    operatorLplus(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L^+ \\right| l,m_l \\right>.
"""
function operatorLplus(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return delta(mlp, ml+1) * sqrt(l*(l+1) - ml*(ml+1))
end

"""
    operatorLminus(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L^- \\right| l,m_l \\right>.
"""
function operatorLminus(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return delta(mlp, ml-1) * sqrt(l*(l+1) - ml*(ml-1))
end

"""
    operatorLx(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L_x \\right| l,m_l \\right>.
"""
function operatorLx(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return 0.5 * (operatorLplus(l,mlp,ml) + operatorLminus(l,mlp,ml))
end

"""
    operatorLy(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< l,m_{lp} \\left| L_y \\right| l,m_l \\right>.
"""
function operatorLy(l::Int64, mlp::Int64, ml::Int64) :: Complex{Float64}
    return -0.5*im * (operatorLplus(l,mlp,ml) - operatorLminus(l,mlp,ml))
end



# export operators
export operatorLx, operatorLy, operatorLz, operatorLplus, operatorLminus


# operator functions for S
# corresponding to matrix elements
# < s,msp | operator | s,ms >
"""
    operatorSz(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S_z \\right| s,m_s \\right>.
"""
function operatorSz(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return delta(msp, ms)   * ms
end

"""
    operatorSplus(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S^+ \\right| s,m_s \\right>.
"""
function operatorSplus(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return delta(msp, ms+1) * sqrt(s*(s+1) - ms*(ms+1))
end

"""
    operatorSminus(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S^- \\right| s,m_s \\right>.
"""
function operatorSminus(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return delta(msp, ms-1) * sqrt(s*(s+1) - ms*(ms-1))
end

"""
    operatorSx(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S_x \\right| s,m_s \\right>.
"""
function operatorSx(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return 0.5 * (operatorSplus(s,msp,ms) + operatorSminus(s,msp,ms))
end

"""
    operatorSy(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}

This function corresponds to the matrix element ``\\left< s,m_{sp} \\left| S_y \\right| s,m_s \\right>.
"""
function operatorSy(s::Rational, msp::Rational, ms::Rational) :: Complex{Float64}
    return -0.5*im * (operatorSplus(s,msp,ms) - operatorSminus(s,msp,ms))
end


# export operators
export operatorSx, operatorSy, operatorSz, operatorSplus, operatorSminus


# Spin operators with respect to local coordinate frames
# ASSUMING in coordinate frame cf, the quantisation is given by S = (Sx, Sy, Sz)
# --> new operator in external coordinates

function operatorSx(s::Rational, msp::Rational, ms::Rational, cf::CoordinateFrame) :: Complex{Float64}
    # find out what the x axis is in this frame
    x_prime = get_in_global_coordinates(cf, [1,0,0])
    # return superposition of operators
    return operatorSx(s,msp,ms)*x_prime[1] + operatorSy(s,msp,ms)*x_prime[2] + operatorSz(s,msp,ms)*x_prime[3]
end
function operatorSy(s::Rational, msp::Rational, ms::Rational, cf::CoordinateFrame) :: Complex{Float64}
    # find out what the y axis is in this frame
    y_prime = get_in_global_coordinates(cf, [0,1,0])
    # return superposition of operators
    return operatorSx(s,msp,ms)*y_prime[1] + operatorSy(s,msp,ms)*y_prime[2] + operatorSz(s,msp,ms)*y_prime[3]
end
function operatorSz(s::Rational, msp::Rational, ms::Rational, cf::CoordinateFrame) :: Complex{Float64}
    # find out what the z axis is in this frame
    z_prime = get_in_global_coordinates(cf, [0,0,1])
    # return superposition of operators
    return operatorSx(s,msp,ms)*z_prime[1] + operatorSy(s,msp,ms)*z_prime[2] + operatorSz(s,msp,ms)*z_prime[3]
end
