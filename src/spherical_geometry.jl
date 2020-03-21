using LinearAlgebra

# NOTE: in all of these functions, the "sky plane" is the x-y plane (unit normal vector aligned with the z-axis), where relevant.

"""
Draw a random unit vector isotropically, in cartesian coordinates.
"""
function draw_random_normal_vector()
    u = 2*rand() - 1 # uniform in [-1,1]
    θ = 2π*rand()

    x = sqrt(1-u^2)*cos(θ)
    y = sqrt(1-u^2)*sin(θ)
    return [x,y,u]
end

Rx(θ) = [1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)] # rotation matrix around x-axis
Ry(θ) = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)] # rotation matrix around y-axis
Rz(θ) = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] # rotation matrix around z-axis

"""
Calculate the rotation matrix that would align vector A to B.
Does not work if A = -B.
"""
function calc_rotation_matrix_A_to_B(A::Vector{Float64}, B::Vector{Float64})
    @assert(length(A) == length(B) == 3)

    v = cross(A, B)
    s = norm(v)
    c = dot(A, B)
    vx = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

    R = I + vx + ((1-c)/s^2)*vx^2
    return R
end

calc_angle_between_vectors(A::Vector{Float64}, B::Vector{Float64}) = acos(dot(A,B)/(norm(A)*norm(B))) # angle between two vectors

"""
Calculate the orientation of a vector in the x-y plane (normal to z-axis).
"""
function calc_Ω_in_sky_plane(h::Vector{Float64})
    n = [-h[2],h[1],0.]
    Ω = acos(n[1]/norm(n))
    n[2] >= 0 ? Ω : 2π - Ω
end

"""
Calculate the angle between two orbits using the spherical law of Cosines.
"""
function calc_incl_spherical_cosine_law(i1::Float64, i2::Float64, ΔΩ::Float64)
    return acos(cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(ΔΩ))
end





"""
    calc_orbit_vector_given_system_vector(i_m, Ω, vec_ref)

Calculate the orbit normal vector given an inclination `i_m` and argument of ascending node `Ω` relative to a system reference plane, in cartesian coordinates.

# Arguments:
- `i_m::Float64`: inclination of orbit (rad) relative to the reference plane.
- `Ω::Float64`: argument of ascending node (rad) relative to the reference plane.
- `vec_ref::Vector{Float64}`: unit normal vector for the reference plane.

# Returns:
- `vec_orb::Vector{Float64}`: unit normal vector for the orbit.
"""
function calc_orbit_vector_given_system_vector(i_m::Float64, Ω::Float64, vec_ref::Vector{Float64})
    @assert(length(vec_ref) == 3)

    # Calculate rotation matrix that rotates z-axis to vec_ref:
    vec_z = [0.,0.,1.]
    inv_R = calc_rotation_matrix_A_to_B(vec_z, vec_ref)

    # Rotate vec_z to get the orbit normal vector:
    vec_orb = Rx(i_m)*vec_z; # rotate by i_m around x-axis
    vec_orb = Rz(Ω)*vec_orb; # rotate by Ω around z-axis
    vec_orb = inv_R*vec_orb; # rotate by inverse rotation matrix

    return vec_orb
end

"""
    calc_sky_incl_Ω_orbits_given_system_vector(i_m_list, Ω_list, vec_ref)

Calculate the inclination and orientation relative to the sky plane for all orbits in the system, given a list of mutual inclinations and arguments of ascending nodes relative to the system reference plane.

# Arguments:
- `i_m_list::Vector{Float64}`: list of orbit inclinations (rad) relative to the reference plane.
- `Ω_list::Vector{Float64}`: list of arguments of ascending node (rad) relative to the reference plane.
- `vec_ref::Vector{Float64}`: unit normal vector for the reference plane.

# Returns:
- `i_sky_list::Vector{Float64}`: list of orbit inclinations (rad) relative to the sky plane.
- `Ω_sky_list::Vector{Float64}`: list of arguments of ascending node (rad) relative to the sky plane.
"""
function calc_sky_incl_Ω_orbits_given_system_vector(i_m_list::Vector{Float64}, Ω_list::Vector{Float64}, vec_ref::Vector{Float64})
    n_pl = length(i_m_list)
    vec_z = [0.,0.,1.]

    i_sky_list = Vector{Float64}(undef, n_pl)
    Ω_sky_list = Vector{Float64}(undef, n_pl)
    for i in 1:n_pl
        vec_orb = calc_orbit_vector_given_system_vector(i_m_list[i], Ω_list[i], vec_ref)
        i_sky_list[i] = calc_angle_between_vectors(vec_z, vec_orb)
        Ω_sky_list[i] = calc_Ω_in_sky_plane(vec_orb)
    end
    return i_sky_list, Ω_sky_list
end
