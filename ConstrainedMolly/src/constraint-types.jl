export 
    get_masses,
    Constraint,
    NoConstraint,
    apply_constraint!,
    ConstrainedSystem,
    bond_constraint,
    angle_constraint,
    ConstrainedVelocityVerlet,
    simulate!,
    BondLogger,
    log_property!


abstract type Constraint end


function get_masses(sys)
    
    masses = []
    for atom in sys.sys.atoms
        push!(masses, atom.mass)
    end
    return masses

end

struct NoConstraint <: Constraint end

function apply_constraint!(sys, coords, new_coords, sim, constr::NoConstraint)
    
    sys.sys.coords = new_coords
end


function bond_constraint(coord_i, coord_j, d)
    return(norm(coord_i - coord_j)^2 - d^2)
end

function angle_constraint(coord_i, coord_j, coord_k, θ)
    angle = dot((coord_k - coord_j), (coord_j, coord_i)) / (norm(coord_k - coord_j) * norm(coord_j - coord_i))

    return angle - θ # TODO : CHECK THIS !!
end    

mutable struct ConstrainedSystem{S, CN}
    sys::S
    constraints::CN     # Only new addition from System
end

function ConstrainedSystem(;
                atoms,
                atoms_data=[],
                pairwise_inters=(),
                specific_inter_lists=(),
                general_inters=(),
                constraints=(),
                coords,
                velocities=zero(coords) * u"ps^-1",
                boundary,
                neighbor_finder=NoNeighborFinder(),
                loggers=(),
                force_units=u"kJ * mol^-1 * nm^-1",
                energy_units=u"kJ * mol^-1",
                k=Unitful.k,
                gpu_diff_safe=isa(coords, CuArray))
    D = n_dimensions(boundary)
    G = gpu_diff_safe
    T = float_type(boundary)
    A = typeof(atoms)
    AD = typeof(atoms_data)
    PI = typeof(pairwise_inters)
    SI = typeof(specific_inter_lists)
    GI = typeof(general_inters)
    CN = typeof(constraints)
    C = typeof(coords)
    V = typeof(velocities)
    B = typeof(boundary)
    NF = typeof(neighbor_finder)
    L = typeof(loggers)
    F = typeof(force_units)
    E = typeof(energy_units)
    K = typeof(k)
    
    s = System{D, G, T, A, AD, PI, SI, GI, C, V, B, NF, L, F, E, K}(atoms, atoms_data, pairwise_inters, specific_inter_lists, general_inters, coords, velocities, boundary, neighbor_finder, loggers, force_units, energy_units, k)
    
    S = typeof(s)

    return ConstrainedSystem{S, CN}(s, constraints)
end


function run_constraints!(sys, coords, new_coords, dt)
    for cons in sys.constraints
        apply_constraint!(sys, coords, new_coords, dt, cons)
    end

end


#function remove_molar(x)
#    fx = first(x)
#    return ([ustrip(x[1]), ustrip(x[2]), ustrip(x[3])]u"nm *  ps^-2")

#end

function remove_molar(x)
    fx = first(x)
    if dimension(fx) == u"𝐋 * 𝐍^-1 * 𝐓^-2"
        T = typeof(ustrip(fx))
        return x / T(Unitful.Na)
    else
        return x
    end
end


struct BondLogger{T}
    n_steps::Int
    is::Vector{Int64}
    js::Vector{Int64}
    bond_lengths::Vector{T}
end

function BondLogger(n_steps::Int64, is, js)
    return BondLogger(n_steps, is, js, [])

end


function log_property!(logger::BondLogger, s::System,neighbors=nothing, step_n::Integer=0; parallel::Bool=true)
    if step_n % logger.n_steps == 0
        lengths = []
        for k in 1:length(logger.is)
            l = norm(s.coords[logger.is[k]] - s.coords[logger.js[k]])
            push!(lengths, l)
        end
        push!(logger.bond_lengths, lengths)
    end
end



mutable struct ConstrainedVelocityVerlet{T, C} 
    dt::T
    coupling::C
    remove_CM_motion::Bool
end

function ConstrainedVelocityVerlet(; dt, coupling=NoCoupling(), remove_CM_motion=true)
    return ConstrainedVelocityVerlet(dt, coupling, remove_CM_motion)
end

function simulate!(sys::ConstrainedSystem,
                    sim::ConstrainedVelocityVerlet,
                    n_steps::Integer;
                    parallel::Bool=true)
    neighbors = find_neighbors(sys.sys, sys.sys.neighbor_finder; parallel=parallel)
    run_loggers!(sys.sys, neighbors, 0; parallel=parallel)
    accels_t = accelerations(sys.sys, neighbors; parallel=parallel)
    accels_t_dt = zero(accels_t)
    sim.remove_CM_motion && remove_CM_motion!(sys.sys)

    for step_n in 1:n_steps
        new_coords = sys.sys.coords + sys.sys.velocities .* sim.dt .+ (remove_molar.(accels_t) .* sim.dt ^ 2) ./ 2

        new_coords = wrap_coords.(new_coords, (sys.sys.boundary,))
        
        #sys.sys.coords = new_coords

        accels_t_dt = accelerations(sys.sys, neighbors; parallel=parallel)

        sys.sys.velocities += remove_molar.(accels_t .+ accels_t_dt) .* sim.dt / 2

        sim.remove_CM_motion && remove_CM_motion!(sys.sys)

        apply_coupling!(sys.sys, sim, sim.coupling)

        run_constraints!(sys, sys.sys.coords, new_coords, sim.dt)

        run_loggers!(sys.sys, neighbors, step_n;parallel=parallel)

        if step_n != n_steps
            neighbors = find_neighbors(sys.sys, sys.sys.neighbor_finder, neighbors, step_n; parallel=parallel)
            accels_t = accels_t_dt
        end
    end
    return sys.sys
end


