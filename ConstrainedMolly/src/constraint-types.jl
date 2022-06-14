export 
    Constraint,
    ConstrainedSystem,
    bond_constraint,
    angle_constraint


abstract type Constraint end

function bond_constraint(coord_i, coord_j, d)
    return(norm(coord_i - coord_j)^2 - d^2)
end

function angle_constraint(coord_i, coord_j, coord_k, θ)
    angle = dot((coord_k - coord_j), (coord_j, coord_i)) / (norm(coord_k - coord_j) * norm(coord_j - coord_i))

    return angle - θ # TODO : CHECK THIS !!
end    

mutable struct ConstrainedSystem{D, G, T, A, AD, PI, SI, GI, CN, C, V, B, NF, L, F, E} <: AbstractSystem{D} 
    atoms::A
    atoms_data::D
    pairwise_inters::PI
    specific_inter_lists::SI
    general_inters::GI
    constraints::CN     # Only new addition from System
    coords::C
    velocities::V
    box_size::B
    neighbor_finder::NF
    loggers::L
    force_units::F
    energy_units::E
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
                box_size,
                neighbor_finder=NoNeighborFinder(),
                loggers=Dict(),
                force_units=u"kJ * mol^-1 * nm^-1",
                energy_units=u"kJ * mol^-1",
                gpu_diff_safe=isa(coords, CuArray))
    D = length(box_size)
    G = gpu_diff_safe
    T = typeof(ustrip(first(box_size)))
    A = typeof(atoms)
    AD = typeof(atoms_data)
    PI = typeof(pairwise_inters)
    SI = typeof(specific_inter_lists)
    GI = typeof(general_inters)
    CN = typeof(constraints)
    C = typeof(coords)
    V = typeof(velocities)
    B = typeof(box_size)
    NF = typeof(neighbor_finder)
    L = typeof(loggers)
    F = typeof(force_units)
    E = typeof(energy_units)
    return System{D, G, T, A, AD, PI, SI, GI, CN, C, V, B, NF, L, F, E}(
                    atoms, atoms_data, pairwise_inters, specific_inter_lists,
                    general_inters, constraints, coords, velocities, box_size, neighbor_finder,
                    loggers, force_units, energy_units)
end







