#= 
    
    Testing SHAKE with triatomic molecules

=#

using Molly
using LinearAlgebra

function CO2_shake()
    
    n_atoms = 30
    atom_mass = 10.0u"u"
    atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]
    boundary = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")

    coords = place_atoms(n_atoms÷3, boundary, min_dist=0.3u"nm")

    for i in 1:(n_atoms÷3)
        push!(coords, coords[i].+[0.13, 0.0, 0.0]u"nm")
    end

    for i in 1:(n_atoms÷3)
        push!(coords, coords[i].+[0.26, 0.0, 0.0]u"nm")
    end 
    
    temp = 100.0u"K"

    velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]

    nb_matrix = trues(n_atoms, n_atoms)

    for i in (1:n_atoms÷3)
        nb_matrix[i, i + (n_atoms÷3)] = false
        nb_matrix[i + (n_atoms÷3), i] = false

        nb_matrix[i + (n_atoms÷3), i + (n_atoms*2÷3)] = false
        nb_matrix[i + (n_atoms÷3) + 1, i + (n_atoms*2÷3)] = false
    end

    neighbor_finder = DistanceNeighborFinder(nb_matrix=nb_matrix, n_steps=10, dist_cutoff=1.5u"nm")

    bond_lengths = [0.1u"nm" for i in 1:2*(n_atoms÷3)]
    
    sh = SHAKE(bond_lengths, [collect(1:n_atoms÷3)..., collect((1+(n_atoms÷3)):2*(n_atoms÷3))...], [collect((1+(n_atoms÷3)):(2*n_atoms÷3))..., collect((1+(2*n_atoms÷3)):n_atoms)...])
    rat = RATTLE(bond_lengths, [collect(1:n_atoms÷3)..., collect((1+(n_atoms÷3)):2*(n_atoms÷3))...], [collect((1+(n_atoms÷3)):(2*n_atoms÷3))..., collect((1+(2*n_atoms÷3)):n_atoms)...])
   

    constraints=(sh,rat,)

sys = System(atoms=atoms,
             pairwise_inters=(LennardJones(nl_only=true),),
             constraints=constraints,
             coords=coords,
             velocities=velocities,
             boundary=boundary,
             neighbor_finder=neighbor_finder,
             loggers=Dict(
                          "coords" => CoordinateLogger(10),
                         )
            )

    for i in 1:length(sys.coords)
        sys.coords[i] += [rand()*0.01, rand()*0.01, rand()*0.01]u"nm"
    end

    old_coords = sys.coords

    apply_constraints!(sys, sh, old_coords, 0.002u"ps")
    apply_constraints!(sys, rat, old_coords, 0.002u"ps")
    lengths = []
    vpset = []

    for r in 1:length(sh.is)
        push!(lengths, norm(vector(sys.coords[sh.is[r]], sys.coords[sh.js[r]], sys.boundary)))
        push!(vpset, dot(sys.velocities[rat.is[r]] -sys.velocities[rat.js[r]], vector(sys.coords[rat.is[r]], sys.coords[rat.js[r]], sys.boundary)))    
    end

    print(lengths.-0.1u"nm")
    print(vpset)

end

CO2_shake()

