using Molly
using BenchmarkTools

function shake_benchmark()
    n_atoms = 100
    atom_mass = 10.0u"u"
    atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]
    boundary = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")

    coords = place_atoms(n_atoms ÷ 2, boundary, min_dist=0.3u"nm")

    for i in 1:length(coords)
        push!(coords, coords[i].+ [0.15, 0.0, 0.0]u"nm")
    end

    temp = 100.0u"K"

    velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]

    nb_matrix = trues(n_atoms, n_atoms)

    for i in (1:n_atoms÷2)
        nb_matrix[i, i + (n_atoms ÷ 2)] = false
        nb_matrix[i + (n_atoms ÷ 2), i] = false
    end

    neighbor_finder = DistanceNeighborFinder(nb_matrix=nb_matrix, n_steps=10, dist_cutoff=1.5u"nm")
    
    bond_lengths = [0.1u"nm" for i in 1:(n_atoms÷2)]

    sh = SHAKE(bond_lengths, collect(1:(n_atoms÷2)), collect((1+n_atoms÷2):n_atoms))
    
    constraint_list = (sh,)

    constrained_sys = System(atoms=atoms,
                 pairwise_inters=(LennardJones(nl_only=true),),
                            constraints=constraint_list,
                            coords=coords,
                            velocities=velocities,
                            boundary=boundary,
                            neighbor_finder=neighbor_finder,
                            loggers=Dict(
                                         "temp" => TemperatureLogger(10),
                                         "coords" => CoordinateLogger(10),
                                        )
                           )  
    
    

    unconstrained_sys = System(atoms=atoms,
                 pairwise_inters=(LennardJones(nl_only=true),),
                            coords=coords,
                            velocities=velocities,
                            boundary=boundary,
                            neighbor_finder=neighbor_finder,
                            loggers=Dict(
                                         "temp" => TemperatureLogger(10),
                                         "coords" => CoordinateLogger(10),
                                        )
                           )  


    simulator = VelocityVerlet(dt=0.002u"ps", coupling=AndersenThermostat(temp, 1.0u"ps"),)
    
    println("BENCHMARKS FOR CONSTRAINED SIMULATION : ")

    @btime simulate!($constrained_sys, $simulator, 1_000)

    println("BENCHMARKS FOR UNCONSTRAINED SIMULATION : ")

    @btime simulate!($unconstrained_sys, $simulator, 1_000)


end

shake_benchmark()
