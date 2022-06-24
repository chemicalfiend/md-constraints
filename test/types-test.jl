using ConstrainedMolly
using GLMakie

function system_test()
    
    # ConstrainedSystem tests

    println("Tests for ConstrainedSystem \n\n")
    
    n_atoms = 100
    atom_mass = 10.0u"u"
    atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]
    
    box_size = SVector(2.0, 2.0, 2.0)u"nm"
    coords = place_atoms(n_atoms, box_size, 0.3u"nm")
    
    temp = 100.0u"K"

    velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]
    pairwise_inters = (LennardJones(),)
    
    constraint_list = (NoConstraint(),)
    
    sys = ConstrainedSystem(atoms=atoms,
                            pairwise_inters=pairwise_inters,
                            constraints=constraint_list,
                            coords=coords,
                            velocities=velocities,
                            box_size=box_size,
                            loggers=(coords = CoordinateLogger(10),)
                           )
    return sys
end

function ConstrainedVerletTest()

    sys = system_test()
    temp = 100.0u"K"
    box_size = SVector(2.0, 2.0, 2.0)u"nm"
    simulator = ConstrainedVelocityVerlet(dt=0.01u"ps", coupling=AndersenThermostat(temp, 1.0u"ps"),)
    simulate!(sys, simulator, 1_000)
    visualize(sys.sys.loggers.coords, box_size, "sim_lj.mp4")
end

#system_test()
ConstrainedVerletTest()
