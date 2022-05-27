using Molly
using ConstrainedMolly

#include("../molly-extend/extend.jl")


inter = LennardJones(force_units=NoUnits, energy_units=NoUnits)
box_size = SVector()

constr_list = (ConstraintListMolecules() )

