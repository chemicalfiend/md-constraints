using Molly
using LinearAlgebra
constraint_list = (SHAKE([0.1u"nm", 0.1u"nm"], [1, 2], [2, 3]),)

coords = [
    SVector(0.88, 1.0, 1.0)u"nm",
    SVector(1.00, 1.0, 1.0)u"nm",
    SVector(1.12, 1.0, 1.0)u"nm",
]
boundary = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")

sys = System(
    atoms=fill(Atom(), 3),
    constraints=constraint_list,
    coords=coords,
    boundary=boundary,
)

apply_constraints!(sys, sys.coords, 0.002u"ps")

ds = distances(sys.coords, boundary)

println(ds[1, 2]) # 0.1150000000000001 nm
println(ds[2, 3]) # 0.09999999999999987 nm

