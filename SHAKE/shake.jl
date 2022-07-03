struct SHAKE{D, B, I, E} <: Constraint
    distance::D
    bond_list::B
    maxiter::I
    tol::E
end


function apply_constraint!(sys, coords, new_coords, dt, constr::SHAKE)
    
    bond_list = constr.bond_list
    masses = get_masses(sys)
    
    m = length(bond_list.is)

    for r in 1:m
        
        # Atoms that are part of the bond

        i0 = bond_list.is[r]
        i1 = bond_list.js[r]
        
        
        # Distance vector between the atoms
        r01 = coords[i0] - coords[i1]
        r = norm(r01)^2

        # Distance vector after unconstrained update
        s01 = new_coords[i0] - new_coords[i1]
        s = norm(s01)^2

        m0 = masses[i0]
        m1 = masses[i1]

        a = (1/m0 + 1/m1)^2 * r
        b = 2.0 * (1/m0 + 1/m1) * dot(r01, s01)

        c = s - constr.d^2

        D = b^2 - 4*a*c

        if (D < 0.0)
            @warn "SHAKE determinant negative. Setting to 0.0"
            D = 0.0

        end
        
        # Quadratic solution for λ

        α1 = (-b + sqrt(D))/(2*a)
        α2 = (-b - sqrt(D))/(2*a)

        λ =(1/dt)^2 * (abs(α1) > abs(α2)) ? α1 : α2

        # update forces
        

end


