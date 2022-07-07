export
    SHAKE,
    apply_constraint!


struct SHAKE{D, B} <: Constraint
    d::D
    bond_list::B
end


function apply_constraint!(sys, coords, new_coords, dt, constr::SHAKE)
    
    bond_list = constr.bond_list
    masses = get_masses(sys)
    
    
   # print("$masses \n")

    m = length(bond_list.is)
    
    sys.sys.coords = new_coords

    for r in 1:m
        # Atoms that are part of the bond
        i0 = bond_list.is[r]
        i1 = bond_list.js[r]
        
        # Distance vector between the atoms
        r01 = coords[i0] - coords[i1]
        r2 = norm(r01)^2

        # Distance vector after unconstrained update
        s01 = new_coords[i0] - new_coords[i1]
        s2 = norm(s01)^2

        m0 = masses[i0]
        m1 = masses[i1]

        a = (1/m0 + 1/m1)^2 * r2
        b = 2.0 * (1/m0 + 1/m1) * dot(r01, s01)

        c = s2 - constr.d^2

        D = (b^2 - 4*a*c)

        print("a : $a \n  b : $b \n c : $c \n D : $D \n")

        if (ustrip(D) < 0.0)
            @warn "SHAKE determinant negative. Setting to 0.0"
            D = 0.0u"nm^4/u^2"

        end
        
        # Quadratic solution for g

        α1 = (-b + sqrt(D))/(2*a)
        α2 = (-b - sqrt(D))/(2*a)

        g = ((abs(α1) <= abs(α2)) ? α1 : α2)

        # update positions
        
        δri0 = r01.*(g/m0)
        δri1 = r01.*(-g/m1)
        
        sys.sys.coords[i0] += δri0
        sys.sys.coords[i1] += δri1
    end
end



struct SHAKE3 <: Constraint
    
    

end


