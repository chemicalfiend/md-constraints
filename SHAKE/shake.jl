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
    
    box = sys.sys.boundary

    for r in 1:m
        # Atoms that are part of the bond
        i0 = bond_list.is[r]
        i1 = bond_list.js[r]
        
        # Distance vector between the atoms
        r01 = vector(coords[i1], coords[i0], box)
        r2 = norm(r01)^2

        # Distance vector after unconstrained update
        s01 = vector(new_coords[i1], new_coords[i0], box)
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



struct SHAKE3{D, B} <: Constraint
    d1::D
    d2::D
    bond_list::B
    max_iter::Int64
    tol::Float64
end


function SHAKE3(d1, d2, B) return SHAKE3(d1, d2, B, 1e10, 1e-6)

function apply_constraint!(sys, coords, new_coords, dt, constr::SHAKE3)

    bond_list = constr.bond_list
    masses = get_masses(sys)

    m = length(bond_list.is)

    sys.sys.coords = new_coords
    
    box = sys.sys.boundary

    for r in 1:m
        
        A = zeros(2, 2)u"nm^2 / u"

        i0 = bond_list.is[r]
        i1 = bond_list.js[r]
        i2 = bond_list.ks[r]

        m0 = masses[i0]
        m1 = masses[i1]
        m2 = masses[i2]

        r01 = vector(coords[i1], coords[i0], box)
        r01sq = norm(r01)^2
        r02 = vector(coords[i2], coords[i0], box)
        r02sq = norm(r02)^2


        s01 = vector(new_coords[i1], new_coords[i0], box)
        s01sq = norm(s01)^2
        s02 = vector(new_coords[i2], new_coords[i0], box)
        s02sq = norm(s02)^2
        
        A[1, 1] = 2.0 * (1/m0 + 1/m1) * dot(r01, s01)
        A[1, 2] = 2.0 * (1/m0) * dot(s01, r02)
        A[2, 1] = 2.0 * (1/m0) * dot(s02, r01)       
        A[2, 2] = 2.0 * (1/m0 + 1/m1) * dot(r02, s02)

        Ainv = inv(A)

        
        r12 = dot(r01, r02)

        quad1 = zeros(3)u"nm^2 / u^2"
        quad2 = zeros(3)u"nm^2 / u^2"
    
        quad1[1] = (1/m0 + 1/m1)^2 * r01sq
        quad1[2] = (1/m0)^2 * r01sq
        quad1[3] = 2.0 * (1/m0 + 1/m1)*(1/m0)*r12


        quad2[1] = (1/m0 + 1/m2)^2 * r02sq
        quad2[2] = (1/m0)^2 * r02sq
        quad2[3] = 2.0 * (1/m0 + 1/m2)*(1/m0)*r12

        g1 = 0.0
        g2 = 0.0
    
        flag = false
        
        q = zeros(2)u"nm^2"

        for i in 1:constr.maxiter
            q[1] = (quad1[1] * (g1^2)) + (quad1[2] * (g2^2)) + (quad1[3] * (g1 * g2))
   
            q[2] = (quad2[1] * (g1^2)) + (quad2[2] * (g2^2)) + (quad2[3] * (g1 * g2))

            b = zeros(2)u"nm^2"

            b[1] = constr.d1^2 - s01sq - q[1]
            b[2] = constr.d2^2 - s02sq - q[2]

            gnew = Ainv * b
       
    
            if ( (abs(g1 - gnew[1]) < constr.tol) && ((abs(g2 - gnew[2]) < constr.tol)) )
                break
            end


            g1 = gnew[1]
            g2 = gnew[2]

        end
        
        δri0 = r01.*g1 + r02.*g2
        δri1 = r01.*g1
        δri2 = r02.*g2
        
        sys.sys.coords[i0] += δri0
        sys.sys.coords[i1] += δri1
        sys.sys.coords[i2] += δri2

    end

end


struct SHAKE4{D, B} <: Constraint
    d1::D
    d2::D
    d3::D
    bond_list::B
end


function apply_constraint!(sys, coords, new_coords, dt, constr::SHAKE4)
    

    bond_list = constr.bond_list
    masses = get_masses(sys)
    m = length(bond_list.is)
    sys.sys.coords = new_coords
    box = sys.sys.boundary
    
    
    for r in 1:m
        i0 = bond_list.is[r]
        i1 = bond_list.js[r]
        i2 = bond_list.ks[r]
        i3 = bond_list.ls[r]
        
        r01 = vector(coords[i0], coords[i1], box)
        r02 = vector(coords[i0], coords[i1], box)
        r03 = vector(coords[i0], coords[i1], box)

        s01 = vector(new_coords[i0], new_coords[i1], box)
        s02 = vector(new_coords[i0], new_coords[i2], box)
        s03 = vector(new_coords[i0], new_coords[i2], box)



        
    end 




end
