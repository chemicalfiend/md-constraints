# SHAKE implementation ignoring the g^2 term but generalised for the multi-atom case


export 
    SHAKE,
    apply_constraint!,



struct SHAKEALL{D,B,O} <: Constraint
    distances::D
    bond_list::B
    ω::O
    max_iter::Int64
end


function SHAKEALL(d, b)
    return SHAKEALL(d, b, 1.0, 1e10)
end


function apply_constraint!(sys, coords, new_coords, dt , constr::SHAKEALL)

    bonds = constr.bond_list
    
    m = length(bond_list.is)
    masses = get_masses(sys)

    box = sys.boundary
    
    tol = 1e-10

    λ = zeros(m)

    for k in 1:max_iter
        for r in 1:m
            i = bonds.is[r]
            j = bonds.js[r]

            r = vector(coords[i], coords[j], box)
            s = vector(new_coords[i], new_coords[j], box)
            
            rsq = norm(r)^2
            ssq = norm(s)^2
            d2 = constr.distances[r]^2

            diff = d2 - rsq

            mi = masses[i]
            mj = masses[j]

            mij = 2 * (1/mi + 1/mj)
            
            if()    # Convergence criterion
                rdots = dot(r, s)

                if (rdots < d2 * tol)
                    @warn "Malformed input at constraint $r"

                else
                    λij = (1/mij) * constr.ω * diff / rdots
                    λ[r] += λij
                    δr = s * λij
                    
                    sys.sys.coords[i] += δr / mi
                    sys.sys.coords[j] += δr / mj

                end
            end

        end
    end

end
