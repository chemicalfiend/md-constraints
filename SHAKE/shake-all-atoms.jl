#=
    
Iterate through all constraints till all atoms' constraints are simultaneously satisfied.

This follows the LAMMPS quadratic scheme, but uses the iterations from GROMACS to keep all atoms constrained. No longer need to constrain by clusters.

=#


export
    SHAKE,
    apply_constraints!

struct SHAKE{D, B, E}
    dists::D
    is::B
    js::B
    tolerance::E
end

SHAKE(dists, is, js, tolerance=1e-10u"nm") = SHAKE{typeof(dists), typeof(tolerance)}(dists, is, js, tolerance)


function apply_constraints!(sys, constraint::SHAKE, old_coords, dt)
    
    converged = false

    while(!converged)
        for r in 1:length(constraint.is)
            # Atoms that are part of the bond
            i0 = constraint.is[r]
            i1 = constraint.js[r]

            # Distance vector between the atoms before unconstraied update
            r01 = vector(old_coords[i1], old_coords[i0], sys.boundary)

            # Distance vector after unconstrained update
            s01 = vector(sys.coords[i1], sys.coords[i0], sys.boundary)
        
            if abs(norm(s01) - constraint.dists[r]) > constraint.tolerance
                m0 = mass(sys.atoms[i0])
                m1 = mass(sys.atoms[i1])
                a = (1/m0 + 1/m1)^2 * norm(r01)^2
                b = 2 * (1/m0 + 1/m1) * dot(r01, s01)
                c = norm(s01)^2 - ((constraint.dists[r])^2)
                D = (b^2 - 4*a*c)
                
                if (ustrip(D) < 0.0)
                    @warn "SHAKE determinant negative. Setting to 0.0"
                    D = zero(D)
                end

                # Quadratic solution for g
                α1 = (-b + sqrt(D))/(2*a)
                α2 = (-b - sqrt(D))/(2*a)
    
                g = ((abs(α1) <= abs(α2)) ? α1 : α2)
    
                # update positions

                δri0 = r01.*(g/m0)
                δri1 = r01.*(-g/m1)

                sys.coords[i0] += δri0
                sys.coords[i1] += δri1
            end
    
        end

        lengths = [abs(norm(vector(sys.coords[constraint.is[r]], sys.coords[constraint.js[r]], sys.boundary)) - constraint.dists[r]) for r in 1:length(constraint.is)]
        
        if maximum(lengths) < constraint.tolerance
            converged = true
        end
    end        
end
