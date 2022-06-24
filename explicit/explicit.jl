# Explicit calculation of matrix inverse for constrained motion of atoms. Refer Ryckaert, Cicottti, Berendsen (1977) Section 2 and 3.


using LinearAlgebra

export
    explicit,
    apply_constraint!,
    constraint


struct explicit{T, B, I, E} <: Constraint
    d::T
    bond_list::B
    maxiter::I
    tol::E
end
 

#function apply_constraint!(bond_list::InteractionList2Atoms, coords, masses, sim, constraint::explicit, cl::CoordinateLogger)

function apply_constraint!(coords, new_coords, dt, constr::explicit)
    
    # not sure if coordinate logger is needed. TODO : Explore better ways to send coords(t) and coords(t + dt) without CL.

    bond_list = constr.bond_list
    #coords = last(cl.coords)
    #new_coords = sys.coords

    σ = []
    m = length(bond_list.is) # Number of constraints = number of bonds
    J = zeros(m, m)  # Initialising the Jacobian
    
    for i in 1:m
        push!(σ, constraint(coords[bond_list.is[i]], coords[bond_list.js[i]], constr.d) )
    end
    
    for l in 1:constr.maxiter
        coord_i = 0
        coord_j = 0
        
        for u in 1:m
            for v in 1:m
                J[u, v] = 2 * norm(coords[bond_list.is[v]] - coords[bond_list.js[v]]) * norm((grad_constraint(coords[bond_list.is[v]], coords[bond_list.js[v]], coords[bond_list.is[u]])/masses[u]) - (grad_constraint(coords[bond_list.js[v]], coords[bond_list.is[v]], coords[bond_list.js[u]])/masses[u])) # CHECK THIS.
            end
        end
        println("\n\n\n\n")

        println(J)

        println(σ)

        λ = - J \ σ      
        
 
        for r in 1:m  # Replace this with enumerate??
            coord_i = new_coords[bond_list.is[r]]
            coord_j = new_coords[bond_list.js[r]]
            
            mi = masses[bond_list.is[r]]
            mj = masses[bond_list.js[r]]

            sumi = zeros(3)
            sumj = zeros(3)
            for k in 1:m
                sumi += (dt^2)* λ[k] * (grad_constraint(coord_i, coord_j, coord_i)/mi) 
                sumj += (dt^2)* λ[k] * (grad_constraint(coord_i, coord_j, coord_j)/mj)
            end
            
            coord_i += sumi
            coord_j += sumj
            
            σnew = []
            
            println("\n\n\n $coord_i \n $coord_j \n\n\n")

            for k in 1:m
                push!(σnew, constraint(coord_i, coord_j, constr.d))
            end

            if (norm(σ - σnew) <= constr.tol)
                coords[bond_list.is[r]] = coord_i
                coords[bond_list.js[r]] = coord_j
                return σnew
            end
            
            σ = σnew

        end


       
        println("\n\n\n")

        println("Iteration $l")
        println(coords)
      
    end

    println("Maximum Iterations Exceeded, σ did not converge")
        
end

function constraint(coord_i, coord_j, d) 
    return (norm(coord_i - coord_j)^2 - (d^2))
end 

function grad_constraint(coord_i1, coord_j1, coord_k)
    
    #coord_i1 and coord_j1 are the coords that the specific constraint is taking care of. coord_k is the coordinate with respect to which we are calculating the gradient.        
    
    if (coord_k == coord_i1)
        return(2*(coord_i1 - coord_j1))
    
    elseif (coord_k == coord_j1)
        return(2*(coord_j1 - coord_i1))

    else
        return 0.0
    
    end

end
