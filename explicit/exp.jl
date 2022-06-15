# Explicit matrix inverse single iteration

struct bL
    is,
    js
end


function one_step():

    coords = [ [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [6.0, 0.0, 0.0], [7.0, 0.0, 0.0] ]
    
    bond_list = bL([1, 3], [2, 4])
    
    masses = [1.0, 1.0]

    d = 1.0

    λ = [0, 0]

    σ = [0, 0]
    
    J = []


end
