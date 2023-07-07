L = 10
R = 4
r = 0.06
nbWires = 24
nbNodesCell = 4
nbCells = 10
nbInters = nbCells*2+1
ϕ = L/nbCells

pos = Vector{Vec3{Float64}}()
nodeID = Vector{Int}()
push!(nodeID, 1)
conn = Vector{Vec2{Int}}()

for n in 1:nbWires
    θ = (n-1)*(2*pi/nbWires)
    for i in 1:nbInters
        for j in 1:nbNodesCell+1
            a(j) = !iseven(i) * (2*r/nbNodesCell*(j-1)) + iseven(i) * (2*r - 2*r/nbNodesCell*(j-1))
            A = R  + a(j)
            Θ = θ + pi/nbWires*(i-1) + pi/(nbWires*nbNodesCell)*(j-1)

            x = A * cos(Θ)
            y = A * sin(Θ)
            z = R*pi*(i-1)/(nbWires*tan(ϕ)) + R*pi*(j-1)/(nbWires*nbNodesCell*tan(ϕ))

            push!(pos, [x, y, z])

            if (i-1)*nbNodesCell + j  > 1
                push!(conn, [nodeID[end-1], nodeID[end]]) 
            end 

            push!(nodeID, nodeID[end]+1)
        end
    end 
end 

for n in 1:nbWires
    θ = (n-1)*(2*pi/nbWires)
    for i in 1:nbInters
        for j in 1:nbNodesCell+1
            a(j) = iseven(i) * (2*r/nbNodesCell*(j-1)) + !iseven(i) * (2*r - 2*r/nbNodesCell*(j-1))
            A = R + a(j)
            Θ = θ - pi/nbWires*(i-1) - pi/(nbWires*nbNodesCell)*(j-1)

            x = A * cos(Θ)
            y = A * sin(Θ)
            z = R*pi*(i-1)/(nbWires*tan(ϕ)) + R*pi*(j-1)/(nbWires*nbNodesCell*tan(ϕ))

            push!(pos, [x, y, z])

            if (i-1)*nbNodesCell + j > 1
                push!(conn, [nodeID[end-1], nodeID[end]]) 
            end 

            push!(nodeID, nodeID[end]+1)
        end
    end 
end 

pos_mat = reshape(reinterpret(Float64, pos), (3, length(pos)))'
max_z = maximum(pos_mat[:,3])
min_z = minimum(pos_mat[:,3])
pos_mat[:,3] = L.*(pos_mat[:,3].-min_z)./(max_z - min_z)
pos = [Vec3(pos_mat[i,1],pos_mat[i,2],pos_mat[i,3]) for i in 1:size(pos_mat,1)]

write_vtk_configuration("stent.vtk", pos, conn)

