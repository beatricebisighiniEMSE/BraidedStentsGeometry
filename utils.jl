using StaticArrays
using DelimitedFiles
using LinearAlgebra

const Vec2{T} = SVector{2,T}
const Vec3{T} = SVector{3,T}

function write_vtk_configuration(filename, positions, connectivity)
    
    fid = open(filename, "w")
    
    write(fid,"# vtk DataFile Version 3.0")
    write(fid,"\nvtk output")
    write(fid,"\nASCII")
    write(fid,"\nDATASET POLYDATA")
    
    numberPoints = length(positions)
    write(fid,"\nPOINTS $numberPoints float")  
    for n in positions
        write(fid,"\n")
        write(fid, string.(n[1]))
        write(fid," ")
        write(fid, string.(n[2]))
        write(fid," ")
        write(fid, string.(n[3])) 
    end 
    
    numberLines = length(connectivity)
    numberElementsPerLine = 2
    numberElements= numberLines*(numberElementsPerLine + 1)
    write(fid,"\nLINES $numberLines $numberElements\n")
    
    for i in 1:numberLines
        
        write(fid, string.(numberElementsPerLine))
        write(fid,"\n")
        
        node = connectivity[i][1] -1
        write(fid, string.(node)) 
        write(fid,"\n")
        
        node = connectivity[i][2] -1
        write(fid, string.(node)) 
        write(fid,"\n")
        
    end 
    
    close(fid)
    
end 

function get_rings(positions)
    
    zvec = Vector{Float64}()
    for p in positions 
        push!(zvec, p[3])
    end 
    # zvec .= round.(zvec, digits = 5)
    unique!(zvec)

    rings = Vector{Vector{Int}}()
    nrings = length(zvec)
    for i in 1:nrings
        ring = Vector{Int}()
        for (k,p) in enumerate(positions)
            if isapprox(p[3], zvec[i])
                push!(ring, k)
            end 
        end 
        push!(rings, ring)
    end
    
    return rings 
    
end 

function get_nodespairs_stent(positions)
    
    constraints = Vector{Vec2{Int}}()
    toll = 0.01
    rings = get_rings(positions)

    for r in 1:length(rings)
        
        for i in 1:length(rings[r])
            
            index_i = rings[r][i]
            pos_i = positions[index_i]
            
            for j in 1:length(rings[r])
                
                index_j = rings[r][j]
                pos_j = positions[index_j]
                
                dist = norm(pos_i-pos_j)
                
                if index_j != index_i && dist < toll && !([index_j, index_i] in constraints)
                    push!(constraints, [index_i, index_j])
                end 
                
            end
        end 
    end 
    
    return constraints
    
end 