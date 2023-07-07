using StaticArrays

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
