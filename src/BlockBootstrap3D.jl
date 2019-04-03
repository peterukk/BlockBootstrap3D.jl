module BlockBootstrap3D

export blockbootstrap_1D
export blockbootstrap_2D
export blockbootstrap_2D_circ

function blockbootstrap_1D(numobs::Int,blocklength::Int ) where F<:Int
    # 1D moving block bootstrap, returning the bootstrapped indices for the data
    inds = zeros(Int, numobs)
    blockstart_ub = max(1, numobs-blocklength+1) # largest index value where block can start
    for n = 1:blocklength:numobs
        inds[n] = rand(1:blockstart_ub) #Start of block
        for s = n+1:min(n+blocklength-1, numobs) #Iterate through block (use of min avoids bounds error)
            inds[s] = inds[s-1] + 1
        end
    end
    return inds
end

"""
   blockbootstrap_2D_circ(datsize::Tuple{F,F},blocklengths::Tuple{F,F}) where F<:Int   

Function for the moving block bootstrap in two dimensions, returning the bootstrapped indices for your data
In the circular bootstrap, the array wraps around itself in both dimensions, meaning there is no "cutoff" at 
the border and all samples are equally as likely to be drawn
On the other hand, blocks drawn at borders are artificial and do not retain their dependency structure

# Examples

```jldoctest 
julia> using BlockBootstrap3D
julia> datsize = (9,9);
julia> blocklengths_xy = (3,3);
julia> blockbootstrap_2D_circ(datsize,blocklengths_xy)
9×9 Array{CartesianIndex{2},2}:
 CartesianIndex(5, 4)  CartesianIndex(5, 5)  CartesianIndex(5, 6)  CartesianIndex(7, 7)  CartesianIndex(7, 8)  CartesianIndex(7, 9)  CartesianIndex(7, 5)  CartesianIndex(7, 6)  CartesianIndex(7, 7)
 CartesianIndex(6, 4)  CartesianIndex(6, 5)  CartesianIndex(6, 6)  CartesianIndex(8, 7)  CartesianIndex(8, 8)  CartesianIndex(8, 9)  CartesianIndex(8, 5)  CartesianIndex(8, 6)  CartesianIndex(8, 7)
 CartesianIndex(7, 4)  CartesianIndex(7, 5)  CartesianIndex(7, 6)  CartesianIndex(9, 7)  CartesianIndex(9, 8)  CartesianIndex(9, 9)  CartesianIndex(9, 5)  CartesianIndex(9, 6)  CartesianIndex(9, 7)
 CartesianIndex(4, 6)  CartesianIndex(4, 7)  CartesianIndex(4, 8)  CartesianIndex(6, 5)  CartesianIndex(6, 6)  CartesianIndex(6, 7)  CartesianIndex(9, 2)  CartesianIndex(9, 3)  CartesianIndex(9, 4)
 CartesianIndex(5, 6)  CartesianIndex(5, 7)  CartesianIndex(5, 8)  CartesianIndex(7, 5)  CartesianIndex(7, 6)  CartesianIndex(7, 7)  CartesianIndex(1, 2)  CartesianIndex(1, 3)  CartesianIndex(1, 4)
 CartesianIndex(6, 6)  CartesianIndex(6, 7)  CartesianIndex(6, 8)  CartesianIndex(8, 5)  CartesianIndex(8, 6)  CartesianIndex(8, 7)  CartesianIndex(2, 2)  CartesianIndex(2, 3)  CartesianIndex(2, 4)
 CartesianIndex(5, 2)  CartesianIndex(5, 3)  CartesianIndex(5, 4)  CartesianIndex(1, 9)  CartesianIndex(1, 1)  CartesianIndex(1, 2)  CartesianIndex(6, 1)  CartesianIndex(6, 2)  CartesianIndex(6, 3)
 CartesianIndex(6, 2)  CartesianIndex(6, 3)  CartesianIndex(6, 4)  CartesianIndex(2, 9)  CartesianIndex(2, 1)  CartesianIndex(2, 2)  CartesianIndex(7, 1)  CartesianIndex(7, 2)  CartesianIndex(7, 3)
 CartesianIndex(7, 2)  CartesianIndex(7, 3)  CartesianIndex(7, 4)  CartesianIndex(3, 9)  CartesianIndex(3, 1)  CartesianIndex(3, 2)  CartesianIndex(8, 1)  CartesianIndex(8, 2)  CartesianIndex(8, 3)
```

"""
function blockbootstrap_2D_circ(datsize::Tuple{F,F},blocklengths::Tuple{F,F}) where F<:Int

    numobs1,numobs2 = datsize
    bl1,bl2 = blocklengths

    inds_old = CartesianIndices(zeros(datsize))
    inds = zero(inds_old) # New indices to be filled in

    # Make a ghost layer to the right and bottom side of the array, of sizes bl1-1 and bl2-1, respectively, using data 
    # from the first rows and columns 
    # This is to avoid a bounds error and ensure that all elements are equally as likely to be sampled 
    inds_old = [inds_old; view(inds_old,1:bl1-1,:)]
    inds_old = [inds_old view(inds_old,:,1:bl2-1)]

    @inbounds for nx = 1:bl1:numobs1 #starting indices in x dir where the block is "pasted"
        @inbounds for ny = 1:bl2:numobs2 #starting indices in the y dir where block is "pasted"
            i0 = rand(1:numobs1) # starting index where block is copied from
            j0 = rand(1:numobs2) #

            # block lengths, min function is necessary to avoid bounds error near the boundaries where a whole block can't be pasted
            dx = min(bl1-1,numobs1-nx)
            dy = min(bl2-1,numobs2-ny)
            
            inds[nx:nx+dx,ny:ny+dy] = inds_old[i0:i0+dx,j0:j0+dy]
        end
    end
    return inds
end

function blockbootstrap_3D_seasonal(datsize::Tuple{F,F,F},blocklengths::Tuple{F,F,F}) where F<:Int
 # TO DO
end

function blockbootstrap_2D(datsize::Tuple{F,F},blocklengths::Tuple{F,F}) where F<:Int
    # Function for the moving block bootstrap in two dimensions, returning the bootstrapped indices for your data
    # No special treatment of borders, blocks can only be sampled as "whole".
    # this means border data are less likely to be sampled. On the other hand, dependency structure is maintained everywhere
    # and the algorithm is also slightly faster compared to circular boostrap.

    numobs1,numobs2 = datsize
    bl1,bl2 = blocklengths

    inds_old = CartesianIndices(zeros(datsize)) # Old indices, represented as CartesianIndices ([i,j]-matrix)
    inds = zero(inds_old) # New indices to be filled in

    blockstart_ub_x = max(1, numobs1-bl1+1) # upper bound (x-dim) of where a block can start
    blockstart_ub_y = max(1, numobs2-bl2+1) # upper bound (y-dim)
    
    @inbounds for nx = 1:bl1:numobs1 #starting indices in x dir where the block is "pasted"
        @inbounds for ny = 1:bl2:numobs2 #starting indices in the y dir where block is "pasted"
            i0 = rand(1:blockstart_ub_x) # starting index in x dir where block is copied from
            j0 = rand(1:blockstart_ub_y) # starting index in y dir where block is copied from

            # block lengths, min function is necessary to avoid bounds error near the boundaries where a whole block can't be pasted
            dx = min(bl1-1,numobs1-nx)
            dy = min(bl2-1,numobs2-ny)
            
            inds[nx:nx+dx,ny:ny+dy] = inds_old[i0:i0+dx,j0:j0+dy]
        end
    end
    return inds
end

end
