# BlockBootstrap3D.jl

Block bootstrap in 3 dimensions suitable e.g. for gridded meteorological data with both spatial and temporal correlation.

   blockbootstrap_2D_circ(datsize::Tuple{F,F},blocklengths::Tuple{F,F}) where F<:Int   

Function for the moving block bootstrap in two dimensions, returning the bootstrapped indices for your data.
In the circular bootstrap, the array wraps around itself in both dimensions, meaning there is no "cutoff" at 
the border and all samples are equally as likely to be drawn.
On the other hand, blocks drawn at borders are artificial and do not retain their dependency structure

# Examples

```
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


TO DO : Add a function for a 3D block bootstrap, where the trend/periodicity in the temporal dimension is respected. See
https://www.authorea.com/users/95958/articles/361741-confidence-intervals-for-forecast-verification-measures-in-meteorology-using-a-three-dimensional-block-bootstrap

This is loosely based on Seasonal Block Bootstrap by Chan et al (2004), where the seasonal trend in an annual time series was preserved by building temporally sequential blocks of length b << p, where p is the annual period, and the blocks are randomly sampled with replacement from different years. 

Here is an exapmle of a 3D variant, based on reconstructing the annual cycle by using spatially bootstrapped blocks which are sampled and placed sequentially in time in blocks which span integer multiples of the diurnal period. Therefore, the temporally ordered blocks are not sampled from other years, but from other locations. 

For now, here is example code:
```
function myblockbootstrap_time(y,blocklengths_xy,blocklength_time)
   # Take an 3D array y = y(time,longitude,latitude) and do a 3D seasonal block bootstrap,
   # where the temporal non-stationarity is retained
    nt,nx,ny = size(y)
    size_xy = (nx,ny)
    nb_t = floor(Int,nt/blocklength_time)

    # The first temporal block is done outside of loop
    ind_start = 1;
    ind_end = bl_time;
    # Spatial block bootstrap
    inds_xy_bb = blockbootstrap_2D_circ(size_xy,blocklengths_xy);

    y_blockboot = y[ind_start:ind_end,inds_xy_bb];

    # Loop over periods nb_t, building blocks sequentially but bootstrapping in spatial dimensions
    for i = 2:nb_t
        ind_start = 1 + (i-1)*bl_time;
        ind_end = i*bl_time;
        inds_xy_bb = blockbootstrap_2D_circ(size_xy,blocklengths_xy);
        
        y_blockboot = cat(y_blockboot,y_blockboot[ind_start:ind_end,inds_xy_bb],dims=1);
    end

    return y_blockboot
end
```
