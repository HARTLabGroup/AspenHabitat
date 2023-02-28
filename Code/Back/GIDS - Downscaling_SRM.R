#######################################################################################
## This code performs a multi-step procedure that spatially downscales gridded climate
  # data using GIDS (Gradient and Inverse Distance-Squared) of Nalder and Wein (1998)and
  # Flint and Flint (2012)

###################################################################################
## Import necessary packages
###################################################################################
package.list <- c("raster", "rgdal", "parallel", "FNN", "maptools", "rgeos",
                  "elevatr")
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###################################################################################
## Creating the functions
###################################################################################
## Quick function for slitting raster to tiles. Useful in broader areas
SplitRas <- function(rast,ntiles_x = 3, ntiles_y = 4){
  h <- ceiling(ncol(rast)/ntiles_x)
  v <- ceiling(nrow(rast)/ntiles_y)
  agg <- aggregate(rast,fact=c(h,v))
  agg[] <- 1:ncell(agg)
  agg_poly <- rasterToPolygons(agg)
  return(agg_poly)
}

## Function to perform GIDS. Called from main function below
gids <- function(coarse_data, fine_data, dist_mat){
  
  ## Reformatting matrices to allow for model fitting and prediction
  coarse_data <- as.data.frame(coarse_data)
  fine_data <- as.data.frame(fine_data)
  worker <- function(i, fine = fine_data, coarse = coarse_data, dists = dist_mat){
    
    ## First, subset data to get local neighbors
    fine_row <- fine[i,]
    ind <- dists[i,,1]
    dist <- dists[i,,2]
    data <- coarse[c(ind),]
    remove(fine, coarse, dists)
    
    ## Fit lm to data using base-level c function. Quite a bit faster than lm()
    x <- model.matrix(~x + y + anc_layer, data = data)
    y <- data$ds_layer
    fit <- .lm.fit(x = x, y = y) # Fit multiple regression model
    coef <- coefficients(fit)[2:4] # Extract coefficients from regression model
    data <- data[5:50,] # Subset data to exclude cells within nugget distance of 1-
      # cell resolution. Remove “goose egg” effect in interpolation
    x_coord <- fine_row$x; y_coord = fine_row$y; anc = fine_row$anc_layer # renaming
      # variables to make the lines (below) shorter
    
    ## Different portions of the GIDS Formula
    d_sq <- dist[5:50]^2
    sum1 <- sum((data$ds_layer + (x_coord - data$x)*coef[1] + (y_coord -
                data$y)*coef[2] + (anc - data$anc_layer)*coef[3])/(d_sq), na.rm = T)
    sum2 <- sum(1/d_sq, na.rm = T) # Doing this step outside the loop would
      # probably increase speed, but also make things harder for me to follow
    ## And returning output
    return(sum1/sum2)
  }
  ## Running worker function on each point and adding as new column to df
  fine_data$ds_layer <- sapply(seq_len(nrow(fine_data)),worker)
  ## Ouptutting df in matrix form for later
  return(as.matrix(fine_data))
}
## Main function. Wrapper for GIDS that formats data and outputs downscaled grid
multiLevelInterpParrallel <- function(boundary, ds_layer, ancillary_list,
                                      clip_dists = c(20000, 5000, 500),file_out = NULL){
  
  ## Define number of cores to be used
  no_cores <- detectCores() - 1
  
  #no_cores <- 3 ## Seems more stable, although slower, to use less than half of available cores
  
  ## Creating projected version of boundary .shp that corresponds with output projection
  boundary_proj <- spTransform(boundary, crs(ancillary_list[[1]]))
  
  ## Looping through each of the DEMs
  for(iter in 1:(length(ancillary_list)-1)){
    if(iter == 1){
      
      ## In first iteration of loop, we need to start with the initial ds_layer and
        # ancillary dataset rather than interpolated data
      ds_layer <- projectRaster(ds_layer,ancillary_list[[1]]) ## Projecting to
        # align with ancillary data grid
      ds_layer <- crop(ds_layer, extend(extent(boundary_proj), clip_dists[iter]))
      
      ## Cropping to study bounds
      anc_layer <- crop(ancillary_list[[1]], extend(extent(boundary_proj),
                                                    clip_dists[iter])) ## Cropping anc data to study bounds
      
      ## Removing cells with NAs in either layer
      ds_layer[is.na(anc_layer)] <- NA
      anc_layer[is.na(ds_layer)] <- NA
      
      ## Converting to non-referenced matrices
      coarse_pts <- rasterToPoints(anc_layer, spatial = F)
      coarse_pts <- cbind(coarse_pts, rasterToPoints(ds_layer, spatial = F)[,3])
      dimnames(coarse_pts)[[2]] <- c("x", "y", "anc_layer", "ds_layer")
      remove(ds_layer, anc_layer)
    
    }else{
      ## More efficient to use previously processed point data in later iterations
      coarse_pts <- out_pts
      remove(out_pts)
    }
    ## Finding nearest neighbors from coarse point layer for each point in finescale ancillary layer
    
    # Cropping raster of finer-scale ancillary data
    new_anc_layer <- crop(ancillary_list[[iter+1]], extend(extent(boundary_proj),
                                                           clip_dists[iter+1]))
    resl <- res(new_anc_layer); refer <- crs(new_anc_layer)
    
    # Converting it to a matrix
    new_pts <- rasterToPoints(new_anc_layer, spatial = F)
    remove(new_anc_layer)
    
    # Defining column names of matrix
    dimnames(new_pts)[[2]] <- c("x", "y", "anc_layer")
    
    # Finding the 50 nearest neighbors and getting distances to each.
      # Based on GIDS typically using 7x7 window
    nns <- get.knnx(coarse_pts[,1:2],new_pts[,1:2], k = 50)
    
    # Reformatting list to 3d array. A little easier to work with later
    nns <- sapply(nns, identity, simplify="array")
    
    ## Splitting data for parrallelization, and initializing cluster
    parts <- split(x = seq_len(nrow(new_pts)), f = 1:no_cores)
    cl <- makeCluster(no_cores)
    clusterExport(cl = cl, varlist = c("coarse_pts", "new_pts", "nns", "parts",
                                       "gids"), envir = environment())
    
    ## Running GIDS, split between number of cores on PC minus 1
    parallelX <- parLapply(cl = cl, X = 1:no_cores,
                           fun = function(t) gids(coarse_data = coarse_pts,
                                                  fine_data = new_pts[parts[[t]],],
                                                  dist_mat = nns[parts[[t]],,]))
    
    ## Terminating cluster
    stopCluster(cl)
    ## Merging parallel output to single matrix
    out_pts <- do.call(rbind, parallelX)
    # Write raster to file if we made it through all iterations and if path for outfile provided
    if(iter == (length(ancillary_list)-1)){
      ## Creating raster from matrix
      downscaled <- rasterFromXYZ(out_pts[,c(1,2,4)], crs = refer,res = resl)
      if(!is.null(file_out)){
        writeRaster(downscaled, file_out, overwrite = T)
      }
    }
    ## Keeping track of progress
    print(paste("Done with downscale", iter, "out of", (length(ancillary_list)-1)))
  }
  ## Return function output to global environment for later use. Comment out if 
    # just writing output to disk
  return(downscaled)
}

###################################################################################
####
## NOTE: The following is a simple example of downscaling with the functions above.
# To use in larger areas and at higher resolutions, data must be split into tiles
# and later merged to deal with memory limits on most computers
###################################################################################
#### Running code in example area
## Boundary of example area
# bounds1<-extent(c(-108.8874,-103.8616,35.28658,42.78432))

bounds1<-extent(c(-106,105,39,40))
# lets making boundry for the SSM areas
# library(sf)
# ssm<-read_sf("C:/Users/Khum/Documents/srm_ecoregion_l3_b100/srm_ecoregion_l3_b100.shp")%>%
#   st_transform(crs=crs(current)) 

## Getting example data from WorldClim and global DEM

current<-raster::stack("C:/Users/paude/Desktop/sdm/climatic_projection/Current_projection.gri")

clim_map<-current$Normal_1991_2020_MAP

## our elevation 
elev<-raster("C:/Users/paude/Desktop/sdm/0_elev_250m/elev250m.tif")

# library(sf)
# library(dplyr)
# srm1<- read_sf("C:/Users/Khum/Downloads/Max_data/Max_data/srm_ecoregion_l3_b100/srm_ecoregion_l3_b100.shp")%>%
#   st_transform(crs(elev1))

elev<-crop(elev,bounds1)

## Projecting both to UTM 13N

clim_map <- projectRaster(clim_map, crs = CRS("+proj=utm +zone=13 +ellps=GRS80 +datum=NAD83
                                      +units=m +no_defs"))
elev <- projectRaster(elev, crs = CRS("+proj=utm +zone=13 +ellps=GRS80 +datum=NAD83
                                      +units=m +no_defs"))

## Creating second DEM at same resolution as climate grid

elev_coarse <- projectRaster(from = elev, to = clim_map)


## NOTE: The goal of the code below is essentially to develop a relationship between
  # the coarse-scale climate data and coarse-scale DEM (aligned and projected to the
  # same grid), and use the higher-resolution DEM and coarse-scale climate data to
  # interpolate higher-resolution climate surfaces. This could also be done with higher-
  # resolution climate data, as we did for annual AET and CWD in this study.
  ## A clipping perimeter is necessary for the function
bound_object <- as(bounds1, "SpatialPolygons")#is this right? Khum dai?
crs(bound_object) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
bound_object <- spTransform(bound_object, CRS("+proj=utm +zone=13 +ellps=GRS80
                                              +datum=NAD83 +units=m +no_defs"))

# Reducing the size of clipping perimeter to reduce edge effects. This is how the function is used
bound_object <- gBuffer(bound_object, width = -20000)

## And actually doing the downscaling, takes a minute or so depending on PC.
  # Warning message happens because number of raster cells is not evenly divisible
  # by number of cores on PC

ds_layer1 <- multiLevelInterpParrallel(boundary = bound_object, ds_layer = clim_map,
                                      ancillary_list = list(elev_coarse, elev),
                                      clip_dists = c(20000, 5000, 500))
ds_layer1


## save the outputs
writeRaster(ds_layer1,"C:/Users/paude/Desktop/sdm/0_elev_250m/map250m1_12b16_01.tif", format="GTiff")

# install.packages("rhr", repos=c("http://R-Forge.R-project.org",
#                                 "http://cran.rstudio.com/"))
