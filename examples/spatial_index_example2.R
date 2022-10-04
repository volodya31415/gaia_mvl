library("RMVL")
#
# In this example we find stars close to a specific point on the sky (ra, dec), specified in degrees
# We use previously constructed spatial index to find points in the Gaia data that are within a spherical distance d, specified in degrees
# Results are printed in the end of the file
#

M<-mvl_open("gaia_dr3.mvl")

Mci<-mvl_open("coordinates_index.mvl")

x<-Mci$gaia[,"x"]
y<-Mci$gaia[,"y"]
z<-Mci$gaia[,"z"]

# point on the sky, degrees
ra <- 45
dec <- 60

# spherical distance, degrees
d <- 0.1



x_poi<-cos(ra*pi/180.0)*cos(dec*pi/180.0)
y_poi<-sin(ra*pi/180.0)*cos(dec*pi/180.0)
z_poi<-sin(dec*pi/180.0)

cos_d<-cos(d*pi/180.0)

# The spatial index uses 12 bits for each x, y and z coordinate, with resolution of 2^-11 rad
# Construct a grid of query points with this resolution:
N<- ceiling(d*pi/180/2^-11)

steps<- ((-N):N)*2^-11
x_grid<- rep(steps, 2*N+1)
y_grid<- rep(steps, each=2*N+1)
z_grid<- 1.0

norm<-1.0/sqrt(x_grid^2+y_grid^2+z_grid^2)

x_grid<- x_grid*norm
y_grid<- y_grid*norm
z_grid<- z_grid*norm

# Rotate the grid so that (0, 0, 1) point becomes (x_poi, y_poi, z_poi)

x_grid0<- x_grid*sin(dec*pi/180.0)+z_grid*cos(dec*pi/180.0)
z_grid0<- -x_grid*cos(dec*pi/180.0)+z_grid*sin(dec*pi/180.0)

x_grid<-x_grid0
z_grid<-z_grid0

x_grid0<- x_grid*cos(ra*pi/180.0)-y_grid*sin(ra*pi/180.0)
y_grid0<- x_grid*sin(ra*pi/180.0)+y_grid*cos(ra*pi/180.0)

x_grid<-x_grid0
y_grid<-y_grid0


f_nearby<-function(i, idx) {
	# Find nearby (on the sky) points to our point of interest
	d<-x[idx]*x_poi+y[idx]*y_poi+z[idx]*z_poi
	return(idx[d>= cos_d])
	}

nearby_indices<-sort(unique(do.call(c, mvl_index_lapply(Mci$coordinates_index, list(x_grid, y_grid, z_grid), f_nearby))))

# Write out results in MVL and CSV formats

Mout<-mvl_open("spatial_index_example2_output.mvl", create=TRUE, append=TRUE)

mvl_write_object(Mout, paste("Gaia sources near ra=", ra, "dec=", dec), name="description")

mvl_write_object(Mout, nearby_indices, name="point_index")

mvl_indexed_copy(Mout, M["gaia", ref=TRUE], nearby_indices, name="gaia_subset")

Mout<-mvl_remap(Mout)

df<-mvl2R(Mout$gaia_subset)

write.table(df, "spatial_index_example2_output.csv", col.names=TRUE, row.names=FALSE, sep="\t")

mvl_close(Mout)

# plot(df[,"ra"], df[,"dec"])

