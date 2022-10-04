library("RMVL")
#
# In this example we use an input file "points_of_interest.csv" which contains coordinates for several points of interest (ra, dec), specified in degrees
# We use previously constructed spatial index to find points in the Gaia data that are closest to points of interest
# Results are printed in the end of the file
#

M<-mvl_open("gaia_dr3.mvl")

Mci<-mvl_open("coordinates_index.mvl")

x<-Mci$gaia[,"x"]
y<-Mci$gaia[,"y"]
z<-Mci$gaia[,"z"]

# User provided file. Example included.
# For very large files, it might be better to store data in MVL format
points_of_interest<-read.table("points_of_interest.csv", header=TRUE)

x_poi<-cos(points_of_interest[,"ra"]*pi/180.0)*cos(points_of_interest[,"dec"]*pi/180.0)
y_poi<-sin(points_of_interest[,"ra"]*pi/180.0)*cos(points_of_interest[,"dec"]*pi/180.0)
z_poi<-sin(points_of_interest[,"dec"]*pi/180.0)



f_closest<-function(i, idx) {
	# Find closest (on the sky) point to our point of interest
	d<-(x[idx]-x_poi[i])^2+(y[idx]-y_poi[i])^2+(z[idx]-z_poi[i])^2
	return(idx[which.min(d)])
	}

closest_indices<-unlist(mvl_index_lapply(Mci$coordinates_index, list(x_poi, y_poi, z_poi), f_closest))

Ffound<-!is.na(closest_indices)
i0<-1:length(x_poi)
i_poi<-i0[Ffound]
i_gaia<-closest_indices[Ffound]

df<-data.frame(points_of_interest[i_poi,,drop=FALSE])
names(df)<-paste0(names(df), "_poi")
df<-data.frame(df, mvl2R(M$gaia[i_gaia, ]))

print(df)

# You can write data to a text file with the following command:
# write.table(df, "spatial_index_example1_output.csv", col.names=TRUE, row.names=FALSE, sep="\t")
#




