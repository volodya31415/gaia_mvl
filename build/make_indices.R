library("RMVL")

M<-mvl_open("gaia_dr3.mvl")


Mout<-mvl_open("parallax_index.mvl", append=TRUE, create=TRUE)

mvl_write_object(Mout, "A vector of indices into gaia[,] in decreasing order of parallax", name="description")

mvl_write_object(Mout, mvl_order_vectors(list(M$gaia[,"parallax"]), decreasing=TRUE), name="parallax_index")

mvl_close(Mout)

Mout<-mvl_open("source_id_index.mvl", append=TRUE, create=TRUE)

mvl_write_object(Mout, "A hash-based index on source_id column created with mvl_write_extent_index() ", name="description")

mvl_write_extent_index(Mout, list(M$gaia[,"source_id"]), name="source_id_index")

mvl_close(Mout)

Mout<-mvl_open("coordinates_index.mvl", append=TRUE, create=TRUE)

mvl_write_object(Mout, "A spatial index on source coordinates on 3-dimensional unit sphere", name="description")

dec_rad<-mvl2R(M$gaia[,"dec"])*pi/180
ra_rad<-mvl2R(M$gaia[,"ra"])*pi/180

cd<-cos(dec_rad)

df<-data.frame(x=cos(ra_rad)*cd, y=sin(ra_rad)*cd, z=sin(dec_rad))

mvl_write_object(Mout, df, name="gaia")

Mout<-mvl_remap(Mout)

mvl_write_spatial_index1(Mout, list(Mout$gaia[,"x"], Mout$gaia[,"y"], Mout$gaia[,"z"]), c(12, 12, 12), name="coordinates_index")

mvl_close(Mout)

