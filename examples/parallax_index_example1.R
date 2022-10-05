library("RMVL")

M<-mvl_open("gaia_dr3.mvl")

Mpi<-mvl_open("parallax_index.mvl")


# Find indices of stars in Gaia within a certain distance

# distance is in parsecs
nearby_stars<-function(distance, parallax=1000/distance) {
	k<-1024
	while(TRUE) {
		x<-M$gaia[Mpi$parallax_index[k], "parallax"]
		if(is.na(x) || x<parallax)break
		k<-k*2
		if(k>length(Mpi$parallax_index))break
		}
	idx<-Mpi$parallax_index[1:k]
	F<-M$gaia[idx, "parallax"]>=parallax
	F[is.na(F)]<-FALSE
	return(idx[F])
	}

# Find nearby stars within 100 pc
nearby_idx<-nearby_stars(100)
	
cat("Found", length(nearby_idx), " stars within 100 pc\n")

cat("Creating table with results\n")

# Write out results in MVL and CSV formats

Mout<-mvl_open("parallax_index_example1.mvl", create=TRUE, append=TRUE)

mvl_write_object(Mout, paste("Gaia sources within 100 pc"), name="description")

mvl_write_object(Mout, nearby_idx, name="nearby_idx")

mvl_indexed_copy(Mout, M["gaia", ref=TRUE], nearby_idx, name="gaia_100pc")

Mout<-mvl_remap(Mout)

high_pm<-mvl2R(Mout$high_pm)

write.table(high_pm, "parallax_index_example1.csv", col.names=TRUE, row.names=FALSE, sep="\t")

mvl_close(Mout)
