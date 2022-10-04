library("RMVL")

M<-mvl_open("gaia_dr3.mvl")

Mtmp<-mvl_open("gaia_tmp1.mvl", append=TRUE, create=TRUE)

Lhid<-list()
Lhid10<-list()

cat("Computing pixel id\n")
for(i in seq(1, dim(M$gaia)[1], by=1e7)) {
	cat(".")
	idx<-i+(1:1e7)-1
	idx<-idx[idx<=dim(M$gaia)[1]]

	hid<- as.numeric(M$gaia[idx,"source_id"]) %/% as.numeric(34359738368)

	plx<- M$gaia[idx, "parallax"]

	Lhid[[length(Lhid)+1]]<-mvl_write_object(Mtmp,  hid)

	Lhid10[[length(Lhid10)+1]]<-mvl_write_object(Mtmp, ifelse(is.na(plx), 3e60, 
							   ifelse(plx>3.7, hid, 
							   ifelse(plx<0, 2e60, 1e60))))
	}
cat("\n")

mvl_write_object(Mtmp, Lhid, name="Lhid")
mvl_write_object(Mtmp, Lhid10, name="Lhid10")

Mtmp<-mvl_remap(Mtmp)

mvl_fused_write_objects(Mtmp, Mtmp$Lhid, name="hid")
mvl_fused_write_objects(Mtmp, Mtmp$Lhid10, name="hid10")

Mtmp<-mvl_remap(Mtmp)

sidx<-mvl_order_vectors(list(Mtmp$hid10, Mtmp$hid))

Mout<-mvl_open("gaia_dr3.sorted.mvl", append=TRUE, create=TRUE)

mvl_write_object(Mout, "Gaia DR3 data, sorted for easy access", name="description")

mvl_indexed_copy(Mout, M$gaia, sidx, name="gaia")

mvl_close(Mout)

mvl_close(Mtmp)
mvl_close(M)

unlink("gaia_tmp1.mvl")


