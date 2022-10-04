library("RMVL")

#
# We find all points with proper motion above a given cutoff by full scan of Gaia data
#

MinPM<-1000

M<-mvl_open("gaia_dr3.mvl")

Mpm<-M$gaia[,"pm"]

N<-dim(M$gaia)[1]

Nchunk<- 16e6

L<-list()

for(i in seq(1, N, Nchunk)) {
	cat(".")
	idx<-i-1+(1:Nchunk)
	idx<-idx[idx<=N]
	
	Fpm<-Mpm[idx]>MinPM
	Fpm[is.na(Fpm)]<-FALSE
	
	if(any(Fpm)) {
		L[[length(L)+1]]<-idx[Fpm]
		}
	}
cat("\n")

pm_idx<-do.call(c, L)

cat("Found", length(pm_idx), "points with proper motion above", MinPM, " mas/yr\n")

cat("Creating table with results\n")

# Write out results in MVL and CSV formats

Mout<-mvl_open("full_scan1_high_pm.mvl", create=TRUE, append=TRUE)

mvl_write_object(Mout, paste("Gaia sources with proper motion above", MinPM, " mas/yr"), name="description")

mvl_write_object(Mout, pm_idx, name="high_pm_idx")

mvl_indexed_copy(Mout, M["gaia", ref=TRUE], pm_idx, name="high_pm")

Mout<-mvl_remap(Mout)

high_pm<-mvl2R(Mout$high_pm)

write.table(high_pm, "full_scan1_high_pm.csv", col.names=TRUE, row.names=FALSE, sep="\t")

mvl_close(Mout)






