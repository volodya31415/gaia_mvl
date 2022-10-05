library("RMVL")
library("lattice")

M<-mvl_open("gaia_dr3.mvl")
Mpi<-mvl_open("parallax_index.mvl")

idx<-Mpi$parallax_index[1:1e6]

F<-M$gaia[idx, "astrometric_sigma5d_max"]<0.5 & M$gaia[idx, "grvs_mag"]<20
F[is.na(F)]<-FALSE

idx<-idx[F]

png("plot_example1.png")

print(xyplot(I(10-grvs_mag-5*log10(parallax))~g_rp, M$gaia[idx, c("g_rp", "grvs_mag", "parallax", "teff_gspphot")], 
	group=round(teff_gspphot/1000)*1000, ylab="10-grvs_mag-5*log10(parallax)-10", auto.key=list(columns=5),
	par.settings=list(superpose.symbol = list(pch=19, cex=1, col=rainbow(18))))
	)

dev.off()
