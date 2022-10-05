library("RMVL")

M<-mvl_open("gaia_dr3.mvl")

Msi<-mvl_open("source_id_index.mvl")

# In this example we find indices (i.e. row numbers of Gaia table) corresponding to particular Gaia source ids
#
# Input: a data.frame with source ids. In this example we just list the ids,
#        but for a real application they can be loaded and the data.frame can have other columns

df<-data.frame(source_id=c("3058016772224", 
			"5192616270720", 
			"7284264691456", 
			"8108898970496",
			"not_an_id"), 
		stringsAsFactors=FALSE)


gsid<-M$gaia[,"source_id"]
sid<-df[,"source_id"]

match_index<-function(i, k) {
	F<-gsid[k]==sid[i]
	F[is.na(F)]<-FALSE
	# Gaia source ids are unique
	df[i, "index"]<<-k[F]
	return(NULL)
	}


mvl_index_lapply(Msi$source_id_index, list(df[,"source_id"]), match_index)

# get some Gaia columns

df[,"gaia_sid"]<-M$gaia[df[,"index"],"source_id"]
df[,"gaia_ra"]<-M$gaia[df[,"index"],"ra"]
df[,"gaia_dec"]<-M$gaia[df[,"index"],"dec"]
df[,"gaia_parallax"]<-M$gaia[df[,"index"],"parallax"]

print(df)


