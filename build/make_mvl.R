library("RMVL")
library("parallel")

files<-sort(Sys.glob("GaiaSource*gz"))

cat("Processing", length(files), "compressed files\n")

setDefaultCluster(makeForkCluster(32))

gaia2mvl<-function(fn) {
	A<-read.table(gzfile(fn), header=TRUE, sep=",", na.strings=c("NA", "null", "NaN"), stringsAsFactors=FALSE, colClasses=list(solution_id="character", designation="character", source_id="character", random_index="numeric", 
				pmra="numeric", pmdec="numeric", parallax="numeric",
				phot_g_mean_flux="numeric", phot_g_mean_flux_error="numeric",
				phot_bp_mean_flux="numeric", phot_bp_mean_flux_error="numeric",
				phot_rp_mean_flux="numeric", phot_rp_mean_flux_error="numeric",
				phot_bp_n_contaminated_transits="integer", phot_bp_n_blended_transits="integer", phot_rp_n_contaminated_transits="integer", phot_rp_n_blended_transits="integer",
				rv_nb_transits="integer", rv_nb_deblended_transits="integer", rv_visibility_periods_used="integer", rv_atm_param_origin="integer", vbroad_nb_transits="integer",
				grvs_mag_nb_transits="integer"))
# 	F<-A[,"parallax"]>1
# 	F[is.na(F)]<-FALSE
# 	if(sum(F)<1)return(invisible(NULL))
# 	A<-A[F,,drop=FALSE]
	# Fixups:
	gc(reset=TRUE)
	
	for(col in c(grep("_corr$", names(A), value=TRUE), 
			grep("_al$", names(A), value=TRUE),
			grep("scan_direction_", names(A), value=TRUE),
			"ra_error", "dec_error", "parallax_error", "parallax_over_error",
			"pm", "pmra_error", "pmdec_error", 
			"astrometric_excess_noise_sig", 
			"astrometric_sigma5d_max",
			"nu_eff_used_in_astrometry",
			"pseudocolour", "pseudocolour_error", 
			"ipd_gof_harmonic_amplitude", "ipd_gof_harmonic_phase", 
			"ruwe", 
			"phot_g_mean_flux_over_error",
			"phot_g_mean_mag",
			"phot_bp_mean_flux_over_error",
			"phot_bp_mean_mag",
			"phot_rp_mean_flux_over_error",
			"phot_rp_mean_mag",
			"phot_bp_rp_excess_factor",
			"bp_rp", "bp_g", "g_rp",
			"radial_velocity", "radial_velocity_error",
			"rv_expected_sig_to_noise", "rv_renormalised_gof", "rv_chisq_pvalue", "rv_time_duration", "rv_amplitude_robust", "rv_template_teff", "rv_template_logg",
			"rv_template_fe_h", "vbroad", "vbroad_error", 
			"grvs_mag", "grvs_mag_error", 
			"rvs_spec_sig_to_noise", 
			"classprob_dsc_combmod_quasar",
			"classprob_dsc_combmod_galaxy",
			"classprob_dsc_combmod_star",
			"teff_gspphot", "teff_gspphot_lower", "teff_gspphot_upper", "logg_gspphot", "logg_gspphot_lower", "logg_gspphot_upper", 
			"mh_gspphot", "mh_gspphot_lower", "mh_gspphot_upper", "distance_gspphot", "distance_gspphot_lower", "distance_gspphot_upper",
			"azero_gspphot", "azero_gspphot_lower", "azero_gspphot_upper", 
			"ag_gspphot", "ag_gspphot_lower", "ag_gspphot_upper", 
			"ebpminrp_gspphot", "ebpminrp_gspphot_lower", "ebpminrp_gspphot_upper", 
			NULL
			)
			) {
		print(col)
		A[,col]<-as.numeric(A[,col])
		attr(A[,col], "MVL_TYPE")<-4
		gc(reset=TRUE)
		}

        for(col in c("duplicated_source", "astrometric_primary_flag",
			"in_qso_candidates", "in_galaxy_candidates", 
			"has_xp_continuous", "has_xp_sampled", "has_rvs", "has_epoch_photometry", "has_epoch_rv", "has_mcmc_gspphot", "has_mcmc_msc",
			"in_andromeda_survey", 
                        NULL
                        ) 
                        ) {
		print(col)
                A[,col]<-as.logical(A[,col])
		gc(reset=TRUE)
                }

        for(col in c("astrometric_params_solved", 
			"ipd_frac_multi_peak", "ipd_frac_odd_win", "phot_proc_mode", 
			"rv_method_used",
                        NULL
                        )
                        ) {
		print(col)
                A[,col]<-as.integer(A[,col])
                attr(A[,col], "MVL_TYPE")<-1
		gc(reset=TRUE)
                }

	
	bn<-basename(fn)
	M<-mvl_open(paste0("tmp_", bn, ".mvl"), append=TRUE, create=TRUE)
	mvl_write_object(M, A, "data")
	mvl_close(M)
	rm(A)
	gc(reset=TRUE)
	return(invisible(NULL))
	}

clusterApplyLB(cl=NULL, files, gaia2mvl)
stopCluster(cl=NULL)

fn<-Sys.glob("tmp*csv.gz.mvl")
Min<-lapply(fn, mvl_open)

L<-lapply(Min, function(M) { return(M["data", ref=TRUE]) })

Mout<-mvl_open("gaia_dr3.mvl", append=TRUE, create=TRUE)

mvl_fused_write_objects(Mout, L, "gaia")

mvl_close(Mout)

lapply(fn, unlink)

