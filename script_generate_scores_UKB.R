#############################################################################
# load UKB data You need to put your R Object name for your dataframe in here
UKB_df <-     
##############################################

covariates <-
  c('bmi0', 'age0', 'sex', 'alcoh_status', 'smoking', 'med')
metVar <- colnames(UKB_df)[!colnames(UKB_df) %in% covariates]

#############################################
# load R model objects
load('./metaboAgeDB_models.Rdata')
#############################################


##############################
#Akker function (coef in original version - specified in van den akker et.al supplementary section)
##############################
Akker_age <- function (df){
  score <- with (df,
                 58.62		+
                   (	(	(	Acetoacetate	-	0.04319	)	/	0.03465	)	*	1.056	)	+
                   (	(	(	Acetate	-	0.0445	)	/	0.01919	)	*	0.6887	)	+
                   (	(	(	Ala	-	0.3006	)	/	0.0767	)	*	-0.3769	)	+
                   (	(	(	Albumin	-	0.08774)	/	0.006863	)	*	-2.843)	+
                   (	(	(	ApoA1	-	1.594	)	/	0.2002	)	*	0.9177	)	+
                   (	(	(	ApoB	-	0.9732	)	/	0.2213	)	*	8.111	)	+
                   (	(	(	Citrate	-	0.09357	)	/	0.02787	)	*	2.501	)	+
                   (	(	(	Creatinine	-	0.07223	)	/	0.01835	)	*	2.168	)	+
                   (	(	(	DHA	-	0.1447	)	/	0.05478	)	*	-0.8051	)	+
                   (	(	(	Omega_3	-	0.4114	)	/	0.134	)	*	7.364	)	+
                   (	(	(	FAW3_FA	-	3.576	)	/	0.9525	)	*	2.75	)	+
                   (	(	(	Omega_6	-	3.871	)	/	0.7695	)	*	63.88	)	+
                   (	(	(	FAW6_FA	-	33.76	)	/	3.567	)	*	-10.21	)	+
                   (	(	(	Glucose	-	4.8	)	/	1.611	)	*	1.394	)	+
                   (	(	(	Gln	-	0.4528	)	/	0.07966	)	*	3.844	)	+
                   (	(	(	GlycA	-	1.359	)	/	0.2039	)	*	0.1866	)	+
                   (	(	(	HDL2_C	-	0.8999	)	/	0.3046	)	*	-161.3	)	+
                   (	(	(	HDL3_C	-	0.4698	)	/	0.06615	)	*	-35.57	)	+
                   (	(	(	HDL_C	-	1.37	)	/	0.3277	)	*	187.4	)	+
                   (	(	(	HDL_size	-	9.972	)	/	0.2498	)	*	1.254	)	+
                   (	(	(	His	-	0.05925	)	/	0.01482	)	*	-2.084	)	+
                   (	(	(	IDL_C	-	0.6855	)	/	0.1948	)	*	-0.04409	)	+
                   (	(	(	IDL_L	-	1.069	)	/	0.2802	)	*	-3.969	)	+
                   (	(	(	Ile	-	0.05285	)	/	0.02023	)	*	-1.844	)	+
                   (	(	(	L_LDL_L	-	1.171	)	/	0.3441	)	*	-23	)	+
                   (	(	(	LA	-	3.099	)	/	0.6966	)	*	-3.273	)	+
                   (	(	(	Lactate	-	1.232	)	/	1.032	)	*	1.868	)	+
                   (	(	(	LDL_C	-	1.488	)	/	0.5062	)	*	15.22	)	+
                   (	(	(	LDL_size	-	23.65	)	/	0.119	)	*	0.2465	)	+
                   (	(	(	Leu	-	0.06167	)	/	0.01617	)	*	-4.118	)	+
                   (	(	(	M_HDL_L	-	0.8082	)	/	0.1624	)	*	-5.544	)	+
                   (	(	(	M_LDL_L	-	0.6606	)	/	0.2066	)	*	33.28	)	+
                   (	(	(	M_VLDL_L	-	0.6664	)	/	0.3765	)	*	-6.233	)	+
                   (	(	(	MUFA	-	2.904	)	/	0.9051	)	*	-8.85	)	+
                   (	(	(	MUFA_FA	-	24.83	)	/	3.645	)	*	-68.24	)	+
                   (	(	(	Phosphatidylc	-	1.985	)	/	0.3749	)	*	-3.59	)	+
                   (	(	(	Phe	-	0.04339	)	/	0.008491	)	*	2.939	)	+
                   (	(	(	PUFA	-	4.282	)	/	0.8429	)	*	-67.3	)	+
                   (	(	(	PUFA_FA	-	37.34	)	/	3.682	)	*	-61.45	)	+
                   (	(	(	S_HDL_L	-	1.004	)	/	0.1017	)	*	4.017	)	+
                   (	(	(	S_LDL_L	-	0.4282	)	/	0.1244	)	*	-10.42	)	+
                   (	(	(	S_VLDL_L	-	0.7053	)	/	0.2318	)	*	-9.983	)	+
                   (	(	(	Total_C	-	4.42	)	/	0.9903	)	*	-28.22	)	+
                   (	(	(	Total_TG	-	1.397	)	/	0.6807	)	*	4.872	)	+
                   (	(	(	SFA	-	4.375	)	/	0.9667	)	*	-19.06	)	+
                   (	(	(	SFA_FA	-	37.84	)	/	1.867	)	*	-31.47	)	+
                   (	(	(	Sphingomyelins	-	0.4548	)	/	0.08845	)	*	1.705	)	+
                   (	(	(	Cholines	-	2.293	)	/	0.3814	)	*	-3.358	)	+
                   (	(	(	Total_FA	-	11.56	)	/	2.482	)	*	23.67	)	+
                   (	(	(	Phosphoglyc	-	1.898	)	/	0.3669	)	*	5.78	)	+
                   (	(	(	Tyr	-	0.06123	)	/	0.01556	)	*	2.209	)	+
                   (	(	(	Unsaturation	-	1.215	)	/	0.0745	)	*	-0.49	)	+
                   (	(	(	Val	-	0.1574	)	/	0.03815	)	*	1.656	)	+
                   (	(	(	VLDL_C	-	0.8767	)	/	0.2799	)	*	15.42	)	+
                   (	(	(	VLDL_size	-	36.79	)	/	1.361	)	*	3.839	)	+
                   (	(	(	XS_VLDL_L	-	0.5663	)	/	0.1263	)	*	5.976	)
  )
  return(score)
  
}

##############################
#Deelen function
##############################
Deelen_age <- function (met_df , age) {
  scaled_log_df <-
    data.frame(scale(
      log(DistributionAwayFromZero(met_df)),
      center = T,
      scale = T
    ))
  
  score <-  with(
    scaled_log_df,
    XXL_VLDL_L * log (0.80) +
      S_HDL_L * log (0.87) +
      VLDL_size * log (0.85) +
      PUFA_FA * log (0.78) +
      Glucose * log (1.16) +
      Lactate * log (1.06) +
      His * log (0.93) +
      Ile * log (1.23) +
      Leu * log (0.82) +
      Val * log (0.87) +
      Phe * log(1.13) +
      Acetoacetate * log(1.08) +
      Albumin * log(0.89) +
      GlycA * log(1.32)
    
  )
  
  age_sd <- sd (age , na.rm = T)
  age_mean <- mean (age , na.rm = T)
  score_sd <- sd(score, na.rm = T)
  output <-  score / score_sd  * age_sd + age_mean
  
  return(output)
}


calculateMortality_ageScore <- function(met_df, age) {
  output <- rep(NA, nrow(met_df))
  
  
  scale_met_df <- data.frame(scale(met_df))
  
  rel_score <-  with(
    scale_met_df,
    
    S_LDL_TG * 0.01168 + Omega_6 * -0.10244 + Glucose *
      0.00103 + His * -0.07158 + Phe * 0.14083 + Val * -0.01078 + GlycA * 0.07416
    
  )
  
  age_sd <- sd (age , na.rm = T)
  if (age_sd == 0) {
    age_sd <- 0.5
  }
  age_mean <- mean (age , na.rm = T)
  score_sd <- sd(rel_score, na.rm = T)
  output <-  rel_score / score_sd  * age_sd + age_mean
  
  return(output)
}

##########################
# UKB data scaling
##########################
UKB_scaling <- read.csv('./UKB_scaling.csv')

mets <- UKB_df[, metVar]

mets_calibrated <- mets

for (j in 1:ncol(mets_calibrated)) {
  if(colnames(mets_calibrated)[j] %in% UKB_scaling$variable){
    mets_calibrated[, colnames(mets_calibrated)[j]] <-
      mets[, colnames(mets_calibrated)[j]] * UKB_scaling[UKB_scaling$variable == colnames(mets_calibrated)[j], "scaling_f"]
  }
}

####################################
#generate model scores
#####################################

mets_prune <- mets_calibrated [, var_prune_list]

UKB_scores <- UKB_df["age0"]

UKB_scores$EN <-
  predict(model_EN_prune_varObject,
          newx = as.matrix(mets_prune),
          s = model_EN_prune_varObject$lambdaOpt)[,1]
UKB_scores$MARS <-
  predict(model_MARS_prune_var, newdata = as.matrix(mets_prune))[,1]

UKB_scores$EN_pheno <-
  predict (model_EN_Penalty , newx = as.matrix(mets_prune))[,1]

UKB_scores$Akker <- Akker_age(mets_calibrated)

UKB_scores$Deelen <- Deelen_age(mets_calibrated, UKB_df$age0)
UKB_scores$MortalityScore <-
  calculateMortality_ageScore(mets_calibrated, UKB_df$age0)


#########################################################
# Output correlation of the scores
#########################################################
write.csv(cor(UKB_scores), './UKB_score_cor.csv')

