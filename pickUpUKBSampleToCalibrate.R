# define input covariate / metabolomics dataset 
UKB_df <-  ######### 
########################################

RequiredList<- read.csv('./SampleListForUKBmatchingToWHII.csv', header = T, fileEncoding = 'UTF-8-BOM')
RequiredList$numberFound <- NA

covariates <- c('bmi0','age0','sex','alcoh_status', 'smoking', 'med')
metVar <- colnames(UKB_df)[!colnames(UKB_df) %in% covariates]

bmi_tol <- 0.5

max_age <- max(RequiredList$age)
min_age <- min(RequiredList$age)

max_BMI <- max(RequiredList$BMI) + bmi_tol
min_BMI <- min(RequiredList$BMI) - bmi_tol 

UKB_ToMatch <- UKB_df[UKB_df$alcoh_status != 2 &
                      UKB_df$smoking != "2" &
                        (is.na(UKB_df$med) | UKB_df$med == 4) &
                        UKB_df$age0 >= min_age & 
                        UKB_df$age0 <= max_age &
                        UKB_df$bmi0 >= min_BMI & 
                        UKB_df$bmi0 <= max_BMI ,]

UKB_ToMatch$pseudoID <- seq(1,nrow(UKB_ToMatch))

match_set <- UKB_ToMatch[,c(covariates, 'pseudoID')]
ID_toPick <- NULL

for (i in 1:nrow(RequiredList)){
  
  #find match in sex
  idx_1 <- which (match_set$sex == RequiredList$sex[i] )
  #find match in age
  idx_2 <- which ( match_set$age0 == RequiredList$age[i])
  #find match in BMI
  idx_3 <- which ( abs(match_set$bmi0 - RequiredList$BMI[i]) < bmi_tol)
  #find which pseudoID is matching requirement
  idx_overlap <- intersect(intersect(idx_1, idx_2), idx_3)
  matchedID <- match_set$pseudoID [idx_overlap]
  
  if (length(matchedID) >= RequiredList$numberNeeded[i] ){
    RequiredList$numberFound[i] <- RequiredList$numberNeeded[i]
    ID_i <- matchedID[1:RequiredList$numberNeeded[i]]
    ID_toPick <- c (ID_toPick, ID_i)
    match_set <- match_set[!match_set$pseudoID %in% ID_toPick, ]
  } else {
    RequiredList$numberFound[i] <- length(matchedID)
    ID_toPick <- c (ID_toPick, matchedID )
    match_set <- match_set[!match_set$pseudoID %in% ID_toPick, ]
  }
}

#print overview

fileConn<-file("./output.txt")
writeLines(paste0('number of samples picked in UKB = ', length(ID_toPick)), 
           fileConn)
close(fileConn)

print(paste0('number of samples picked in UKB = ', length(ID_toPick)))
write.csv(RequiredList,'./CharacterisitcsOfPickedUKB_Samples.csv', row.names = F)

#print metabolite statistics

ImputeZeroWithLowestDivideBy2 <- function(x) {
  x <- as.matrix (x)
  y <- x
  for (i in 1:ncol(x)){
    if(min(x, na.rm = T) <= 0){
      dum <- x[,i]
      dum[dum <= 0] <-NA
      min_x <- min(dum, na.rm =T)
      dum[is.na(dum)] <-min_x /2
      y[,i] <- dum
    }
  }
  y <- data.frame (y)
  colnames (y) <- colnames (x)
  row.names (y) <- row.names (x)
  return(y)
}

metabolite_data <- ImputeZeroWithLowestDivideBy2(UKB_ToMatch[, metVar])
metabolite_data <- metabolite_data[UKB_ToMatch$pseudoID %in% ID_toPick,]

mean_x <- colMeans(metabolite_data , na.rm = T)
mean_logx <- colMeans(log(metabolite_data) , na.rm = T)
NA_x <- apply(metabolite_data, 2, function(x) sum(is.na(x)))

output_statistics <- data.frame('variable' = metVar, 
                                'mean' = mean_x,
                                'meanLog' = mean_logx,
                                'NumNA' = NA_x)

write.csv(output_statistics,'./output_statisticsOfPickedUKB_Samples.csv', row.names = F)

#####
# calculate scaling factors
#####

CharacterisitcsOfPickedUKB_Samples <- read.csv(paste0(file_dir,'/CharacterisitcsOfPickedUKB_Samples.csv') , sep = ';')

#-----

RequiredList <- CharacterisitcsOfPickedUKB_Samples[CharacterisitcsOfPickedUKB_Samples$numberFound > 0, ]

covariate_to_match <- c("age", "sex", "bmi")

cohort1 <- 'WHII'

refset1 <- db [ db$cohort == cohort1 & 
                  db$ethnicity_2 == 1  &
                  (db$hyp_2 == 1 | is.na(db$hyp_2) ) &
                  (db$diab_2 == 1 | is.na(db$diab_2) ) &
                  (db$alc_best_3 != 3 | is.na(db$alc_best_3) ) &
                  (db$smoking_3 != 3 | is.na(db$smoking_3) ) 
                , c("sampleid", covariate_to_match)]

refset1$sex <- as.character(refset1$sex)
refset1$sex [refset1$sex == 1] <- 'Male'
refset1$sex [refset1$sex == 2] <- 'Female'

RequiredList$numberFoundWHII <- NA

bmi_tol <- 0.5
age_tol <- 0.5

match_set <- refset1

ID_toPick <- NULL

for (i in 1:nrow(RequiredList)){
  
  #find match in sex
  idx_1 <- which (match_set$sex == RequiredList$sex[i] )
  #find match in age
  idx_2 <- which ( abs(match_set$age - RequiredList$age[i]) <= age_tol)
  #find match in BMI
  idx_3 <- which ( abs(match_set$bmi - RequiredList$BMI[i]) <= bmi_tol)
  #find which ID is matching requirement
  idx_overlap <- intersect(intersect(idx_1, idx_2), idx_3)
  matchedID <- match_set$sampleid [idx_overlap]
  
  if (length(matchedID) >= RequiredList$numberFound[i] ){
    RequiredList$numberFoundWHII[i] <- RequiredList$numberFound[i]
    ID_i <- matchedID[1:RequiredList$numberFound[i]]
    ID_toPick <- c (ID_toPick, ID_i)
    match_set <- match_set[!match_set$sampleid %in% ID_toPick, ]
  } else {
    RequiredList$numberFoundWHII[i] <- length(matchedID)
    ID_toPick <- c (ID_toPick, matchedID )
    match_set <- match_set[!match_set$sampleid %in% ID_toPick, ]
  }
}
sum(RequiredList$numberFound)
sum(RequiredList$numberFoundWHII)

print(paste0('number of samples picked in WHII = ', length(ID_toPick)))


#print metabolite statistics
ImputeZeroWithLowestDivideBy2 <- function(x) {
  x <- as.matrix (x)
  y <- x
  for (i in 1:ncol(x)){
    if(min(x, na.rm = T) <= 0){
      dum <- x[,i]
      dum[dum <= 0] <-NA
      min_x <- min(dum, na.rm =T)
      dum[is.na(dum)] <-min_x /2
      y[,i] <- dum
    }
  }
  y <- data.frame (y)
  colnames (y) <- colnames (x)
  row.names (y) <- row.names (x)
  return(y)
}
colnames(mets)
new_add_to_calibrate <- c("ACACE", "FAW3_FA", "FAW6_FA",
                          "HDL2_C","HDL3_C" ,"LA" ,
                          "LDL_D","MUFA_FA","PUFA_FA","SFA_FA","UNSAT")
metabolite_data <- mets[db$sampleid %in% ID_toPick, c(var_to_use, new_add_to_calibrate)]
mean_x <- colMeans(metabolite_data , na.rm = T)
mean_logx <- colMeans(log(ImputeZeroWithLowestDivideBy2(metabolite_data)) , na.rm = T)
NA_x <- apply(metabolite_data, 2, function(x) sum(is.na(x)))

output_statistics <- data.frame('variable' = c(var_to_use,new_add_to_calibrate), 
                                'mean' = mean_x,
                                'meanLog' = mean_logx,
                                'NumNA' = NA_x)

write.csv(output_statistics,'./output_statisticsOfPickedWHII_Samples.csv', row.names = F)

###############################################################
UKB_Metstatistics <- read.csv('./output_statisticsOfPickedUKB_Samples.csv') 
WHII_Metstatistics <- read.csv( './output_statisticsOfPickedWHII_Samples.csv') 

which(UKB_Metstatistics$variable %in% WHII_Metstatistics$variable )

nightingale_name_conversion <- read.csv('./nightingale_name_conversion.csv', fileEncoding = 'UTF-8-BOM')
for (i in 1:nrow(nightingale_name_conversion)){
  ind <- which(UKB_Metstatistics$variable %in% nightingale_name_conversion[i,2] )
  UKB_Metstatistics$variable[ind] <- nightingale_name_conversion[i,1]
}

which(UKB_Metstatistics$variable %in% WHII_Metstatistics$variable )
UKB_Metstatistics <- UKB_Metstatistics[match(WHII_Metstatistics$variable, UKB_Metstatistics$variable),]

ratio_mean <- UKB_Metstatistics$mean/ WHII_Metstatistics$mean
scaling_f <-  exp(WHII_Metstatistics$meanLog - UKB_Metstatistics$meanLog) # need to multiple this factor in the UKB dataset


UKB_scaling <- data.frame ('variable' = WHII_Metstatistics$variable, ratio_mean, scaling_f)
write.csv(UKB_scaling,'./UKB_scaling.csv', row.names = F)

