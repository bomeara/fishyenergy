# Fork from https://github.com/troiamj/GB_BEM_v01, making into full package

#' Parameters for Bioenergetics Model (BEM) from FB4
#' 
#' @description This dataset comes from the Fish Bioenergetics 4 (FB4) model, version 1.17, and is from https://github.com/jim-breck/FB4/blob/master/Parameters_official.csv. If you use this, please cite the original source (see fb4_data$Source) as well as Deslauriers et al. (1997), https://doi.org/10.1080/03632415.2017.1377558. It contains common names for aquatic species, scientific names and taxonomy, the life stage modeled, source, and parameters for the BEM. 
#' 
#'  @format A data frame with 105 rows and 44 variables
#' 
#' @source: Deslauriers et al. (1997), https://doi.org/10.1080/03632415.2017.1377558 and https://github.com/jim-breck/FB4/blob/master/Parameters_official.csv, as well as the original source for each row of data.
#' 
"fb4_data"



#' Create a BEM object from FB4 data
#' @param species Scientific name of the species (e.g. "Micropterus salmoides")
#' @param lifestage Life stage of the species (e.g. "Adult"). If NULL, uses the first life stage found for the species.
#' @param CP proportion of maximum consumption (default is 0.5)
#' @param FA fish body mass (grams) (default is 0.1018428)
#' @param UA activity (mg/kg-day) (default is 0.08743042)
#' @param SDA salinity (psu) (default is 0.150819)
#' @param ED energy density (J/kg) (default is 3983.646)
#' @return BEM object with parameters for the specified species and life stage
#' @export
BEM_from_FB4 <- function(species, lifestage=NULL, CP=0.5, FA=0.1018428, UA=0.08743042, SDA=0.150819, ED=3983.646) {
	focal_info <- fb4_data[fb4_data$Sci_Name == gsub("_", " ", species),]
	if(lifestage != NULL) {
  		focal_info <- focal_info[focal_info$Life_Stage == lifestage,]
	}
	if(nrow(focal_info) == 0) {
   		stop(paste("No BEM data found for", species, "with life stage", lifestage))
 	}
	if(nrow(focal_info) > 1) {
	 warning(paste("Multiple BEM data rows found for", species, "with life stage", lifestage, ". Using the first row."))
	 focal_info <- focal_info[1,]
 	}
	bem_object <- data.frame(CP=CP,
		CA=focal_info$CA,
		CB=focal_info$CB,
		CTM=focal_info$CTM,
		CTO=focal_info$CTO,
		CQ=focal_info$CQ,
		ACT=focal_info$ACT,
		RA=focal_info$RA,
		RB=focal_info$RB,
		RTM=focal_info$RTM,
		RTO=focal_info$RTO,
		RQ=focal_info$RQ,
		FA=FA,  
		UA=UA,  
		SDA=SDA, 
		ED=ED
	)
	class(bem_object) <- c('BEM', 'data.frame')
	return(bem_object)
}

#' Get available life stages for a species in the FB4 dataset
#' @param species Scientific name of the species (e.g. "Micropterus salmoides")
#' @return A vector of unique life stages for the species
#' @export
LifeStages_from_FB4 <- function(species) {
 focal_info <- fb4_data[fb4_data$Sci_Name == gsub("_", " ", species),]
 if(nrow(focal_info) == 0) {
	 stop(paste("No BEM data found for", species))
  }
 return(unique(focal_info$Life_Stage))
}


#' Find closest relatives of a species in the FB4 dataset
#' @param species Scientific name of the species (e.g. "Micropterus salmoides")
#' @return A vector of scientific names of the closest relatives in the FB4 dataset
#' @export
#' 
#' @description This uses the Open Tree of Life to find relatives of a species.
Find_Closest_FB4_Species <- function(species) {
  species <- gsub("_", " ", species)
  if(species %in% fb4_data$Sci_Name) {
	return(species) # If the species is already in the dataset, return it directly
  }
  fb4_names_unique <- unique(fb4_data$Sci_Name)
  fb4_names <- gsub(" spp.", "", fb4_names_unique) # Remove " spp." suffix for matching
  fb4_names <- gsub(" X .*$", "", fb4_names, ignore.case=TRUE) # Remove hybrid names
  fb4_names <- unique(fb4_names) # Ensure unique names
  fb4_names <- fb4_names[!fb4_names %in% "Zander zander"] # Remove "Zander zander" as it is not in Open Tree

  resolved_names <- rotl::tnrs_match_names(c(fb4_names, species), context_name = "Animals")
  phy <- suppressWarnings(tol_induced_subtree(ott_ids = ott_id(resolved_names), label_format="name")) # suppress because I really don't care about singleton node suppression
  phy$tip.label <- gsub("_\\(species_in_domain_Eukaryota\\)", "", phy$tip.label)
  phy$tip.label <- gsub("_", " ", phy$tip.label) # Replace underscores with spaces
  parent <- phangorn::Ancestors(phy, species, type="parent")
  children <- phy$tip.label[unlist(phangorn::Descendants(phy, parent, type="tips"))]
  other_relatives <- children[!children %in% species] # Get all other relatives
  final_names <- c()
  for (relative in other_relatives) {
	if (relative %in% fb4_names_unique) {
	  final_names <- c(final_names, relative)
	} else {
		final_names <- c(final_names, fb4_names_unique[which.min( 
	adist(relative, fb4_names_unique))]) #to match spp etc
	}
  }
  return(final_names) 
}

#' Create a new object with BEM data
#' @param CP proportion of maximum consumption
#' @param CA intercept of allometric mass function
#' @param CB slope of allometric mass function; should be negative bc bigger fish consume less per gram of body mass than smaller fish
#' @param CTM critical thermal maximum (degrees C)
#' @param CTO laboratory temperature preferendum (degrees C)=
#' @param CQ approximates a Q10; the rate at which the function increases over relatively low water temperatures
#' @param ACT activity coefficient
#' @param RA intercept of allometric activity function
#' @param RB slope of allometric activity function
#' @param RTM critical thermal maximum (degrees C)
#' @param RTO laboratory temperature preferendum (degrees C)
#' @param RQ approximates a Q10; the rate at which the function increases over relatively low water temperatures
#' @param FA fish body mass (grams)
#' @param UA activity (mg/kg-day)
#' @param SDA salinity (psu)
#' @param ED energy density (J/kg)
#' @return BEM object
#' @export

new_BEM <- function(CP, CA, CB, CTM, CTO, CQ, ACT, RA, RB, RTM, RTO, RQ, FA, UA, SDA, ED) {
	result <- c(CP, CA, CB, CTM, CTO, CQ, ACT, RA, RB, RTM, RTO, RQ, FA, UA, SDA, ED)
	names(result) <- c('CP', 'CA', 'CB', 'CTM', 'CTO', 'CQ', 'ACT', 'RA', 'RB', 'RTM', 'RTO', 'RQ', 'FA', 'UA', 'SDA', 'ED')
	result_df <- data.frame(t(as.matrix(result)))
	class(result_df) <- c('BEM', 'data.frame')
	return(result_df)
}

#' Consumption of energy
#' @param T temperature (degrees C) at which consumption is calculated
#' @param W weight (grams) of fish
#' @param BEM object
#' @return specific consumption rate (grams of food consumed per gram of fish mass per day)
#' @export

consumption <- function(T, W, BEM) {
	Y <- log(BEM$CQ)*(BEM$CTM-BEM$CTO+2)
	Z <- log(BEM$CQ)*(BEM$CTM-BEM$CTO)
	X <- (Z^2*(1+(1+40/Y)^0.5)^2)/400
	V <- (BEM$CTM-T)/(BEM$CTM-BEM$CTO)
	fT_C <- V^X*exp(X*(1-V))
	Cmax <- BEM$CA*W^BEM$CB
	C <- Cmax*BEM$CP*fT_C
	return(C)	
}

#' Respiration 1 equation
#' @param T temperature (degrees C) at which consumption is calculated
#' @param W weight (grams) of fish
#' @param BEM BEM object
respiration <- function(T, W, BEM)
{
  fT <- exp(BEM$RQ*T)
  ACTIVITY <- BEM$ACT
  R <- BEM$RA*W^BEM$RB*fT*ACTIVITY
  return(R)
}

#' Compute growth for a year at a station
#' @param T_vector Temperatures (degrees C), starting at Julian day 1
#' @param BEM BEM object
#' @param starting_weight Starting weight (grams) from Groeschel-Taylor et al
#' @param prey_ED joules per gram of wet mass
#' @param oxycal_coeff joules per gram of oxygen
#' @return data.frame with growth over year
#' @export
compute_single_station <- function(T_vector, BEM, starting_weight=6.382417, prey_ED=3698.0, oxycal_coeff=13560.0) {
	results <- data.frame(julian=sequence(365), temp=T_vector, C1_ins=NA, C2_ins=NA, R1_ins=NA, R2_ins=NA, F1_ins=NA, F2_ins=NA, U1_ins=NA, U2_ins=NA, SDA1_ins=NA, SDA2_ins=NA, W1_ins=NA, W2_ins=NA, W1_cum=NA, W2_cum=NA)
	results[1,] <- 0
	results$W2_cum[1] <- starting_weight
	for (day in 2:nrow(results)) {
		# simulate consumption --> grams of prey
		results$C1_ins[day] <- consumption(T=results$temp[day], W=results$W2_cum[day-1], BEM=BEM) * results$W2_cum[day-1]
		if(is.na(results$C1_ins[day])) {
			results$C1_ins[day] <- 0
		}
		
		# simulate consumption --> joules of energy
		results$C2_ins[day] <- results$C1_ins[day] * prey_ED  # convert from grams of food to joules with prey energy density parameter
		
		# simulate respiration --> grams of oxygen
		results$R1_ins[day] <- respiration(T=results$temp[day], W=results$W2_cum[day-1], BEM=BEM) *  results$W2_cum[day-1]
		
		# simulate respiration --> joules of energy
		results$R2_ins[day] <- results$R1_ins[day] * oxycal_coeff # convert from grams of oxygen to joules with oxycalorific coefficient
		
		# simulate egestion --> grams of prey
		results$F1_ins[day] <- results$C1_ins[day] * BEM$UA # assume egestion is a constant proportion of consumption, for now...
		
		# simulate egestion --> joules of energy
		results$F2_ins[day] <- results$C2_ins[day] * BEM$FA	# assume excretion is a constant proportion of consumption, for now...
		
		# simulate excretion --> grams of prey
		results$U1_ins[day] <- results$C1_ins[day] * BEM$UA # assume excretion is a constant proportion of consumption, for now...  
		
		# simulate excretion --> joules of energy
		results$U2_ins[day] <- results$C2_ins[day] * BEM$UA	# assume excretion is a constant proportion of consumption, for now...
		
		# simulate specific dynamic action --> grams of prey
		results$SDA1_ins[day] <- BEM$SDA * (results$C1_ins[day] - results$F1_ins[day]) # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
		
		# simulate specific dynamic action --> joules of energy
		results$SDA2_ins[day] <- BEM$SDA * (results$C2_ins[day] - results$F2_ins[day]) # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
		
		# simulate daily weight change, that is: C-(M+E+U+SDA) --> joules of energy
		results$W1_ins[day] <- results$C2_ins[day] - (results$R2_ins[day]+ results$F2_ins[day]+ results$U2_ins[day]+ results$SDA2_ins[day])
		
		# simulate daily weight change, that is: C-(M+E+U+SDA) --> grams of body mass
		results$W2_ins[day] <- results$W1_ins[day] / BEM$ED
		
		# convert from joules to grams of body mass with predator energy density parameter
		
		# simulate cumulative weight --> joules of energy
		results$W1_cum[day] <- results$W1_cum[day-1] + results$W1_ins[day]
		
		
		# simulate cumulative weight --> grams of body mass
		results$W2_cum[day] <- results$W2_cum[day-1] + results$W2_ins[day]
		
	}
	results$W1_cum <- cumsum(results$W1_ins)
	return(results)
}