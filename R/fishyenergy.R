### TODO

# Fish Bioenergetics PDF from 1997 from version 3 has appendices with info on prey energy. 
# Check FB4 for egestion, etc. parameters.
# Some params can be vectors
# New work has condition factor: weight to body mass. FishBase has length/weight relationships, can use that to convert weight to length up to a point. Fish can have a minimum viable body condition (length can go up, weight can go up or down). Could be for an individual based model
# Check for fish tree of life: Jonathan Chang. Fish Tree of Life. Fish Phylomaker can insert species that aren't in the phylogeny -- could be better than OToL. https://cran.r-project.org/web/packages/fishtree/fishtree.pdf
# Punzet 2012 -- for converting air temp to water temperature for different climate zones. Could add that to this package. Punzet, M., Voß, F., Voß, A., Kynast, E., & Bärlund, I. (2012). A global approach to assess the potential impact of climate change on stream water temperatures and related in-stream first-order decay rates. Journal of Hydrometeorology, 13(3), 1052-1065.
# Put in error checking for the temperature, also handling missing data, especially stuff like sensor removed when there's going to be ice. Maybe figure out way to determine these sorts of error, maybe compare with air. 
# More consumption functions possible, could put in. 
# https://cran.r-project.org/web/packages/FishPhyloMaker/readme/README.html
# Add matt as collaborator. 

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
#' @param lifestage Life stage of the species (e.g. "Adult"). If NA, uses the first life stage found for the species.
#' @param CP proportion of maximum consumption (default is 0.5)
#' @param FA egestion proportion (what proportion of food gets passed through) (default is 0.1018428)
#' @param UA excretion proportion (what proportion is assimilated but not allocated towards metabolism or gonad growth) (default is 0.08743042)
#' @param SDA specific dynamic action, proportion of energy going to digestion (default is 0.150819)
#' @param ED energy density (J/kg) (default is 3983.646)
#' @return BEM object with parameters for the specified species and life stage
#' @export
bem_from_fb4 <- function(species, lifestage=NULL, CP=0.5, FA=0.1018428, UA=0.08743042, SDA=0.150819, ED=3983.646) {
	focal_info <- fb4_data[fb4_data$Sci_Name == gsub("_", " ", species),]
	if(!is.null(lifestage)) {
  		focal_info <- focal_info[focal_info$LifeStage == lifestage,]
	}
	if(nrow(focal_info) == 0) {
   		stop(paste("No BEM data found for", species, "with life stage", lifestage))
 	}
	if(nrow(focal_info) > 1) {
	 warning(paste("Multiple BEM data rows found for", species, "with life stage", lifestage, ". Using the first row, but look at all the life stages and specify next time: ", paste(lifestages_from_fb4(species), collapse=", ")))
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
lifestages_from_fb4 <- function(species) {
 focal_info <- fb4_data[fb4_data$Sci_Name == gsub("_", " ", species),]
 if(nrow(focal_info) == 0) {
	 stop(paste("No BEM data found for", species))
  }
 return(unique(focal_info$LifeStage))
}


#' Find closest relatives of a species in the FB4 dataset
#' @param species Scientific name of the species (e.g. "Micropterus salmoides")
#' @return A vector of scientific names of the closest relatives in the FB4 dataset
#' @export
#' 
#' @description This uses the Open Tree of Life to find relatives of a species.
fb4_closest_species <- function(species) {
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
  phy <- suppressWarnings(rotl::tol_induced_subtree(ott_ids = rotl::ott_id(resolved_names), label_format="name")) # suppress because I really don't care about singleton node suppression
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
  sources <- fb4_data$Source[match(final_names, fb4_data$Sci_Name)]
  print(paste("Sources for", paste(final_names, collapse=", "), "are:", paste(sources, collapse=", ")))
  return(final_names) 
}

#' Create a new object with BEM data
#' @param CP proportion of maximum consumption
#' @param CA intercept of allometric mass function
#' @param CB slope of allometric mass function; should be negative bc bigger fish consume less per gram of body mass than smaller fish
#' @param CTM critical thermal maximum (degrees C)
#' @param CTO laboratory temperature preferendum (degrees C)
#' @param CQ approximates a Q10; the rate at which the function increases over relatively low water temperatures
#' @param CTL temperature at which consumption is XXXX (degrees C)
#' @param CK1 temperature at which consumption is XXXX (degrees C)
#' @param CK4 temperature at which consumption is XXXX (degrees C)
#' @param ACT activity coefficient
#' @param RA intercept of allometric activity function
#' @param RB slope of allometric activity function
#' @param RTM critical thermal maximum (degrees C)
#' @param RTO laboratory temperature preferendum (degrees C)
#' @param RQ approximates a Q10; the rate at which the function increases over relatively low water temperatures
#' @param FA fish body mass (grams)
#' @param UA activity (mg/kg-day)
#' @param SDA specific dynamic action (proportion)
#' @param ED energy density (J/kg)
#' @return BEM object
#' @export
bem_new <- function(CP, CA, CB, CTM, CTO, CQ, CTL, CK1, CK4, ACT, RA, RB, RTM, RTO, RQ, FA, UA=0.1, SDA=0.15, ED=3000) {
	result <- c(CP, CA, CB, CTM, CTO, CQ, CTL, CK1, CK4, ACT, RA, RB, RTM, RTO, RQ, FA, UA, SDA, ED)
	names(result) <- c('CP', 'CA', 'CB', 'CTM', 'CTO', 'CQ', 'CTL', 'CK1', 'CK4', 'ACT', 'RA', 'RB', 'RTM', 'RTO', 'RQ', 'FA', 'UA', 'SDA', 'ED')
	result_df <- data.frame(t(as.matrix(result)))
	class(result_df) <- c('BEM', 'data.frame')
	return(result_df)
}

#' Get daily temperature data from monitoring stations
#' 
#' @description This uses the dataRetrieval package to fetch daily temperature data 
#' @param site_code WQP site code (e.g. "USGS-12345678")
#' @param start_date Start date for the data retrieval (e.g. "2020-01-01")
#' @param end_date End date for the data retrieval (e.g. "2020-12-31")
#' @return A data frame with daily temperature data
#' @export
daily_temperature_raw <- function(site_code, start_date, end_date) {
	results <- dataRetrieval::readWQPqw(siteNumbers=site_code, startDate=start_date, endDate=end_date, parameterCd="00010")	
	return(results)
}

#' Interpolate daily temperature data
#' 
#' @description From a raw set of daily temperature data from daily_temperature_raw, this interpolates the temperature data to fill in missing dates. You may want to have raw data before and after the period of interest so that dates within can be interpolated. If the start and end dates are not provided, it will use the min and max dates from the temperature data.
#' @param temperature_data Data frame with daily temperature data from daily_temperature_raw()
#' @param start_date Start date for the interpolation (e.g. "2020-01-01")
#' @param end_date End date for the interpolation (e.g. "2020-12-31")
#' @param minimum_fraction_missing Minimum fraction of data that must be present to do the interpolation; otherwise, it will return NA for all temperatures. Default is 0.8, meaning at least 80% of the data must be present.
#' @return A data frame with interpolated daily temperature data, including date, julian day, and temperature
#' @export
#' 
daily_temperature_interpolate <- function(temperature_data, start_date=NULL, end_date=NULL, minimum_fraction_missing=0.8) {
  # Ensure the date column is in Date format
  temperature_data$ActivityStartDate <- as.Date(temperature_data$ActivityStartDateTime)
  
  # Create a sequence of dates from start_date to end_date
  all_dates <- seq(min(c(temperature_data$ActivityStartDate, as.Date(start_date))), max(c(temperature_data$ActivityStartDate, as.Date(end_date))), by="day")
  
  number_observed_days <- length(unique(temperature_data$ActivityStartDate))
  number_total_days <- length(all_dates)
  fraction_observed <- number_observed_days / number_total_days
  if(fraction_observed < minimum_fraction_missing) {
	warning(paste("Only", round(fraction_observed * 100, 2), "% of the data is present. Returning NA for all temperatures."))
	return(data.frame(date=all_dates, julian=as.numeric(format(all_dates, "%j")), temp=NA))
  }
  
  # Interpolate the temperature data to fill in missing dates
  interpolated_temps <- zoo::na.approx(zoo::zoo(temperature_data$ResultMeasureValue, temperature_data$ActivityStartDate), xout=all_dates, na.rm=FALSE)
  
  # Create a data frame with the interpolated temperatures
  result_df <- data.frame(date=all_dates, julian=as.numeric(format(all_dates, "%j")), temp=as.numeric(interpolated_temps))
  if(!is.null(start_date)) {
	result_df <- result_df[result_df$date >= as.Date(start_date), ]
  }
  if(!is.null(end_date)) {
	 result_df <- result_df[result_df$date <= as.Date(end_date), ]
  }
  return(result_df)
}


#' Get water quality monitoring stations
#' 
#' @description This can  retrieve water quality monitoring stations based on state or bounding box coordinates. It uses the dataRetrieval package to fetch the data.
#' @param state State code (e.g. "CA" for California) to filter stations by state
#' @param min_lon Minimum longitude for bounding box (optional)
#' @param max_lon Maximum longitude for bounding box (optional)
#' @param min_lat Minimum latitude for bounding box (optional)
#' @param max_lat Maximum latitude for bounding box (optional)
#' @return A data frame with water quality monitoring stations, including site number, site name, latitude, longitude, and other relevant information
#' @export
water_stations_find <- function(state=NULL, min_lon=NULL, max_lon=NULL, min_lat=NULL, max_lat=NULL) {
  # Get all water quality monitoring stations
  if (!is.null(state)) {
	siteList <- dataRetrieval::whatNWISsites(stateCd = state, parameterCd = "00010", hasDataTypeCd=c('iv',"dv"))  
  } else { # tile lat and long if needed
   width <- abs(max_lon - min_lon)
   height <- abs(max_lat - min_lat)
   if(width*height <=25) {
	siteList <- dataRetrieval::whatNWISsites(bBox = c(min_lon, min_lat, max_lon, max_lat), parameterCd = "00010", hasDataTypeCd=c('iv',"dv"))
   } else {
		stop("The bounding box is too large -- it must have an area of less than 25 square degrees (yes, this is a weird restriction, from the underlying data provider)")
   }	
  }
  return(siteList)
}

#' Consumption equation 2
#' 
#' @description ADD LATER
#' @param T temperature (degrees C) at which consumption is calculated
#' @param W weight (grams) of fish
#' @param BEM object
#' @return specific consumption rate (grams of food consumed per gram of fish mass per day)
#' @export
consumption2 <- function(T, W, BEM) {
	Y <- log(BEM$CQ)*(BEM$CTM-BEM$CTO+2)
	Z <- log(BEM$CQ)*(BEM$CTM-BEM$CTO)
	X <- (Z^2*(1+(1+40/Y)^0.5)^2)/400
	V <- (BEM$CTM-T)/(BEM$CTM-BEM$CTO)
	fT_C <- V^X*exp(X*(1-V))
	Cmax <- BEM$CA*W^BEM$CB
	C <- Cmax*BEM$CP*fT_C
	return(C)	
}

#' Consumption equation 3
#' 
#' @description ADD LATER
#' @param T temperature (degrees C) at which consumption is calculated
#' @param W weight (grams) of fish
#' @param BEM object
#' @return specific consumption rate (grams of food consumed per gram of fish mass per day)
#' @export
consumption3 <- function(T, W, BEM)
{
  G2 = (1/(BEM$CTL-BEM$CTM))*log((0.98*(1-BEM$CK4))/(BEM$CK4*0.02))
  L2 = exp(G2*(BEM$CTL-T))
  KB = (BEM$CK4*L2)/(1+BEM$CK4*(L2-1))
  G1 = (1/(BEM$CTO-BEM$CQ))* log((0.98*(1-BEM$CK1))/(BEM$CK1*0.02))
  L1 = exp(G1)*(T-BEM$CQ)
  KA = (BEM$CK1*L1)/(1+BEM$CK1*(L1-1))
  fT_C = KA*KB
  Cmax = CA*M^CB
  C = Cmax*BEM$CP*fT_C                         
  return(C)
}

#' Respiration equation 1
#' 
#' @description ADD LATER
#' @param T temperature (degrees C) at which consumption is calculated
#' @param W weight (grams) of fish
#' @param BEM BEM object
respiration1 <- function(T, W, BEM)
{
  fT <- exp(BEM$RQ*T)
  ACTIVITY <- BEM$ACT
  R <- BEM$RA*W^BEM$RB*fT*ACTIVITY
  return(R)
}

#' Respiration equation 2
#' 
#' @description ADD LATER
#' @param T temperature (degrees C) at which consumption is calculated
#' @param W weight (grams) of fish
#' @param BEM BEM object
respiration2 <- function(T, W, BEM)
{
  Y=log(BEM$RQ)*(BEM$RTM-BEM$RTO+2)
  Z=log(BEM$RQ)*(BEM$RTM-BEM$RTO)
  X=(Z^2*(1+(1+40/Y)^0.5)^2)/400
  V=(BEM$RTM-T)/(BEM$RTM-BEM$RTO)
  ACTIVITY=BEM$ACT
  fT=V^X*exp(X*(1-V))
  R=BEM$RA*BEM$W^BEM$RB*fT*ACTIVITY
  return(R)
}

#' Compute growth for a year at a station
#' 
#' @description ADD LATER
#' @param T_vector Temperatures (degrees C), starting at Julian day 1
#' @param BEM BEM object
#' @param starting_weight Starting weight (grams) from Groeschel-Taylor et al
#' @param prey_ED joules per gram of wet mass
#' @param oxycal_coeff joules per gram of oxygen
#' @return data.frame with growth over year
#' @export
single_station_compute <- function(T_vector, BEM, starting_weight=6.382417, prey_ED=3698.0, oxycal_coeff=13560.0) {
	results <- data.frame(julian=sequence(365), temp=T_vector, C1_ins=NA, C2_ins=NA, R1_ins=NA, R2_ins=NA, F1_ins=NA, F2_ins=NA, U1_ins=NA, U2_ins=NA, SDA1_ins=NA, SDA2_ins=NA, W1_ins=NA, W2_ins=NA, W1_cum=NA, W2_cum=NA)
	results[1,] <- 0
	results$W2_cum[1] <- starting_weight
	for (day in 2:nrow(results)) {
		# simulate consumption --> grams of prey
		results$C1_ins[day] <- consumption2(T=results$temp[day], W=results$W2_cum[day-1], BEM=BEM) * results$W2_cum[day-1]
		if(is.na(results$C1_ins[day])) {
			results$C1_ins[day] <- 0
		}
		
		# simulate consumption --> joules of energy
		results$C2_ins[day] <- results$C1_ins[day] * prey_ED  # convert from grams of food to joules with prey energy density parameter
		
		# simulate respiration --> grams of oxygen
		results$R1_ins[day] <- respiration1(T=results$temp[day], W=results$W2_cum[day-1], BEM=BEM) *  results$W2_cum[day-1]
		
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
	results$ProportionStartingWeight <- results$W2_cum / starting_weight
	return(results)
}

#' Loop through species across a station
#' 
#' @description ADD LATER
#' @param T_vector Temperatures (degrees C), starting at Julian day 1
#' @param species Scientific name of the species (e.g. "Micropterus salmoides")
#' @param starting_weight Starting weight (grams) for the species; can be a single value or a vector with one value per organism
#' @param prey_ED joules per gram of wet mass for the species' prey; can be a single value or a vector with one value per organism
#' @param oxycal_coeff joules per gram of oxygen for the species; can be a single value or a vector with one value per organism
#' @return A data.frame with growth results for each species and life stage, including species name, life stage, and cumulative weight over the year
#' @export
loop_species_across_station <- function(T_vector, starting_weight=6.382417, prey_ED=3698.0, oxycal_coeff=13560.0) {
  n_organisms <- nrow(fb4_data)
  results_list <- list()

  for (organism_index in sequence(n_organisms)) {
	focal_species <- fb4_data$Sci_Name[organism_index]
	focal_lifestage <- fb4_data$LifeStage[organism_index]
	focal_common <- fb4_data$Species[organism_index]
	BEM <- bem_from_fb4(focal_species, focal_lifestage)
	actual_starting_weight <- ifelse(length(starting_weight) == 1, starting_weight, starting_weight[organism_index])
	actual_prey_ED <- ifelse(length(prey_ED) == 1, prey_ED, prey_ED[organism_index])
	actual_oxycal_coeff <- ifelse(length(oxycal_coeff) == 1, oxycal_coeff, oxycal_coeff[organism_index])
	weight_result <- single_station_compute(T_vector, BEM, actual_starting_weight, actual_prey_ED, actual_oxycal_coeff)
	weight_result$species <- focal_species
	weight_result$life_stage <- focal_lifestage
	weight_result$organism_index <- organism_index
	weight_result$common_name <- focal_common
	results_list[[organism_index]] <- weight_result
  }
  weight_results_df <- do.call(rbind, results_list)
  return(weight_results_df)
}