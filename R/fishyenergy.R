### TODO

#' Intrinsic bioenergetics parameters from FB4
#'
#' This dataset contains intrinsic bioenergetic parameters for five species.
#'
#' @format A data frame with 5 rows and 27 columns:
#' \describe{
#'   \item{genus_species}{(character) description}
#'   \item{LifeStage}{(character) description}
#'   \item{Source}{(character) description}
#'   \item{lwA}{(numeric) y-intercept for length-weight relationship on log-log scale}
#'   \item{lwB}{(numeric) slope for length-weight relationship on log-log scale}
#'   \item{CEQ}{(integer) consumption equation 1, 2, or 3}
#'   \item{CA}{(numeric) intercept estimate for mass-dependent consumption}
#'   \item{CB}{(numeric) slope estimate for mass-dependent consumption}
#'   \item{CQ}{(numeric) Q10 for temperature-dependent consumption}
#'   \item{CTO}{(numeric) thermal optimum for consumption}
#'   \item{CTM}{(numeric) thermal maximum for consumption}
#'   \item{CTL}{(numeric) temperature-dependent parameter for consumption}
#'   \item{CK1}{(numeric) temperature-dependent parameter for consumption}
#'   \item{CK4}{(numeric) temperature-dependent parameter for consumption}
#'   \item{REQ}{(integer) respiration equation 1 or 2}
#'   \item{RA}{(numeric) intercept estimate for mass-dependent respiration}
#'   \item{RB}{(numeric) slope estimate for mass-dependent respiration}
#'   \item{RQ}{(numeric) Q10 for temperature-dependent respiration}
#'   \item{RTO}{(numeric) thermal optimum for respiration}
#'   \item{RTM}{(numeric) thermal maximum for respiration}
#'   \item{RTL}{(numeric) temperature-dependent parameter for respiration}
#'   \item{RK1}{(numeric) temperature-dependent parameter for respiration}
#'   \item{RK4}{(numeric) temperature-dependent parameter for respiration}
#'   \item{SDA}{(numeric) specific dynamic action}
#'   \item{FA}{(numeric) egestion}
#'   \item{UA}{(numeric) excretion}
#'   \item{pred_ED}{(numeric) energy density of predator (i.e., focal species)}
#'   ...
#' }
#' @source Fish Bioenergetics 4 (Deslauriers et al. 2017).
#' 
"parms_fb4"

#' Daily water temperature for one stream
#'
#' This dataset contains mean daily water temperature for every day of the year 2025 from the USGS gage on the Little River above Townsend, TN (03497300).
#'
#' @format A data frame with 365 rows and 2 columns:
#' \describe{
#'   \item{date}{(date) dates from 1 January 2025 to 31 December 2025}
#'   \item{WT_mean}{(numeric) mean daily water temperature in degrees C}
#'   ...
#' }
#' @source National Water Information System.
#' 
"temps_nwis"

#' Daily water temperature for many lakes
#'
#' This dataset contains modeled maximum daily water temperature for every day of the year 2022 for 8,046 lakes in 15 southeastern United States plus Puerto Rico.
#'
#' @format A data frame with 365 rows and 8,046 columns:
#' \describe{
#'   \item{date}{(date) dates from 1 January 2021 to 31 December 2021; note that temps are the mean from years 2013 to 2021}
#'   \item{Hylak_id}{(numeric) mean daily water temperature in degrees C for all 8046 unique identifier of lake, from HydroLAKES (Messager et al., 2016)}
#'   ...
#' }
#' @source LakeTEMP; Korver et al. 2024. Remote Sensing of Environment. https://doi.org/10.1016/j.rse.2024.114164.
#' 
"temps_lake"

#' Daily water temperature for many streams
#'
#' This dataset contains modeled maximum daily water temperature for every day of the year 2022 for all interconfluence reaches of the Little River, TN.
#'
#' @format A data frame with 365 rows and 515 columns:
#' \describe{
#'   \item{date}{(date) dates from 1 January 2022 to 31 December 2022}
#'   \item{COMID}{(numeric) mean daily water temperature in degrees C for all 513 COMIDs}
#'   ...
#' }
#' @source Southeastern Climate Adaptation Science Center (CASC).
#' 
"temps_casc"

#' Daily water temperature many bay/estuary grid cells
#'
#' This dataset provides daily ocean/estuary surface temperatures along the western coast of the Gulf of Mexico. In fishyenergy, we provide a 1-km resolution raster layer of these temperatures.
#'
#' @format A 3-dimensional array with spatial dimensions of 542 vertical pixels by 410 horizontal pixels and 365 days in the temporal dimension. Note that only 8,862 of the 222,220 spatial pixels represent a bay/estuary and have a temperature value.:
#' \describe{
#'   ...
#' }
#' @source Texas BAYCAST.
#' 
"temps_bayc"

#' Environmental covariates for many lakes
#'
#' This dataset contains geographic and bathymetric covariates for all 8,046 lakes extracted from the lakeTEMP dataset.
#'
#' @format A data frame with 8,046 rows and 10 columns:
#' \describe{
#'   \item{Hylak_id}{(numeric) mean daily water temperature in degrees C for all 8046 unique identifier of lake, from HydroLAKES (Messager et al., 2016)}
#'   \item{center_long}{(numeric) longitude of lake center in decimal degrees}
#'   \item{center_lat}{(numeric) latitude of lake center in decimal degrees}
#'   \item{stat_method}{(numeric) see Korver et al. 2024 for description}
#'   \item{n_obs}{(numeric) see Korver et al. 2024 for description}
#'   \item{intermittency}{(numeric) see Korver et al. 2024 for description}
#'   \item{Lake_type}{(numeric) 1 = natural lake, 2 = reservoir, 3 = natural lake with regulation structure}
#'   \item{Elevation}{(numeric) lake surface elevation in meters above sea level}
#'   \item{Lake_area}{(numeric) lake surface area in square kilometers}
#'   \item{Depth_avg}{(numeric) mean lake depth in meters}
#'   ...
#' }
#' @source LakeTEMP; Korver et al. 2024. Remote Sensing of Environment. https://doi.org/10.1016/j.rse.2024.114164.
#' 
"envir_lake"

#' Environmental covariates for many stream reaches
#'
#' This dataset contains geographic and hydrographic covariates for all 513 interconfluence stream reaches in the Little River, TN.
#'
#' @format A vector shapefile with an attribute table containing 513 rows and 5 columns:
#' \describe{
#'   \item{COMID}{(character) Common Identifier for NHD reaches}
#'   \item{GNIS_NAME}{(character) Stream name from NHD after Geographic Names Information System}
#'   \item{FCode}{(integer) Flowline type from NHD}
#'   \item{TotDASqKM}{(numeric) Catchment area draining to the interconfluence stream reach in square kilometers}
#'   \item{ElevCat}{(numeric) Median elevation of the interconfluence stream reach in meters above sea level}
#'   ...
#' }
#' @source Southeastern Climate Adaptation Science Center (CASC), National Hydrography Dataset (NHD), StreamCat
#' 
"envir_casc"

#' Temporally varying parameters
#'
#' This dataset contains bioenergetics parameters that potentially vary each day of a 365-day calendar year.
#'
#' @format A data frame with 365 rows and 5 columns:
#' \describe{
#'   \item{date}{(date) dates from 1 to 365}
#'   \item{CP}{(numeric) Proportion of consumption; set to 1.0 assumes prey is unlimited}
#'   \item{ACT}{(numeric) Activity multiplier; set to 1.0 assumes predator (i.e., focal fish) always is at rest}
#'   \item{gsi_f}{(numeric) Gonadosomatic index for females; set to 0 assumes reproductive immaturity}
#'   \item{gsi_m}{(numeric) Gonadosomatic index for males; set to 0 assumes reproductive immaturity}
#'   ...
#' }
#' @source n/a.
#' 
"parms_temporal_DEFAULT"

#' Consumption equation 1
#' 
#' @description consumption equation 1 from Hanson et al. 1997
#' @param M mass (grams) of fish
#' @param T temperature (degrees C) at which consumption is calculated
#' @param CP proportion of maximum consumption (default is 1.0)
#' @param parms.intrinsic A parms.intrinsic object
#' @return specific consumption rate (grams of food consumed per gram of fish mass per day)
#' @export
<<<<<<< Updated upstream
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
bem_new <- function(CP, CA, CB, CTM, CTO, CQ, CTL, CK1, CK4, ACT, RA, RB, RTM, RTO, RQ, FA, UA, SDA, ED) {
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
=======
consumption1 <- function(M, T, CP= 1.0, parms.intrinsic)
{
  fT_C <- exp(parms.intrinsic$CQ*T)
  Cmax <- parms.intrinsic$CA*M^parms.intrinsic$CB
  C <- Cmax*CP*fT_C
  return(C)
>>>>>>> Stashed changes
}

#' Consumption equation 2
#' 
#' @description consumption equation 2 from Hanson et al. 1997
#' @param M mass (grams) of fish
#' @param T temperature (degrees C) at which consumption is calculated
#' @param CP proportion of maximum consumption (default is 1.0)
#' @param parms.intrinsic A parms.intrinsic object
#' @return specific consumption rate (grams of food consumed per gram of fish mass per day)
#' @export
consumption2 <- function(M, T, CP = 1.0, parms.intrinsic)
{
  Y <- log(parms.intrinsic$CQ)*(parms.intrinsic$CTM-parms.intrinsic$CTO+2)
  Z <- log(parms.intrinsic$CQ)*(parms.intrinsic$CTM-parms.intrinsic$CTO)
  X <- (Z^2*(1+(1+40/Y)^0.5)^2)/400
  V <- (parms.intrinsic$CTM-T)/(parms.intrinsic$CTM-parms.intrinsic$CTO)
  fT_C <- V^X*exp(X*(1-V))
  Cmax <- parms.intrinsic$CA*M^parms.intrinsic$CB
  C <- Cmax*CP*fT_C
  return(C)	
}

#' Consumption equation 3
#' 
#' @description Consumption equation 3 from Hanson et al. 1997
#' @param M Mass (grams) of fish
#' @param T Temperature (degrees C) at which consumption is calculated
#' @param CP Proportion of maximum consumption (default is 1.0)
#' @param parms.intrinsic A parms.intrinsic object
#' @return Specific consumption rate (grams of food consumed per gram of fish mass per day)
#' @export
consumption3 <- function(M, T, CP= 1.0, parms.intrinsic)
{
  G2 <- (1/(parms.intrinsic$CTL-parms.intrinsic$CTM))*log((0.98*(1-parms.intrinsic$CK4))/(parms.intrinsic$CK4*0.02))
  L2 <- exp(G2*(parms.intrinsic$CTL-T))
  KB <- (parms.intrinsic$CK4*L2)/(1+parms.intrinsic$CK4*(L2-1))
  G1 <- (1/(parms.intrinsic$CTO-parms.intrinsic$CQ))* log((0.98*(1-parms.intrinsic$CK1))/(parms.intrinsic$CK1*0.02))
  L1 <- exp(G1)*(T-parms.intrinsic$CQ)
  KA <- (parms.intrinsic$CK1*L1)/(1+parms.intrinsic$CK1*(L1-1))
  fT_C <- KA*KB
  Cmax <- parms.intrinsic$CA*M^parms.intrinsic$CB
  C <- Cmax*CP*fT_C                         
  return(C)
}

#' Respiration equation 1
#' 
#' @description Respiration equation 1 from Hanson et al. 1997
#' @param M Mass (grams) of fish
#' @param T Temperature (degrees C) at which respiration is calculated
#' @param ACT Activity multiplier (default is 1.0)
#' @param parms.intrinsic A parms.intrinsic object
#' @return Specific respiration rate (grams of oxygen per gram of fish mass per day)
#' @export
respiration1 <- function(M, T, ACT= 1.0, parms.intrinsic)
{
  fT_R <- exp(parms.intrinsic$RQ*T)
  Rrest <- parms.intrinsic$RA*M^parms.intrinsic$RB*fT_R
  R <- Rrest*ACT
  return(R)
}

#' Respiration equation 2
#' 
#' @description Respiration equation 2 from Hanson et al. 1997
#' @param M Mass (grams) of fish
#' @param T Temperature (degrees C) at which respiration is calculated
#' @param ACT Activity multiplier (default is 1.0)
#' @param parms.intrinsic A parms.intrinsic object
#' @return Specific respiration rate (grams of oxygen per gram of fish mass per day)
#' @export
respiration2 <- function(M, T, ACT= 1.0, parms.intrinsic)
{
  Y <- log(parms.intrinsic$RQ)*(parms.intrinsic$RTM-parms.intrinsic$RTO+2)
  Z <- log(parms.intrinsic$RQ)*(parms.intrinsic$RTM-parms.intrinsic$RTO)
  X <- (Z^2*(1+(1+40/Y)^0.5)^2)/400
  V <- (parms.intrinsic$RTM-T)/(parms.intrinsic$RTM-parms.intrinsic$RTO)
  fT_R <- V^X*exp(X*(1-V))
  Rrest <- parms.intrinsic$RA*M^parms.intrinsic$RB*fT_R
  R <- Rrest*ACT
  return(R)
}

#' Plot temperature dependent curves
#' 
#' @description Plot the thermal reaction norms for consumption, respiration, and other intrinsic rates
#' @param M Mass (grams) of fish
#' @param T Vector of temperatures representing the temperature range to be plotted
#' @param CP Proportion of maximum consumption (default is 1.0)
#' @param ACT Activity multiplier (default is 1.0)
#' @param prey_ED Prey energy density in joules per gram of wet mass (default is 4000)
#' @param parms.intrinsic A parms.intrinsic object
#' @param C_eq Specify consumption equation 1, 2, or 3 from Hanson et al. 1997
#' @param R_eq Specify respiration equation 1 or 2 from Hanson et al. 1997
#' @return A plot of temperature dependent rates
#' @export
<<<<<<< Updated upstream
single_station_compute <- function(T_vector, BEM, starting_weight=6.382417, prey_ED=3698.0, oxycal_coeff=13560.0) {
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
	results$ProportionStartingWeight <- results$W2_cum / starting_weight
	return(results)
=======
bem_curve <- function(M, T, CP = 1.0, ACT = 1.0, prey_ED = 4000, parms.intrinsic, C_eq, R_eq)       
{
  # consumption curve --> grams of prey and joules of energy
  if(C_eq == 1){
    C1_curve <- consumption1(M, T, CP, parms.intrinsic)
    C2_curve = C1_curve * prey_ED
  }
  else if(C_eq == 2){
    C1_curve <- consumption2(M, T, CP, parms.intrinsic)
    C2_curve = C1_curve * prey_ED
  }
  else if (C_eq == 3){
    C1_curve <- consumption3(M, T, CP, parms.intrinsic)
    C2_curve = C1_curve * prey_ED
  }
  else{
    stop("Consumption and/or respiration equation(s) not specified")
  }
  
  # respiration curve --> grams of oxygen and joules of energy
  if(R_eq == 1){
    R1_curve <- respiration1(M, T, ACT, parms.intrinsic)
    R2_curve = R1_curve * 13560.0
  }
  else if(R_eq == 2){
    R1_curve <- respiration2(M, T, ACT, parms.intrinsic)
    R2_curve = R1_curve * 13560.0
  }
  else{
    stop("Consumption and/or respiration equation(s) not specified")
  }
  
  # egestion curve --> grams and joules of prey; assume egestion is a constant proportion of consumption
  F1_curve = C1_curve * parms.intrinsic$FA
  F2_curve = C2_curve * parms.intrinsic$FA
  
  # excretion curve --> grams and joules of prey; assume excretion is a constant proportion of consumption
  U1_curve = C1_curve * parms.intrinsic$UA
  U2_curve = C2_curve * parms.intrinsic$UA
  
  # SDA curve --> grams and joules of prey; assume specific dynamic action is a constant proportion of consumption
  SDA1_curve = parms.intrinsic$SDA * (C1_curve - F1_curve)
  SDA2_curve = parms.intrinsic$SDA * (C2_curve - F2_curve)
  
  # mass curve --> grams and joules of predator
  M1_curve <- C1_curve - (R1_curve + F1_curve + U1_curve + SDA1_curve)
  M2_curve <- M1_curve * (1/parms.intrinsic$pred_ED)
  
  # combine curves into a dataframe
  dframe.C1_curve <- data.frame(T, rep(1,length(C1_curve)), rep("Consumption",length(C1_curve)), C1_curve); colnames(dframe.C1_curve) <- c("temp","units","Curve","value")
  dframe.C2_curve <- data.frame(T, rep(2,length(C2_curve)), rep("Consumption",length(C2_curve)), C2_curve); colnames(dframe.C2_curve) <- c("temp","units","Curve","value")
  dframe.R1_curve <- data.frame(T, rep(1,length(R1_curve)), rep("Respiration",length(R1_curve)), R1_curve); colnames(dframe.R1_curve) <- c("temp","units","Curve","value")
  dframe.R2_curve <- data.frame(T, rep(2,length(R2_curve)), rep("Respiration",length(R2_curve)), R2_curve); colnames(dframe.R2_curve) <- c("temp","units","Curve","value")
  dframe.F1_curve <- data.frame(T, rep(1,length(F1_curve)), rep("Egestion",length(F1_curve)), F1_curve); colnames(dframe.F1_curve) <- c("temp","units","Curve","value")
  dframe.F2_curve <- data.frame(T, rep(2,length(F2_curve)), rep("Egestion",length(F2_curve)), F2_curve); colnames(dframe.F2_curve) <- c("temp","units","Curve","value")
  dframe.U1_curve <- data.frame(T, rep(1,length(U1_curve)), rep("Excretion",length(U1_curve)), U1_curve); colnames(dframe.U1_curve) <- c("temp","units","Curve","value")
  dframe.U2_curve <- data.frame(T, rep(2,length(U2_curve)), rep("Excretion",length(U2_curve)), U2_curve); colnames(dframe.U2_curve) <- c("temp","units","Curve","value")
  dframe.SDA1_curve <- data.frame(T, rep(1,length(SDA1_curve)), rep("SDA",length(SDA1_curve)), SDA1_curve); colnames(dframe.SDA1_curve) <- c("temp","units","Curve","value")
  dframe.SDA2_curve <- data.frame(T, rep(2,length(SDA2_curve)), rep("SDA",length(SDA2_curve)), SDA2_curve); colnames(dframe.SDA2_curve) <- c("temp","units","Curve","value")
  dframe.M1_curve <- data.frame(T, rep(1,length(M1_curve)), rep("Surplus",length(M1_curve)), M1_curve); colnames(dframe.M1_curve) <- c("temp","units","Curve","value")
  dframe.M2_curve <- data.frame(T, rep(2,length(M2_curve)), rep("Surplus",length(M2_curve)), M2_curve); colnames(dframe.M2_curve) <- c("temp","units","Curve","value")
  dframe.curves <- rbind(dframe.C1_curve,
                         dframe.C2_curve,
                         dframe.R1_curve,
                         dframe.R2_curve,
                         dframe.F1_curve,
                         dframe.F2_curve,
                         dframe.U1_curve,
                         dframe.U2_curve,
                         dframe.SDA1_curve,
                         dframe.SDA2_curve,
                         dframe.M1_curve,
                         dframe.M2_curve)
  
  dframe.curves$value <- ifelse(dframe.curves$value < 0, 0 ,dframe.curves$value)
  dframe.curves$value <- ifelse(is.na(dframe.curves$value), 0 ,dframe.curves$value)
  dframe.curves$Curve <- factor(dframe.curves$Curve, levels = c("Surplus","Consumption","Respiration","SDA","Egestion","Excretion"))
  
  # plot curves
  ggplot(NULL) + 
    theme(plot.tag.position = c(0.15, 0.02)) +
    theme_classic() +
    labs(x = "Temperature (°C)",
         y = "Energy (joules/g/d)") +
    scale_size_manual(values = c("Surplus" = 2, 
                                 "Consumption" = 1, 
                                 "Respiration" = 1, 
                                 "Excretion" = 1, 
                                 "Egestion" = 1, 
                                 "SDA" = 1)) +
    scale_color_manual(values = c("Surplus" = "black", 
                                  "Consumption" = "blue", 
                                  "Respiration" = "red", 
                                  "Excretion" = "gray65", 
                                  "Egestion" = "gray50", 
                                  "SDA" = "gray35")) +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0, max(dframe.curves[dframe.curves$units == 1,"value"])*1.05)) +
    geom_line(aes(x = temp, y = value, color = Curve, size = Curve), dframe.curves[dframe.curves$units == 1,])
>>>>>>> Stashed changes
}

#' Solve daily energy budgets to simulate growth over a calendar year
#' 
#' @description the Wisconsin bioenergetics model from Hanson et al. 1997
#' @param start_M2 (default is 100 grams)
#' @param temperature a dataframe populated with a time series of mean daily water temperature (degrees C) of a habitat patch (default is 20 degrees C)
#' @param parms.intrinsic A parms.intrinsic object
#' @param parms.temporal a dataframe populated with a time series of intrinsic (e.g., GSI) and extrinsic (e.g., prey energy density) biological parameters
#' @param C_eq Specify consumption equation 1, 2, or 3 from Hanson et al. 1997
#' @param R_eq Specify respiration equation 1 or 2 from Hanson et al. 1997
#' @return dataframe populated with output parameters (columns) for time series of projected days (rows)
#' @export
bem_grow <- function(start_M2 = 100, temperature = 20, parms.intrinsic, parms.temporal, C_eq, R_eq)
{
  dframe.sim_parms <- temperature
  dframe.sim_parms$pred_ED = parms.intrinsic$pred_ED
  
  # add temporal factors
  dframe.sim_parms$gsi_m = parms.temporal$gsi_m
  dframe.sim_parms$gsi_f = parms.temporal$gsi_f
  dframe.sim_parms$CP = parms.temporal$CP
  dframe.sim_parms$ACT = parms.temporal$ACT
  dframe.sim_parms$prey_ED = parms.temporal$prey_ED
  
  # add parm columns with NAs as fillers
  dframe.sim_parms$C1_ins = NA
  dframe.sim_parms$C2_ins = NA
  dframe.sim_parms$R1_ins = NA
  dframe.sim_parms$R2_ins = NA
  dframe.sim_parms$F1_ins = NA
  dframe.sim_parms$F2_ins = NA
  dframe.sim_parms$U1_ins = NA
  dframe.sim_parms$U2_ins = NA
  dframe.sim_parms$SDA1_ins = NA
  dframe.sim_parms$SDA2_ins = NA
  dframe.sim_parms$MS1_ins = NA
  dframe.sim_parms$MS2_ins = NA
  dframe.sim_parms$MG1_ins = NA
  dframe.sim_parms$MG2_ins = NA
  dframe.sim_parms$M1_ins = NA
  dframe.sim_parms$M2_ins = NA
  dframe.sim_parms$C1_cum = NA
  dframe.sim_parms$C2_cum = NA
  dframe.sim_parms$R1_cum = NA
  dframe.sim_parms$R2_cum = NA
  dframe.sim_parms$F1_cum = NA
  dframe.sim_parms$F2_cum = NA
  dframe.sim_parms$U1_cum = NA
  dframe.sim_parms$U2_cum = NA
  dframe.sim_parms$SDA1_cum = NA
  dframe.sim_parms$SDA2_cum = NA
  dframe.sim_parms$MS1_cum = NA
  dframe.sim_parms$MS2_cum = NA
  dframe.sim_parms$MG1_cum = NA
  dframe.sim_parms$MG2_cum = NA
  dframe.sim_parms$M1_cum = NA
  dframe.sim_parms$M2_cum = NA
  
  # add starting weight in first row; add zeros to all other first rows
  dframe.sim_parms[1,c("C1_ins")] = 0
  dframe.sim_parms[1,c("C2_ins")] = 0
  dframe.sim_parms[1,c("R1_ins")] = 0
  dframe.sim_parms[1,c("R2_ins")] = 0
  dframe.sim_parms[1,c("F1_ins")] = 0
  dframe.sim_parms[1,c("F2_ins")] = 0
  dframe.sim_parms[1,c("U1_ins")] = 0
  dframe.sim_parms[1,c("U2_ins")] = 0
  dframe.sim_parms[1,c("SDA1_ins")] = 0
  dframe.sim_parms[1,c("SDA2_ins")] = 0
  dframe.sim_parms[1,c("MS1_ins")] = 0
  dframe.sim_parms[1,c("MS2_ins")] = 0
  dframe.sim_parms[1,c("MG1_ins")] = 0
  dframe.sim_parms[1,c("MG2_ins")] = 0
  dframe.sim_parms[1,c("M1_ins")] = 0
  dframe.sim_parms[1,c("M2_ins")] = 0
  dframe.sim_parms[1,c("C1_cum")] = 0
  dframe.sim_parms[1,c("C2_cum")] = 0
  dframe.sim_parms[1,c("R1_cum")] = 0
  dframe.sim_parms[1,c("R2_cum")] = 0
  dframe.sim_parms[1,c("F1_cum")] = 0
  dframe.sim_parms[1,c("F2_cum")] = 0
  dframe.sim_parms[1,c("U1_cum")] = 0
  dframe.sim_parms[1,c("U2_cum")] = 0
  dframe.sim_parms[1,c("SDA1_cum")] = 0
  dframe.sim_parms[1,c("SDA2_cum")] = 0
  dframe.sim_parms[1,c("MS1_cum")] = 0
  dframe.sim_parms[1,c("MS2_cum")] = start_M2 * (1 - parms.temporal[1,"gsi_f"])               # set initial weight of somatic tissue
  dframe.sim_parms[1,c("MG1_cum")] = 0
  dframe.sim_parms[1,c("MG2_cum")] = start_M2 * parms.temporal[1,"gsi_f"]                     # set initial weight of gonadal tissue
  dframe.sim_parms[1,c("M1_cum")] = 0
  dframe.sim_parms[1,c("M2_cum")] = start_M2                                                  # set initial weight total (somatic + gonadal)
  dframe.sim_parms[1,c("L_cum")] = (start_M2 / parms.intrinsic$lwA)^(1 / parms.intrinsic$lwB) # calculate total length
  
  if(C_eq == 1 & R_eq == 1){
    for(i in 2:nrow(dframe.sim_parms)){
      dframe.sim_parms[i,c("C1_ins")] = ifelse(is.nan(consumption1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"CP"],                                  
                                                                   parms.intrinsic)),
                                               0,                                                                                              # ** if temp > CTM, then consumption is zero
                                               consumption1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"CP"],                                         
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],                           
                                                                   parms.intrinsic)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("C2_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.temporal[i,"prey_ED"]                                         # convert from grams of food to joules with prey energy density parameter
      dframe.sim_parms[i,c("R2_ins")] = dframe.sim_parms[i,c("R1_ins")] * 13560.0                                                              # convert from grams of oxygen to joules with oxycalorific coefficient
      dframe.sim_parms[i,c("F1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("F2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("U1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption
      dframe.sim_parms[i,c("U2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption 
      dframe.sim_parms[i,c("SDA1_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C1_ins")] - dframe.sim_parms[i,c("F1_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("SDA2_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C2_ins")] - dframe.sim_parms[i,c("F2_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("MS1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (1- parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MG1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MS2_ins")] = dframe.sim_parms[i,c("MS1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MG2_ins")] = dframe.sim_parms[i,c("MG1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MS1_cum")] = dframe.sim_parms[i-1,c("MS1_cum")] + dframe.sim_parms[i,c("MS1_ins")]
      dframe.sim_parms[i,c("MG1_cum")] = dframe.sim_parms[i-1,c("MG1_cum")] + dframe.sim_parms[i,c("MG1_ins")]
      dframe.sim_parms[i,c("MS2_cum")] = dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i,c("MS2_ins")]
      dframe.sim_parms[i,c("MG2_cum")] = dframe.sim_parms[i-1,c("MG2_cum")] + dframe.sim_parms[i,c("MG2_ins")]
      dframe.sim_parms[i,c("L_cum")] = (dframe.sim_parms[i-1,c("MS2_cum")] / parms.intrinsic$lwA)^(1 / parms.intrinsic$lwB)
      dframe.sim_parms[i,c("L_cum")] = ifelse(dframe.sim_parms[i,c("L_cum")] < dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i,c("L_cum")])
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3)) 
      dframe.sim_parms[i,c("Surv")]<- ifelse(dframe.sim_parms[i,c("K")] <= 0.4, 0, 1)}
  }
  else if(C_eq == 2 & R_eq == 1){
    for(i in 2:nrow(dframe.sim_parms)){
      dframe.sim_parms[i,c("C1_ins")] = ifelse(is.nan(consumption2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"CP"],                                  
                                                                   parms.intrinsic)),
                                               0,                                                                                              # ** if temp > CTM, then consumption is zero
                                               consumption2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"CP"],                                         
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],                           
                                                                   parms.intrinsic)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("C2_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.temporal[i,"prey_ED"]                                         # convert from grams of food to joules with prey energy density parameter
      dframe.sim_parms[i,c("R2_ins")] = dframe.sim_parms[i,c("R1_ins")] * 13560.0                                                              # convert from grams of oxygen to joules with oxycalorific coefficient
      dframe.sim_parms[i,c("F1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("F2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("U1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption
      dframe.sim_parms[i,c("U2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption 
      dframe.sim_parms[i,c("SDA1_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C1_ins")] - dframe.sim_parms[i,c("F1_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("SDA2_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C2_ins")] - dframe.sim_parms[i,c("F2_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("MS1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (1- parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MG1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MS2_ins")] = dframe.sim_parms[i,c("MS1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MG2_ins")] = dframe.sim_parms[i,c("MG1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MS1_cum")] = dframe.sim_parms[i-1,c("MS1_cum")] + dframe.sim_parms[i,c("MS1_ins")]
      dframe.sim_parms[i,c("MG1_cum")] = dframe.sim_parms[i-1,c("MG1_cum")] + dframe.sim_parms[i,c("MG1_ins")]
      dframe.sim_parms[i,c("MS2_cum")] = dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i,c("MS2_ins")]
      dframe.sim_parms[i,c("MG2_cum")] = dframe.sim_parms[i-1,c("MG2_cum")] + dframe.sim_parms[i,c("MG2_ins")]
      dframe.sim_parms[i,c("L_cum")] = (dframe.sim_parms[i-1,c("MS2_cum")] / parms.intrinsic$lwA)^(1 / parms.intrinsic$lwB)
      dframe.sim_parms[i,c("L_cum")] = ifelse(dframe.sim_parms[i,c("L_cum")] < dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i,c("L_cum")])
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3)) 
      dframe.sim_parms[i,c("Surv")]<- ifelse(dframe.sim_parms[i,c("K")] < 0.4, 0, 1)}
    
  }
  else if(C_eq == 3 & R_eq == 1){
    for(i in 2:nrow(dframe.sim_parms)){
      dframe.sim_parms[i,c("C1_ins")] = ifelse(is.nan(consumption3((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"CP"],                                  
                                                                   parms.intrinsic)),
                                               0,                                                                                              # ** if temp > CTM, then consumption is zero
                                               consumption3((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"CP"],                                         
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],                           
                                                                   parms.intrinsic)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("C2_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.temporal[i,"prey_ED"]                                         # convert from grams of food to joules with prey energy density parameter
      dframe.sim_parms[i,c("R2_ins")] = dframe.sim_parms[i,c("R1_ins")] * 13560.0                                                              # convert from grams of oxygen to joules with oxycalorific coefficient
      dframe.sim_parms[i,c("F1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("F2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("U1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption
      dframe.sim_parms[i,c("U2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption 
      dframe.sim_parms[i,c("SDA1_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C1_ins")] - dframe.sim_parms[i,c("F1_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("SDA2_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C2_ins")] - dframe.sim_parms[i,c("F2_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("MS1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (1- parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MG1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MS2_ins")] = dframe.sim_parms[i,c("MS1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MG2_ins")] = dframe.sim_parms[i,c("MG1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MS1_cum")] = dframe.sim_parms[i-1,c("MS1_cum")] + dframe.sim_parms[i,c("MS1_ins")]
      dframe.sim_parms[i,c("MG1_cum")] = dframe.sim_parms[i-1,c("MG1_cum")] + dframe.sim_parms[i,c("MG1_ins")]
      dframe.sim_parms[i,c("MS2_cum")] = dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i,c("MS2_ins")]
      dframe.sim_parms[i,c("MG2_cum")] = dframe.sim_parms[i-1,c("MG2_cum")] + dframe.sim_parms[i,c("MG2_ins")]
      dframe.sim_parms[i,c("L_cum")] = (dframe.sim_parms[i-1,c("MS2_cum")] / parms.intrinsic$lwA)^(1 / parms.intrinsic$lwB)
      dframe.sim_parms[i,c("L_cum")] = ifelse(dframe.sim_parms[i,c("L_cum")] < dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i,c("L_cum")])
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3)) 
      dframe.sim_parms[i,c("Surv")]<- ifelse(dframe.sim_parms[i,c("K")] < 0.4, 0, 1)}
    
  }
  else if(C_eq == 1 & R_eq == 2){
    for(i in 2:nrow(dframe.sim_parms)){
      dframe.sim_parms[i,c("C1_ins")] = ifelse(is.nan(consumption1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"CP"],                                  
                                                                   parms.intrinsic)),
                                               0,                                                                                              # ** if temp > CTM, then consumption is zero
                                               consumption1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"CP"],                                         
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],                           
                                                                   parms.intrinsic)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("C2_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.temporal[i,"prey_ED"]                                         # convert from grams of food to joules with prey energy density parameter
      dframe.sim_parms[i,c("R2_ins")] = dframe.sim_parms[i,c("R1_ins")] * 13560.0                                                              # convert from grams of oxygen to joules with oxycalorific coefficient
      dframe.sim_parms[i,c("F1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("F2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("U1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption
      dframe.sim_parms[i,c("U2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption 
      dframe.sim_parms[i,c("SDA1_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C1_ins")] - dframe.sim_parms[i,c("F1_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("SDA2_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C2_ins")] - dframe.sim_parms[i,c("F2_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("MS1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (1- parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MG1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MS2_ins")] = dframe.sim_parms[i,c("MS1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MG2_ins")] = dframe.sim_parms[i,c("MG1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MS1_cum")] = dframe.sim_parms[i-1,c("MS1_cum")] + dframe.sim_parms[i,c("MS1_ins")]
      dframe.sim_parms[i,c("MG1_cum")] = dframe.sim_parms[i-1,c("MG1_cum")] + dframe.sim_parms[i,c("MG1_ins")]
      dframe.sim_parms[i,c("MS2_cum")] = dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i,c("MS2_ins")]
      dframe.sim_parms[i,c("MG2_cum")] = dframe.sim_parms[i-1,c("MG2_cum")] + dframe.sim_parms[i,c("MG2_ins")]
      dframe.sim_parms[i,c("L_cum")] = (dframe.sim_parms[i-1,c("MS2_cum")] / parms.intrinsic$lwA)^(1 / parms.intrinsic$lwB)
      dframe.sim_parms[i,c("L_cum")] = ifelse(dframe.sim_parms[i,c("L_cum")] < dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i,c("L_cum")])
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3)) 
      dframe.sim_parms[i,c("Surv")]<- ifelse(dframe.sim_parms[i,c("K")] < 0.4, 0, 1)}
    
  }
  else if(C_eq == 2 & R_eq == 2){
    for(i in 2:nrow(dframe.sim_parms)){
      dframe.sim_parms[i,c("C1_ins")] = ifelse(is.nan(consumption2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"CP"],                                  
                                                                   parms.intrinsic)),
                                               0,                                                                                              # ** if temp > CTM, then consumption is zero
                                               consumption2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"CP"],                                         
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],                           
                                                                   parms.intrinsic)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("C2_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.temporal[i,"prey_ED"]                                         # convert from grams of food to joules with prey energy density parameter
      dframe.sim_parms[i,c("R2_ins")] = dframe.sim_parms[i,c("R1_ins")] * 13560.0                                                              # convert from grams of oxygen to joules with oxycalorific coefficient
      dframe.sim_parms[i,c("F1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("F2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("U1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption
      dframe.sim_parms[i,c("U2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption 
      dframe.sim_parms[i,c("SDA1_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C1_ins")] - dframe.sim_parms[i,c("F1_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("SDA2_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C2_ins")] - dframe.sim_parms[i,c("F2_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("MS1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (1- parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MG1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MS2_ins")] = dframe.sim_parms[i,c("MS1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MG2_ins")] = dframe.sim_parms[i,c("MG1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MS1_cum")] = dframe.sim_parms[i-1,c("MS1_cum")] + dframe.sim_parms[i,c("MS1_ins")]
      dframe.sim_parms[i,c("MG1_cum")] = dframe.sim_parms[i-1,c("MG1_cum")] + dframe.sim_parms[i,c("MG1_ins")]
      dframe.sim_parms[i,c("MS2_cum")] = dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i,c("MS2_ins")]
      dframe.sim_parms[i,c("MG2_cum")] = dframe.sim_parms[i-1,c("MG2_cum")] + dframe.sim_parms[i,c("MG2_ins")]
      dframe.sim_parms[i,c("L_cum")] = (dframe.sim_parms[i-1,c("MS2_cum")] / parms.intrinsic$lwA)^(1 / parms.intrinsic$lwB)
      dframe.sim_parms[i,c("L_cum")] = ifelse(dframe.sim_parms[i,c("L_cum")] < dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i,c("L_cum")])
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3)) 
      dframe.sim_parms[i,c("Surv")]<- ifelse(dframe.sim_parms[i,c("K")] < 0.4, 0, 1)}
    
  }
  else if(C_eq == 3 & R_eq == 2){
    for(i in 2:nrow(dframe.sim_parms)){
      dframe.sim_parms[i,c("C1_ins")] = ifelse(is.nan(consumption3((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"CP"],                                  
                                                                   parms.intrinsic)),
                                               0,                                                                                              # ** if temp > CTM, then consumption is zero
                                               consumption3((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"CP"],                                         
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],                           
                                                                   parms.intrinsic)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration2((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])
      
      dframe.sim_parms[i,c("C2_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.temporal[i,"prey_ED"]                                         # convert from grams of food to joules with prey energy density parameter
      dframe.sim_parms[i,c("R2_ins")] = dframe.sim_parms[i,c("R1_ins")] * 13560.0                                                              # convert from grams of oxygen to joules with oxycalorific coefficient
      dframe.sim_parms[i,c("F1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("F2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$FA                                                   # assume egestion is a constant proportion of consumption
      dframe.sim_parms[i,c("U1_ins")] = dframe.sim_parms[i,c("C1_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption
      dframe.sim_parms[i,c("U2_ins")] = dframe.sim_parms[i,c("C2_ins")] * parms.intrinsic$UA                                                   # assume excretion is a constant proportion of consumption 
      dframe.sim_parms[i,c("SDA1_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C1_ins")] - dframe.sim_parms[i,c("F1_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("SDA2_ins")] = parms.intrinsic$SDA * (dframe.sim_parms[i,c("C2_ins")] - dframe.sim_parms[i,c("F2_ins")])            # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
      dframe.sim_parms[i,c("MS1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (1- parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MG1_ins")] = (dframe.sim_parms[i,c("C2_ins")] - (dframe.sim_parms[i,c("R2_ins")] + dframe.sim_parms[i,c("F2_ins")] + dframe.sim_parms[i,c("U2_ins")] + dframe.sim_parms[i,c("SDA2_ins")])) * (parms.temporal[i,"gsi_f"])
      dframe.sim_parms[i,c("MS2_ins")] = dframe.sim_parms[i,c("MS1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MG2_ins")] = dframe.sim_parms[i,c("MG1_ins")] * (1/dframe.sim_parms[i,"pred_ED"])                                   # convert from joules to grams of body mass with predator energy density parameter
      dframe.sim_parms[i,c("MS1_cum")] = dframe.sim_parms[i-1,c("MS1_cum")] + dframe.sim_parms[i,c("MS1_ins")]
      dframe.sim_parms[i,c("MG1_cum")] = dframe.sim_parms[i-1,c("MG1_cum")] + dframe.sim_parms[i,c("MG1_ins")]
      dframe.sim_parms[i,c("MS2_cum")] = dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i,c("MS2_ins")]
      dframe.sim_parms[i,c("MG2_cum")] = dframe.sim_parms[i-1,c("MG2_cum")] + dframe.sim_parms[i,c("MG2_ins")]
      dframe.sim_parms[i,c("L_cum")] = (dframe.sim_parms[i-1,c("MS2_cum")] / parms.intrinsic$lwA)^(1 / parms.intrinsic$lwB)
      dframe.sim_parms[i,c("L_cum")] = ifelse(dframe.sim_parms[i,c("L_cum")] < dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i-1,c("L_cum")], dframe.sim_parms[i,c("L_cum")])
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3)) 
      dframe.sim_parms[i,c("Surv")]<- ifelse(dframe.sim_parms[i,c("K")] < 0.4, 0, 1)}
    
  }
  else{
    stop("Consumption and/or respiration equation(s) not specified")
  }
  
  # compute total weights (somatic + gonadal)
  dframe.sim_parms$M1_ins = (dframe.sim_parms$MS1_ins + dframe.sim_parms$MG1_ins)
  dframe.sim_parms$M2_ins = (dframe.sim_parms$MS2_ins + dframe.sim_parms$MG2_ins)
  dframe.sim_parms$M1_cum = (dframe.sim_parms$MS1_cum + dframe.sim_parms$MG1_cum)
  dframe.sim_parms$M2_cum = (dframe.sim_parms$MS2_cum + dframe.sim_parms$MG2_cum)
  
  # compute remaining cumulative parameters
  dframe.sim_parms$C1_cum = cumsum(dframe.sim_parms$C1_ins)
  dframe.sim_parms$C2_cum = cumsum(dframe.sim_parms$C2_ins)
  dframe.sim_parms$R1_cum = cumsum(dframe.sim_parms$R1_ins)
  dframe.sim_parms$R2_cum = cumsum(dframe.sim_parms$R2_ins)
  dframe.sim_parms$F1_cum = cumsum(dframe.sim_parms$F1_ins)
  dframe.sim_parms$F2_cum = cumsum(dframe.sim_parms$F2_ins)
  dframe.sim_parms$U1_cum = cumsum(dframe.sim_parms$U1_ins)
  dframe.sim_parms$U2_cum = cumsum(dframe.sim_parms$U2_ins)
  dframe.sim_parms$SDA1_cum = cumsum(dframe.sim_parms$SDA1_ins)
  dframe.sim_parms$SDA2_cum = cumsum(dframe.sim_parms$SDA2_ins)
  
  # reduce significant digits to save disk space
  dframe.sim_parms$WT_mean <- round(dframe.sim_parms$WT_mean, digits = 1)
  dframe.sim_parms$gsi_m <- round(dframe.sim_parms$gsi_m, digits = 3)
  dframe.sim_parms$gsi_f <- round(dframe.sim_parms$gsi_f, digits = 3)
  dframe.sim_parms$CP <- round(dframe.sim_parms$CP, digits = 1)
  dframe.sim_parms$ACT <- round(dframe.sim_parms$ACT, digits = 1)
  dframe.sim_parms$prey_ED <- round(dframe.sim_parms$prey_ED, digits = 0)
  dframe.sim_parms$pred_ED <- round(dframe.sim_parms$pred_ED, digits = 0)
  dframe.sim_parms$C1_ins <- round(dframe.sim_parms$C1_ins, digits = 3)
  dframe.sim_parms$C2_ins <- round(dframe.sim_parms$C2_ins, digits = 3)
  dframe.sim_parms$R1_ins <- round(dframe.sim_parms$R1_ins, digits = 3)
  dframe.sim_parms$R2_ins <- round(dframe.sim_parms$R2_ins, digits = 3)
  dframe.sim_parms$F1_ins <- round(dframe.sim_parms$F1_ins, digits = 3)
  dframe.sim_parms$F2_ins <- round(dframe.sim_parms$F2_ins, digits = 3)
  dframe.sim_parms$U1_ins <- round(dframe.sim_parms$U1_ins, digits = 3)
  dframe.sim_parms$U2_ins <- round(dframe.sim_parms$U2_ins, digits = 3)
  dframe.sim_parms$SDA1_ins <- round(dframe.sim_parms$SDA1_ins, digits = 3)
  dframe.sim_parms$SDA2_ins <- round(dframe.sim_parms$SDA2_ins, digits = 3)
  dframe.sim_parms$MS1_ins <- round(dframe.sim_parms$MS1_ins, digits = 3)
  dframe.sim_parms$MS2_ins <- round(dframe.sim_parms$MS2_ins, digits = 3)
  dframe.sim_parms$MG1_ins <- round(dframe.sim_parms$MG1_ins, digits = 3)
  dframe.sim_parms$MG2_ins <- round(dframe.sim_parms$MG2_ins, digits = 3)
  dframe.sim_parms$M1_ins <- round(dframe.sim_parms$M1_ins, digits = 3)
  dframe.sim_parms$M2_ins <- round(dframe.sim_parms$M2_ins, digits = 3)
  dframe.sim_parms$C1_cum <- round(dframe.sim_parms$C1_cum, digits = 3)
  dframe.sim_parms$C2_cum <- round(dframe.sim_parms$C2_cum, digits = 3)
  dframe.sim_parms$R1_cum <- round(dframe.sim_parms$R1_cum, digits = 3)
  dframe.sim_parms$R2_cum <- round(dframe.sim_parms$R2_cum, digits = 3)
  dframe.sim_parms$F1_cum <- round(dframe.sim_parms$F1_cum, digits = 3)
  dframe.sim_parms$F2_cum <- round(dframe.sim_parms$F2_cum, digits = 3)
  dframe.sim_parms$U1_cum <- round(dframe.sim_parms$U1_cum, digits = 3)
  dframe.sim_parms$U2_cum <- round(dframe.sim_parms$U2_cum, digits = 3)
  dframe.sim_parms$SDA1_cum <- round(dframe.sim_parms$SDA1_cum, digits = 3)
  dframe.sim_parms$SDA2_cum <- round(dframe.sim_parms$SDA2_cum, digits = 3)
  dframe.sim_parms$MS1_cum <- round(dframe.sim_parms$MS1_cum, digits = 3)
  dframe.sim_parms$MS2_cum <- round(dframe.sim_parms$MS2_cum, digits = 3)
  dframe.sim_parms$MG1_cum <- round(dframe.sim_parms$MG1_cum, digits = 3)
  dframe.sim_parms$MG2_cum <- round(dframe.sim_parms$MG2_cum, digits = 3)
  dframe.sim_parms$M1_cum <- round(dframe.sim_parms$M1_cum, digits = 3)
  dframe.sim_parms$M2_cum <- round(dframe.sim_parms$M2_cum, digits = 3)
  
  class(dframe.sim_parms) <- c('parms.output', 'data.frame')
  return(dframe.sim_parms)
}

#' Project growth to validate a bioenergetics model with empirical end-of-year mass
#' 
#' @description simulate growth with resampled CP and ACT parameters to identify parameter sets that match empirical end-of-year mass
#' @param start_M2 (default is 100 grams)
#' @param end_M2.empirical mass of fish at end of year based on empirically-derived field or lab size-at-age estimate
#' @param temperature a dataframe populated with a time series of mean daily water temperature (degrees C) of a habitat patch (default is 20 degrees C)
#' @param parms.intrinsic A parms.intrinsic object
#' @param parms.temporal a dataframe populated with a time series of intrinsic (e.g., GSI) and extrinsic (e.g., prey energy density) biological parameters
#' @param C_eq Specify consumption equation 1, 2, or 3 from Hanson et al. 1997
#' @param R_eq Specify respiration equation 1 or 2 from Hanson et al. 1997
#' @param resamp.n specify the number of CP-ACT parameter sets to draw from uniform distributions (default is 1000)
#' @return a list of length 5; raw bem_grow output for each patch and year-end performance indices for each patch; also plotted validation output
#' @export
bem_validate <- function(start_L = 10, end_L.empirical, temperature = 20, parms.intrinsic, parms.temporal, C_eq, R_eq, resamp.n = 1000)
{
  time.start <- Sys.time()
  # resample CP and ACT
  set.seed(1); parms.resample <- data.frame(1:resamp.n, round(runif(resamp.n, 0.0, 1.0), digits = 2), round(runif(resamp.n, 1.0, 4.0), digits = 2))
  colnames(parms.resample) <- c("resampID","CP","ACT")
  
  list.resamp <- list()
  for(i in 1:resamp.n){
    list.resamp[[i]] <- parms.temporal
    list.resamp[[i]]$resampID <- parms.resample[i,"resampID"]
    list.resamp[[i]]$CP <- parms.resample[i,"CP"]
    list.resamp[[i]]$ACT <- parms.resample[i,"ACT"]
  }
  
  # grow fish; loop through n CP-ACT parameter sets
  start_M2 <- parms.intrinsic$lwA * (start_L) ^ parms.intrinsic$lwB
  list.sim_parms <- list()
  for(i in 1:length(list.resamp)){
    list.sim_parms[[i]] <- bem_grow(start_M2, temperature, parms.intrinsic, list.resamp[[i]], C_eq, R_eq)
    list.sim_parms[[i]]$resampID <- parms.resample[i,"resampID"]
    list.sim_parms[[i]]$julian <- 1:365
  }
  names(list.sim_parms) <- parms.resample$resampID
  
  # extract year-end length and mass for n CP-ACT parameter sets
  resampID <- parms.resample$resampID
  CP <- parms.resample$CP
  ACT <- parms.resample$ACT
  M2_cum <- vector()
  L_cum <- vector()
  for(i in 1:resamp.n){
    M2_cum[i] <- round(list.sim_parms[[i]][365,"M2_cum"], digits = 2)
    L_cum[i] <- round(list.sim_parms[[i]][365,"L_cum"], digits = 1)
  }
  dframe.perform <- data.frame(resampID, CP, ACT, M2_cum, L_cum)
  dframe.perform$end_L <- end_L.empirical
  dframe.perform$accurate <- ifelse(dframe.perform$L_cum >= dframe.perform$end_L*0.95 & dframe.perform$L_cum <= dframe.perform$end_L*1.05,
                                    TRUE,
                                    FALSE)
  dframe.perform$accurate <- ifelse(is.na(dframe.perform$accurate), FALSE, dframe.perform$accurate)
  
  df.sim_parms.all <- do.call(rbind, list.sim_parms)
  df.sim_parms.yes <- df.sim_parms.all[df.sim_parms.all$resampID %in% dframe.perform[dframe.perform$accurate == TRUE,"resampID"],]
  
  time.stop <- Sys.time()
  run.time <- time.stop - time.start
  
  run.time <- ifelse(as.numeric(as.duration(run.time)) < 60,
                     noquote(paste("Run time is ", round(as.numeric(as.duration(run.time)), digits = 1), " seconds for ", length(list.sim_parms), " parameter sets", sep = "")),
                     noquote(paste("Run time is ", round(as.numeric(as.duration(run.time) / 60), digits = 2), " minutes for ", length(list.sim_parms), " parameter sets", sep = "")))
  
  # make plot
  plot1 <- ggplot(NULL) +
    theme_classic() +
    labs(x = "Day of year", y = "Total length (cm)") +
    scale_x_continuous(expand = c(0,0), limits = c(-5,367), breaks = seq(0,360,60)) +
    #scale_y_continuous(expand = c(0,0)) +
    geom_line(aes(x = julian, y = L_cum, group = resampID), df.sim_parms.all, alpha = 0.5, color = "gray30") +
    geom_line(aes(x = julian, y = L_cum, group = resampID), df.sim_parms.yes, alpha = 1.0, color = "red") +
    geom_hline(yintercept = end_L.empirical, color = "blue", linetype = "dashed")
  
  # make plot
  plot2 <- ggplot(NULL) +
    theme_classic() +
    labs(title = paste(sum(!is.na(dframe.perform$M2_cum)), " sets completed, ", sum(dframe.perform$accurate == TRUE, na.rm = TRUE), " sets validated", sep = ""),
         x = "Activity multiplier (ACT)", y = "Proportion of consumption (CP)") +
    scale_x_continuous(expand = c(0,0), limits = c(0.96,4.04), breaks = seq(1,4,1)) +
    scale_y_continuous(expand = c(0,0), limits = c(-0.01,1.01), breaks = seq(0,1,0.2)) +
    geom_point(aes(x = ACT, y = CP), dframe.perform, shape = 16, color = "gray") +
    geom_point(aes(x = ACT, y = CP), dframe.perform[dframe.perform$accurate == TRUE,], shape = 16, color = "red")
  
  # combine outputs into a list
  list.results <- list(run.time,
                       list.sim_parms,
                       dframe.perform,
                       plot1,
                       plot2)
  names(list.results) <- c("run.time","daily.output","year.end.perform","plot_growth","plot_parms")
  return(list.results)
}

#' Project growth and performance indices across heterogeneous habitat patches 
#' 
#' @description solve daily energy budgets for n habitat patches
#' @param start_M2 (default is 100 grams)
#' @param temperature a dataframe with 365 rows; the first column is the date; subsequent columns are WT_mean for n habitat patches OR a raster layer with x, y, and z values representing latitudes, longitudes, and 365 days of temperature
#' @param temperature a list of dataframes populated with a time series of mean daily water temperature (degrees C); each list element is a different habitat patch (default is 20 degrees C)
#' @param parms.intrinsic A parms.intrinsic object
#' @param parms.temporal a list of dataframes populated with a time series of intrinsic (e.g., GSI) and extrinsic (e.g., prey energy density) biological parameters; each list element is a different habitat patch
#' @param C_eq Specify consumption equation 1, 2, or 3 from Hanson et al. 1997
#' @param R_eq Specify respiration equation 1 or 2 from Hanson et al. 1997
#' @param is.raster FALSE if temperature data are in tabular format; TRUE if temperature data are in raster format. Default is FALSE
#' @return a list of length 2; raw bem_grow output for each patch and year-end performance indices for each patch
#' @export
bem_project <- function(start_M2 = 100, temperature, parms.intrinsic, parms.temporal, C_eq, R_eq, is.raster = FALSE)
{
  if(isFALSE(is.raster)){
    time.start <- Sys.time()
    
    # reformat temps to list
    list.results <- list()
    for(i in 1:nrow(temperature)){
      list.results[[i]] <- data.frame(colnames(temperature)[-1], as.numeric(temperature[i,-1]))
      colnames(list.results[[i]]) <- c("date","WT_mean")
    }
    names(list.results) <- temperature[,1]
    temperature <- list.results
    
    # grow fish; loop through n habitat patches
    list.sim_parms <- list()
    for(i in 1:length(temperature)){
      list.sim_parms[[i]] <- bem_grow(start_M2, temperature[[i]], parms.intrinsic, parms.temporal, C_eq, R_eq)
    }
    names(list.sim_parms) <- names(temperature)
    
    # extract year-end performance for n habitat patches
    patchID <- names(temperature)
    M2_cum <- vector()
    for(i in 1:length(temperature)){
      M2_cum[i] <- list.sim_parms[[i]][365,"M2_cum"]
    }
    dframe.perform <- data.frame(patchID, M2_cum)
    
    time.stop <- Sys.time()
    run.time <- time.stop - time.start
    run.time <- ifelse(as.numeric(as.duration(run.time)) < 60,
                       noquote(paste("Run time is ", round(as.numeric(as.duration(run.time)), digits = 1), " seconds for ", length(list.sim_parms), " habitat patches", sep = "")),
                       noquote(paste("Run time is ", round(as.numeric(as.duration(run.time) / 60), digits = 2), " minutes for ", length(list.sim_parms), " habitat patches", sep = "")))
    
    # combine outputs into a list
    list.results <- list(run.time,
                         list.sim_parms,
                         dframe.perform)
    names(list.results) <- c("run.time","daily.output","year.end.perform")
    return(list.results)
  }
  else{
    #################################################
    time.start <- Sys.time()
    
    # reformat temps to list
    df.temperature <- data.frame(as.points(temperature))
    colnames(df.temperature) <- paste("julian", str_pad(1:365, 3, pad = 0), sep = "")
    df.temperature <- round(df.temperature, digits = 0)/10
    gridID <- paste("cell", str_pad(1:nrow(df.temperature), nchar(nrow(df.temperature)), pad = 0), sep = "")
    df.temperature <- data.frame(gridID, crds(as.points(temperature)), df.temperature)
    
    df.latlon <- df.temperature[,names(df.temperature) %in% c("gridID","x","y")]
    temperature <- df.temperature[,!names(df.temperature) %in% c("x","y")]
    
    list.results <- list()
    for(i in 1:nrow(temperature)){
      list.results[[i]] <- data.frame(colnames(temperature)[-1], as.numeric(temperature[i,-1]))
      colnames(list.results[[i]]) <- c("date","WT_mean")
    }
    names(list.results) <- temperature[,1]
    temperature <- list.results
    
    # grow fish; loop through n habitat patches
    list.sim_parms <- list()
    for(i in 1:length(temperature)){
      list.sim_parms[[i]] <- bem_grow(start_M2, temperature[[i]], parms.intrinsic, parms.temporal, C_eq, R_eq)
    }
    names(list.sim_parms) <- names(temperature)
    
    # extract year-end performance for n habitat patches
    patchID <- names(temperature)
    M2_cum <- vector()
    for(i in 1:length(temperature)){
      M2_cum[i] <- list.sim_parms[[i]][365,"M2_cum"]
    }
    dframe.perform <- data.frame(patchID, M2_cum)
    
    tmp <- as.data.frame(temps_bayc, xy = TRUE)[,c("x","y")]
    dframe.perform <- cbind(dframe.perform, tmp)
    colnames(dframe.perform) <- c("gridID", "M2_cum", "lon", "lat")
    dframe.perform <- dframe.perform[,c("lon", "lat", "M2_cum")]
    rast.perform <- rast(dframe.perform, crs = "local")
    
    time.stop <- Sys.time()
    run.time <- time.stop - time.start
    run.time <- ifelse(as.numeric(as.duration(run.time)) < 60,
                       noquote(paste("Run time is ", round(as.numeric(as.duration(run.time)), digits = 1), " seconds for ", length(list.sim_parms), " habitat patches", sep = "")),
                       noquote(paste("Run time is ", round(as.numeric(as.duration(run.time) / 60), digits = 2), " minutes for ", length(list.sim_parms), " habitat patches", sep = "")))
    
    # combine outputs into a list
    list.results <- list(run.time,
                         list.sim_parms,
                         rast.perform)
    names(list.results) <- c("run.time","daily.output","year.end.perform")
    return(list.results)
    #################################################
  }
}
