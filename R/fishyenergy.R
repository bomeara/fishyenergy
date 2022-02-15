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
	class(result) <- c('BEM', 'numeric')
	return(result)
}

#' Consumption of energy
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
	C <- Cmax*p*fT_C
	return(C)	
}

#' Respiration 1 equation
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
	return(results)
}