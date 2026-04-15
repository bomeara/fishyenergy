
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
#'   \item{CEQ}{(integer) consumption equation 1, 2, 3, or 4}
#'   \item{CA}{(numeric) intercept estimate for mass-dependent consumption}
#'   \item{CB}{(numeric) slope estimate for mass-dependent consumption}
#'   \item{CQ}{(numeric) Q10 for temperature-dependent consumption}
#'   \item{CTO}{(numeric) thermal optimum for consumption}
#'   \item{CTM}{(numeric) thermal maximum for consumption}
#'   \item{CTL}{(numeric) temperature-dependent parameter for consumption}
#'   \item{CK1}{(numeric) temperature-dependent parameter for consumption}
#'   \item{CK4}{(numeric) temperature-dependent parameter for consumption}
#'   \item{REQ}{(integer) respiration equation 1, 2, or 4}
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
consumption1 <- function(M, T, CP = 1.0, parms.intrinsic)
{
  fT_C <- exp(parms.intrinsic$CQ*round(T, digits = 1))
  Cmax <- parms.intrinsic$CA*M^parms.intrinsic$CB
  C <- Cmax*CP*fT_C
  return(C)
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
  V <- (parms.intrinsic$CTM-round(T, digits = 1))/(parms.intrinsic$CTM-parms.intrinsic$CTO)
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
consumption3 <- function(M, T, CP = 1.0, parms.intrinsic)
{
  G2 <- (1/(parms.intrinsic$CTL-parms.intrinsic$CTM))*log((0.98*(1-parms.intrinsic$CK4))/(parms.intrinsic$CK4*0.02))
  L2 <- exp(G2*(parms.intrinsic$CTL-round(T, digits = 1)))
  KB <- (parms.intrinsic$CK4*L2)/(1+parms.intrinsic$CK4*(L2-1))
  G1 <- (1/(parms.intrinsic$CTO-parms.intrinsic$CQ))* log((0.98*(1-parms.intrinsic$CK1))/(parms.intrinsic$CK1*0.02))
  L1 <- exp(G1)*(round(T, digits = 1)-parms.intrinsic$CQ)
  KA <- (parms.intrinsic$CK1*L1)/(1+parms.intrinsic$CK1*(L1-1))
  fT_C <- KA*KB
  Cmax <- parms.intrinsic$CA*M^parms.intrinsic$CB
  C <- Cmax*CP*fT_C
  return(C)
}

#' Consumption equation 4
#'
#' @description Consumption equation 4; new equation for FishyEnergy that simply matches Cmax for a given temperature
#' @param M Mass (grams) of fish
#' @param T Temperature (degrees C) at which consumption is calculated
#' @param CP Proportion of maximum consumption (default is 1.0)
#' @param parms.intrinsic A parms.intrinsic object; note that temperature dependent parameters are bypassed
#' @param match.table temperature dependent parameters formatted as a table with three columns named WT_mean, C_match, and R_match
#' @return Specific consumption rate (grams of food consumed per gram of fish mass per day)
#' @export
consumption4 <- function(M, T, CP = 1.0, parms.intrinsic, match.table)
{
  fT_C <- match.table[match.table$WT_mean %in% round(T, digits = 1),"C_match"]
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
respiration1 <- function(M, T, ACT = 1.0, parms.intrinsic)
{
  fT_R <- exp(parms.intrinsic$RQ*round(T, digits = 1))
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
respiration2 <- function(M, T, ACT = 1.0, parms.intrinsic)
{
  Y <- log(parms.intrinsic$RQ)*(parms.intrinsic$RTM-parms.intrinsic$RTO+2)
  Z <- log(parms.intrinsic$RQ)*(parms.intrinsic$RTM-parms.intrinsic$RTO)
  X <- (Z^2*(1+(1+40/Y)^0.5)^2)/400
  V <- (parms.intrinsic$RTM-round(T, digits = 1))/(parms.intrinsic$RTM-parms.intrinsic$RTO)
  fT_R <- V^X*exp(X*(1-V))
  Rrest <- parms.intrinsic$RA*M^parms.intrinsic$RB*fT_R
  R <- Rrest*ACT
  return(R)
}

#' Respiration equation 4
#'
#' @description Respiration equation 4; new equation for FishyEnergy that simply matches Rrest for a given temperature
#' @param M Mass (grams) of fish
#' @param T Temperature (degrees C) at which consumption is calculated
#' @param ACT Activity multiplier (default is 1.0)
#' @param parms.intrinsic A parms.intrinsic object; note that temperature dependent parameters are bypassed
#' @param match.table temperature dependent parameters formatted as a table with three columns named WT_mean, C_match, and R_match
#' @return Specific respiration rate (grams of oxygen per gram of fish mass per day)
#' @export
respiration4 <- function(M, T, ACT = 1.0, parms.intrinsic, match.table)
{
  fT_R <- match.table[match.table$WT_mean %in% round(T, digits = 1),"R_match"]
  Rrest <- parms.intrinsic$RA*M^parms.intrinsic$RB
  R <- Rrest*ACT*fT_R
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
#' @param C_eq Specify consumption equation 1, 2, or 3 from Hanson et al. 1997 or new equation 4
#' @param R_eq Specify respiration equation 1 or 2 from Hanson et al. 1997 or new equation 4
#' @return A plot of temperature dependent rates
#' @export
bem_curve <- function(M, T, CP = 1.0, ACT = 1.0, prey_ED = 4000, parms.intrinsic, C_eq, R_eq, match.table)
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
  else if (C_eq == 4){
    C1_curve <- consumption4(M, T, CP, parms.intrinsic, match.table)
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
  else if(R_eq == 4){
    R1_curve <- respiration4(M, T, ACT, parms.intrinsic, match.table)
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
  ggplot2::ggplot(NULL) +
    ggplot2::theme(plot.tag.position = c(0.15, 0.02)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Temperature (°C)",
                  y = "Energy (joules/g/d)") +
    ggplot2::scale_linewidth_manual(values = c("Surplus" = 2,
                                               "Consumption" = 1,
                                               "Respiration" = 1,
                                               "Excretion" = 1,
                                               "Egestion" = 1,
                                               "SDA" = 1)) +
    ggplot2::scale_color_manual(values = c("Surplus" = "black",
                                           "Consumption" = "blue",
                                           "Respiration" = "red",
                                           "Excretion" = "gray65",
                                           "Egestion" = "gray50",
                                           "SDA" = "gray35")) +
    ggplot2::scale_y_continuous(expand = c(0,0),
                                limits = c(0, max(dframe.curves[dframe.curves$units == 1,"value"])*1.05)) +
    ggplot2::geom_line(ggplot2::aes(x = temp, y = value, color = Curve, linewidth = Curve), dframe.curves[dframe.curves$units == 1,])
}

#' Solve daily energy budgets to simulate growth over a calendar year
#'
#' @description the Wisconsin bioenergetics model from Hanson et al. 1997
#' @param start_M2 (default is 100 grams)
#' @param temperature a dataframe populated with a time series of mean daily water temperature (degrees C x 10) of a habitat patch
#' @param parms.intrinsic A parms.intrinsic object; note that temperature dependent parameters are bypassed if C_eq = 4 and/or R_eq = 4
#' @param parms.temporal a dataframe populated with a time series of intrinsic (e.g., GSI) and extrinsic (e.g., prey energy density) biological parameters
#' @param C_eq Specify consumption equation 1, 2, or 3 from Hanson et al. 1997 or equation 4
#' @param R_eq Specify respiration equation 1 or 2 from Hanson et al. 1997
#' @param match.table temperature dependent parameters formatted as a table with three columns named WT_mean, C_match, and R_match; only applicable if C_eq = 4 and/or R_eq = 4
#' @return dataframe populated with output parameters (columns) for time series of projected days (rows)
#' @export
bem_grow <- function(start_M2 = 100, temperature, parms.intrinsic, parms.temporal, C_eq, R_eq, match.table)
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
  dframe.sim_parms[1,c("K")] = 1.0

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}
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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

  }
  else if(C_eq == 4 & R_eq == 1){
    for(i in 2:nrow(dframe.sim_parms)){
      dframe.sim_parms[i,c("C1_ins")] = ifelse(is.nan(consumption4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"CP"],
                                                                   parms.intrinsic,
                                                                   match.table)),
                                               0,                                                                                              # ** if temp > CTM, then consumption is zero
                                               consumption4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"CP"],
                                                            parms.intrinsic)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])

      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration1((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],
                                                                   parms.intrinsic,
                                                                   match.table)),
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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

  }
  else if(C_eq == 4 & R_eq == 2){
    for(i in 2:nrow(dframe.sim_parms)){
      dframe.sim_parms[i,c("C1_ins")] = ifelse(is.nan(consumption4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"CP"],
                                                                   parms.intrinsic,
                                                                   match.table)),
                                               0,                                                                                              # ** if temp > CTM, then consumption is zero
                                               consumption4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"CP"],
                                                            parms.intrinsic,
                                                            match.table)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])

      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],
                                                                   parms.intrinsic,
                                                                   match.table)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic,
                                                            match.table)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

  }
  else if(C_eq == 1 & R_eq == 4){
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

      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],
                                                                   parms.intrinsic,
                                                                   match.table)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic,
                                                            match.table)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

  }
  else if(C_eq == 2 & R_eq == 4){
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

      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],
                                                                   parms.intrinsic,
                                                                   match.table)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic,
                                                            match.table)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

  }
  else if(C_eq == 3 & R_eq == 4){
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

      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],
                                                                   parms.intrinsic,
                                                                   match.table)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic,
                                                            match.table)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

  }
  else if(C_eq == 4 & R_eq == 4){
    for(i in 2:nrow(dframe.sim_parms)){
      dframe.sim_parms[i,c("C1_ins")] = ifelse(is.nan(consumption4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"CP"],
                                                                   parms.intrinsic,
                                                                   match.table)),
                                               0,                                                                                              # ** if temp > CTM, then consumption is zero
                                               consumption4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"CP"],
                                                            parms.intrinsic,
                                                            match.table)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])

      dframe.sim_parms[i,c("R1_ins")] = ifelse(is.nan(respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),  # mass on the preceding day (somatic + gonadal)
                                                                   dframe.sim_parms[i,c("WT_mean")],                                           # temperature on the ith julian day
                                                                   parms.temporal[i,"ACT"],
                                                                   parms.intrinsic,
                                                                   match.table)),
                                               NA,                                                                                             # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                               respiration4((dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")]),         # mass on the preceding day (somatic + gonadal)
                                                            dframe.sim_parms[i,c("WT_mean")],                                                  # temperature on the ith julian day
                                                            parms.temporal[i,"ACT"],
                                                            parms.intrinsic,
                                                            match.table)) * (dframe.sim_parms[i-1,c("MS2_cum")] + dframe.sim_parms[i-1,c("MG2_cum")])

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
      dframe.sim_parms[i,c("K")] = 100*(dframe.sim_parms[i,c("MS2_cum")]/(dframe.sim_parms[i,c("L_cum")]^3))}

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
#' @param start_L starting total length in cm
#' @param end_L.empirical total length (cm) of fish at end of year based on empirically-derived field or lab size-at-age estimate
#' @param temperature a dataframe populated with a time series of mean daily water temperature (degrees C x 10) of a habitat patch
#' @param parms.intrinsic A parms.intrinsic object; note that temperature dependent parameters are bypassed if C_eq = 4 and/or R_eq = 4
#' @param parms.temporal a dataframe populated with a time series of intrinsic (e.g., GSI) and extrinsic (e.g., prey energy density) biological parameters
#' @param C_eq Specify consumption equation 1, 2, or 3 from Hanson et al. 1997 or equation 4
#' @param R_eq Specify respiration equation 1 or 2 from Hanson et al. 1997
#' @param match.table temperature dependent parameters formatted as a table with three columns named WT_mean, C_match, and R_match; only applicable if C_eq = 4 and/or R_eq = 4
#' @param resamp.n specify the number of CP-ACT parameter sets to draw from uniform distributions (default is 1000)
#' @return a list of length 5; raw bem_grow output for each patch and year-end performance indices for each patch; also plotted validation output
#' @export
bem_validate <- function(start_L = 10, end_L.empirical, temperature, parms.intrinsic, parms.temporal, C_eq, R_eq, match.table, resamp.n = 1000)
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
    list.sim_parms[[i]] <- bem_grow(start_M2, temperature, parms.intrinsic, list.resamp[[i]], C_eq, R_eq, match.table)
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

  run.time <- ifelse(as.numeric(lubridate::as.duration(run.time)) < 60,
                     noquote(paste("Run time is ", round(as.numeric(lubridate::as.duration(run.time)), digits = 1), " seconds for ", length(list.sim_parms), " parameter sets", sep = "")),
                     noquote(paste("Run time is ", round(as.numeric(lubridate::as.duration(run.time) / 60), digits = 2), " minutes for ", length(list.sim_parms), " parameter sets", sep = "")))

  # make plot
  plot1 <- ggplot2::ggplot(NULL) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Day of year", y = "Total length (cm)") +
    ggplot2::scale_x_continuous(expand = c(0,0), limits = c(-5,367), breaks = seq(0,360,60)) +
    #ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::geom_line(ggplot2::aes(x = julian, y = L_cum, group = resampID), df.sim_parms.all, alpha = 0.5, color = "gray30") +
    ggplot2::geom_line(ggplot2::aes(x = julian, y = L_cum, group = resampID), df.sim_parms.yes, alpha = 1.0, color = "red") +
    ggplot2::geom_hline(yintercept = end_L.empirical, color = "blue", linetype = "dashed")

  # make plot
  plot2 <- ggplot2::ggplot(NULL) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste(sum(!is.na(dframe.perform$M2_cum)), " sets completed, ", sum(dframe.perform$accurate == TRUE, na.rm = TRUE), " sets validated", sep = ""),
                  x = "Activity multiplier (ACT)", y = "Proportion of consumption (CP)") +
    ggplot2::scale_x_continuous(expand = c(0,0), limits = c(0.96,4.04), breaks = seq(1,4,1)) +
    ggplot2::scale_y_continuous(expand = c(0,0), limits = c(-0.01,1.01), breaks = seq(0,1,0.2)) +
    ggplot2::geom_point(ggplot2::aes(x = ACT, y = CP), dframe.perform, shape = 16, color = "gray") +
    ggplot2::geom_point(ggplot2::aes(x = ACT, y = CP), dframe.perform[dframe.perform$accurate == TRUE,], shape = 16, color = "red")

  # combine outputs into a list
  list.results <- list(run.time,
                       list.sim_parms,
                       dframe.perform,
                       plot1,
                       plot2)
  names(list.results) <- c("run.time","daily.output","year.end.perform","plot_growth","plot_parms")
  return(list.results)
}

#' Project performance indices across heterogeneous habitat patches
#'
#' @description solve daily energy budgets for n habitat patches
#' @param start_M2 (default is 100 grams)
#' @param temperature a dataframe with 365 rows; the first column is the date; subsequent columns are WT_mean (degrees C x 10) for n habitat patches OR a raster layer with x, y, and z values representing latitudes, longitudes, and 365 days of temperature
#' @param parms.intrinsic A parms.intrinsic object
#' @param parms.temporal a list of dataframes populated with a time series of intrinsic (e.g., GSI) and extrinsic (e.g., prey energy density) biological parameters; each list element is a different habitat patch
#' @param C_eq Specify consumption equation 1, 2, or 3 from Hanson et al. 1997 or equation 4
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
      list.results[[i]]$WT_mean <- list.results[[i]]$WT_mean / 10
    }
    names(list.results) <- temperature[,1]
    temperature <- list.results

    # grow fish; loop through n habitat patches
    list.sim_parms <- list()
    for(i in 1:length(temperature)){
      list.sim_parms[[i]] <- bem_grow(start_M2, temperature[[i]], parms.intrinsic, parms.temporal, C_eq, R_eq, match.table)
    }
    names(list.sim_parms) <- names(temperature)

    # extract year-end performance for n habitat patches
    patchID <- names(temperature)
    M2_end <- vector()
    M2_dif <- vector()
    M1_m01 <- vector()
    M1_m02 <- vector()
    M1_m03 <- vector()
    M1_m04 <- vector()
    M1_m05 <- vector()
    M1_m06 <- vector()
    M1_m07 <- vector()
    M1_m08 <- vector()
    M1_m09 <- vector()
    M1_m10 <- vector()
    M1_m11 <- vector()
    M1_m12 <- vector()
    M1_cold <- vector()
    M1_warm <- vector()
    def_days <- vector()
    K_min <- vector()
    K_yday <- vector()

    for(i in 1:length(temperature)){
      M2_end[i] <- list.sim_parms[[i]][365,"M2_cum"]
      M2_dif[i] <- ((list.sim_parms[[i]][365,"M2_cum"] - list.sim_parms[[i]][001,"M2_cum"]) / list.sim_parms[[i]][001,"M2_cum"]) * 100
      M1_m01[i] <- sum(list.sim_parms[[i]][001:031,"M1_ins"])
      M1_m02[i] <- sum(list.sim_parms[[i]][032:059,"M1_ins"])
      M1_m03[i] <- sum(list.sim_parms[[i]][060:090,"M1_ins"])
      M1_m04[i] <- sum(list.sim_parms[[i]][091:120,"M1_ins"])
      M1_m05[i] <- sum(list.sim_parms[[i]][121:151,"M1_ins"])
      M1_m06[i] <- sum(list.sim_parms[[i]][152:181,"M1_ins"])
      M1_m07[i] <- sum(list.sim_parms[[i]][182:212,"M1_ins"])
      M1_m08[i] <- sum(list.sim_parms[[i]][213:243,"M1_ins"])
      M1_m09[i] <- sum(list.sim_parms[[i]][244:273,"M1_ins"])
      M1_m10[i] <- sum(list.sim_parms[[i]][274:304,"M1_ins"])
      M1_m11[i] <- sum(list.sim_parms[[i]][305:334,"M1_ins"])
      M1_m12[i] <- sum(list.sim_parms[[i]][335:365,"M1_ins"])
      WT_month <- c(mean(list.sim_parms[[i]][001:031,"WT_mean"]), mean(list.sim_parms[[i]][032:059,"WT_mean"]), mean(list.sim_parms[[i]][060:090,"WT_mean"]), mean(list.sim_parms[[i]][091:120,"WT_mean"]), mean(list.sim_parms[[i]][121:151,"WT_mean"]), mean(list.sim_parms[[i]][152:181,"WT_mean"]), mean(list.sim_parms[[i]][182:212,"WT_mean"]), mean(list.sim_parms[[i]][213:243,"WT_mean"]), mean(list.sim_parms[[i]][244:273,"WT_mean"]), mean(list.sim_parms[[i]][274:304,"WT_mean"]), mean(list.sim_parms[[i]][305:334,"WT_mean"]), mean(list.sim_parms[[i]][335:365,"WT_mean"]))
      m1_month <- c(sum(list.sim_parms[[i]][001:031,"M1_ins"]), sum(list.sim_parms[[i]][032:059,"M1_ins"]), sum(list.sim_parms[[i]][060:090,"M1_ins"]), sum(list.sim_parms[[i]][091:120,"M1_ins"]), sum(list.sim_parms[[i]][121:151,"M1_ins"]), sum(list.sim_parms[[i]][152:181,"M1_ins"]), sum(list.sim_parms[[i]][182:212,"M1_ins"]), sum(list.sim_parms[[i]][213:243,"M1_ins"]), sum(list.sim_parms[[i]][244:273,"M1_ins"]), sum(list.sim_parms[[i]][274:304,"M1_ins"]), sum(list.sim_parms[[i]][305:334,"M1_ins"]), sum(list.sim_parms[[i]][335:365,"M1_ins"]))
      M1_cold[i] <- m1_month[WT_month == min(WT_month)][1]
      M1_warm[i] <- m1_month[WT_month == max(WT_month)][1]
      K_min[i]  <- min(list.sim_parms[[i]]$K, na.rm = TRUE)
      K_yday[i]  <- as.numeric(gsub("julian", "", list.sim_parms[[i]][list.sim_parms[[i]]$K %in% min(list.sim_parms[[i]]$K, na.rm = TRUE),"date"][1]))
      def_days[i] <- ifelse(min(list.sim_parms[[i]]$M1_ins) >= 0, 0, max(rle(list.sim_parms[[i]]$M1_ins < 0)$lengths[rle(list.sim_parms[[i]]$M1_ins < 0)$values]))
    }
    #rm(WT_month)
    #rm(m1_month)
    dframe.perform <- data.frame(patchID,
                                 M2_end, M2_dif,
                                 M1_m01, M1_m02, M1_m03, M1_m04, M1_m05, M1_m06, M1_m07, M1_m08, M1_m09, M1_m10, M1_m11, M1_m12,
                                 M1_cold, M1_warm,
                                 K_min, K_yday,
                                 def_days)

    time.stop <- Sys.time()
    run.time <- time.stop - time.start
    run.time <- ifelse(as.numeric(lubridate::as.duration(run.time)) < 60,
                       noquote(paste("Run time is ", round(as.numeric(lubridate::as.duration(run.time)), digits = 1), " seconds for ", length(list.sim_parms), " habitat patches", sep = "")),
                       noquote(paste("Run time is ", round(as.numeric(lubridate::as.duration(run.time) / 60), digits = 2), " minutes for ", length(list.sim_parms), " habitat patches", sep = "")))

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
    df.temperature <- data.frame(terra::as.points(temperature)) / 10
    colnames(df.temperature) <- paste("julian", stringr::str_pad(1:365, 3, pad = 0), sep = "")
    gridID <- paste("cell", stringr::str_pad(1:nrow(df.temperature), nchar(nrow(df.temperature)), pad = 0), sep = "")
    df.temperature <- data.frame(gridID, terra::crds(terra::as.points(temperature)), df.temperature)

    df.latlon <- df.temperature[,names(df.temperature) %in% c("gridID","x","y")]
    df.temperature <- df.temperature[,!names(df.temperature) %in% c("x","y")]

    list.results <- list()
    for(i in 1:nrow(df.temperature)){
      list.results[[i]] <- data.frame(colnames(df.temperature)[-1], as.numeric(df.temperature[i,-1]))
      colnames(list.results[[i]]) <- c("date","WT_mean")
    }
    names(list.results) <- df.temperature[,1]
    df.temperature <- list.results

    # grow fish; loop through n habitat patches
    list.sim_parms <- list()
    for(i in 1:length(df.temperature)){
      list.sim_parms[[i]] <- bem_grow(start_M2, df.temperature[[i]], parms.intrinsic, parms.temporal, C_eq, R_eq, match.table)
    }
    names(list.sim_parms) <- names(df.temperature)

    # extract year-end performance for n habitat patches
    patchID <- names(df.temperature)
    M2_end <- vector()
    M2_dif <- vector()
    M1_m01 <- vector()
    M1_m02 <- vector()
    M1_m03 <- vector()
    M1_m04 <- vector()
    M1_m05 <- vector()
    M1_m06 <- vector()
    M1_m07 <- vector()
    M1_m08 <- vector()
    M1_m09 <- vector()
    M1_m10 <- vector()
    M1_m11 <- vector()
    M1_m12 <- vector()
    M1_cold <- vector()
    M1_warm <- vector()
    def_days <- vector()
    K_min <- vector()
    K_yday <- vector()

    for(i in 1:length(list.sim_parms)){
      M2_end[i] <- list.sim_parms[[i]][365,"M2_cum"]
      M2_dif[i] <- ((list.sim_parms[[i]][365,"M2_cum"] - list.sim_parms[[i]][001,"M2_cum"]) / list.sim_parms[[i]][001,"M2_cum"]) * 100
      M1_m01[i] <- sum(list.sim_parms[[i]][001:031,"M1_ins"])
      M1_m02[i] <- sum(list.sim_parms[[i]][032:059,"M1_ins"])
      M1_m03[i] <- sum(list.sim_parms[[i]][060:090,"M1_ins"])
      M1_m04[i] <- sum(list.sim_parms[[i]][091:120,"M1_ins"])
      M1_m05[i] <- sum(list.sim_parms[[i]][121:151,"M1_ins"])
      M1_m06[i] <- sum(list.sim_parms[[i]][152:181,"M1_ins"])
      M1_m07[i] <- sum(list.sim_parms[[i]][182:212,"M1_ins"])
      M1_m08[i] <- sum(list.sim_parms[[i]][213:243,"M1_ins"])
      M1_m09[i] <- sum(list.sim_parms[[i]][244:273,"M1_ins"])
      M1_m10[i] <- sum(list.sim_parms[[i]][274:304,"M1_ins"])
      M1_m11[i] <- sum(list.sim_parms[[i]][305:334,"M1_ins"])
      M1_m12[i] <- sum(list.sim_parms[[i]][335:365,"M1_ins"])
      WT_month <- c(mean(list.sim_parms[[i]][001:031,"WT_mean"]), mean(list.sim_parms[[i]][032:059,"WT_mean"]), mean(list.sim_parms[[i]][060:090,"WT_mean"]), mean(list.sim_parms[[i]][091:120,"WT_mean"]), mean(list.sim_parms[[i]][121:151,"WT_mean"]), mean(list.sim_parms[[i]][152:181,"WT_mean"]), mean(list.sim_parms[[i]][182:212,"WT_mean"]), mean(list.sim_parms[[i]][213:243,"WT_mean"]), mean(list.sim_parms[[i]][244:273,"WT_mean"]), mean(list.sim_parms[[i]][274:304,"WT_mean"]), mean(list.sim_parms[[i]][305:334,"WT_mean"]), mean(list.sim_parms[[i]][335:365,"WT_mean"]))
      m1_month <- c(sum(list.sim_parms[[i]][001:031,"M1_ins"]), sum(list.sim_parms[[i]][032:059,"M1_ins"]), sum(list.sim_parms[[i]][060:090,"M1_ins"]), sum(list.sim_parms[[i]][091:120,"M1_ins"]), sum(list.sim_parms[[i]][121:151,"M1_ins"]), sum(list.sim_parms[[i]][152:181,"M1_ins"]), sum(list.sim_parms[[i]][182:212,"M1_ins"]), sum(list.sim_parms[[i]][213:243,"M1_ins"]), sum(list.sim_parms[[i]][244:273,"M1_ins"]), sum(list.sim_parms[[i]][274:304,"M1_ins"]), sum(list.sim_parms[[i]][305:334,"M1_ins"]), sum(list.sim_parms[[i]][335:365,"M1_ins"]))
      M1_cold[i] <- m1_month[WT_month == min(WT_month)][1]
      M1_warm[i] <- m1_month[WT_month == max(WT_month)][1]
      K_min[i]  <- min(list.sim_parms[[i]]$K, na.rm = TRUE)
      K_yday[i]  <- as.numeric(gsub("julian", "", as.numeric(gsub("julian", "", list.sim_parms[[i]][list.sim_parms[[i]]$K %in% min(list.sim_parms[[i]]$K, na.rm = TRUE),"date"][1]))))
      def_days[i] <- ifelse(min(list.sim_parms[[i]]$M1_ins) >= 0, 0, max(rle(list.sim_parms[[i]]$M1_ins < 0)$lengths[rle(list.sim_parms[[i]]$M1_ins < 0)$values]))
    }
    #rm(WT_month)
    #rm(m1_month)
    dframe.perform <- data.frame(patchID,
                                 M2_end, M2_dif,
                                 M1_m01, M1_m02, M1_m03, M1_m04, M1_m05, M1_m06, M1_m07, M1_m08, M1_m09, M1_m10, M1_m11, M1_m12,
                                 M1_cold, M1_warm,
                                 K_min, K_yday,
                                 def_days)

    tmp <- as.data.frame(terra::crds(temperature), xy = TRUE)[,c("x","y")]
    dframe.perform <- cbind(dframe.perform, tmp)
    colnames(dframe.perform) <- c("gridID",
                                  "M2_end", "M2_dif",
                                  "M1_m01", "M1_m02", "M1_m03", "M1_m04", "M1_m05", "M1_m06", "M1_m07", "M1_m08", "M1_m09", "M1_m10", "M1_m11", "M1_m12",
                                  "M1_cold", "M1_warm",
                                  "K_min", "K_yday",
                                  "def_days",
                                  "lon", "lat")

    rast.perform <- c(terra::rast(dframe.perform[, c("lon", "lat", "M2_end")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M2_dif")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m01")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m02")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m03")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m04")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m05")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m06")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m07")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m08")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m09")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m10")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m11")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_m12")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_cold")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "M1_warm")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "K_min")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "K_yday")], crs = "local"),
                      terra::rast(dframe.perform[, c("lon", "lat", "def_days")], crs = "local"))

    time.stop <- Sys.time()
    run.time <- time.stop - time.start
    run.time <- ifelse(as.numeric(lubridate::as.duration(run.time)) < 60,
                       noquote(paste("Run time is ", round(as.numeric(lubridate::as.duration(run.time)), digits = 1), " seconds for ", length(list.sim_parms), " habitat patches", sep = "")),
                       noquote(paste("Run time is ", round(as.numeric(lubridate::as.duration(run.time) / 60), digits = 2), " minutes for ", length(list.sim_parms), " habitat patches", sep = "")))

    # combine outputs into a list
    list.results <- list(run.time,
                         list.sim_parms,
                         rast.perform)
    names(list.results) <- c("run.time","daily.output","year.end.perform")
    return(list.results)
    #################################################
  }
}
