# This function takes persondata 
# persondata has following colums:
# time, n, optional: age, gender, MET, CO2rate in L/s at 1 atm & 0 Celsius
# time must be in HOURS since start of measurement
# temp in Celsius
# frequency is frequency of measurements in seconds

# time needs to be in hours, starting at 0 (i.e., it can represent
# time since start of measurements), and frequency should match the 
# frequency of measurements in the CO2 data

persondata_to_emission <- function(persondata, temp=25, freq){
  
  # external CO2ratedata
  # use CO2 generation rates https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5666301/ Table 4
  CO2ratedata <- data.frame(
    gender = c(rep('M', 13*7), rep('F', 13*7)),
    minage = c(rep(c(0, 1, 3, 6, 11, 16, 21, 30, 40, 50, 60, 70, 80), each=7, times=2)),
    minMET = c(rep(c(1, 1.2, 1.4, 1.6, 2, 3, 4), 26)),
    CO2rate = c(.0009, .0011, .0013, .0014, .0018, 0.0027,	0.0036,
                .0015, .0018, .0021, .0024, .0030, 0.0044,	0.0059,
                0.0019, 0.0023,	0.0026,	0.0030,	0.0038, 0.0057, 0.0075,
                0.0025,	0.0030,	0.0035,	0.0040,	0.0050, 0.0075, 0.0100,
                0.0034,	0.0041,	0.0048,	0.0054,	0.0068, 0.0102,	0.0136,
                0.0037,	0.0045,	0.0053,	0.0060,	0.0075, 0.0113,	0.0150,
                0.0039,	0.0048,	0.0056,	0.0064,	0.0080, 0.0120,	0.0160,
                0.0037, 0.0046, 0.0053, 0.0061, 0.0076, 0.0114,	0.0152,
                0.0038, 0.0046, 0.0054, 0.0062, 0.0077, 0.0116,	0.0155,
                0.0038, 0.0046, 0.0054, 0.0062, 0.0077, 0.0116,	0.0154,
                0.0033, 0.0040, 0.0046, 0.0053, 0.0066, 0.0099,	0.0133,
                0.0031, 0.0038, 0.0045, 0.0051, 0.0064, 0.0095,	0.0127,
                0.0030, 0.0036, 0.0042, 0.0048, 0.0060, 0.0090,	0.0120,
                0.0008, 0.0010, 0.0012, 0.0014, 0.0017, 0.0025,	0.0034,
                0.0014, 0.0017, 0.0020, 0.0022, 0.0028, 0.0042,	0.0056,
                0.0017, 0.0021, 0.0024, 0.0028, 0.0035, 0.0052,	0.0070,
                0.0023, 0.0027, 0.0032, 0.0037, 0.0046, 0.0069,	0.0092,
                0.0029, 0.0035, 0.0041, 0.0047, 0.0058, 0.0088,	0.0117,
                0.0029, 0.0036, 0.0042, 0.0047, 0.0059, 0.0089,	0.0119,
                0.0031, 0.0038, 0.0044, 0.0050, 0.0063, 0.0094,	0.0126,
                0.0029, 0.0035, 0.0041, 0.0047, 0.0059, 0.0088,	0.0118,
                0.0029, 0.0036, 0.0042, 0.0048, 0.0060, 0.0090,	0.0119,
                0.0030, 0.0036, 0.0042, 0.0048, 0.0060, 0.0090,	0.0120,
                0.0027, 0.0033, 0.0038, 0.0044, 0.0055, 0.0082,	0.0110,
                0.0026, 0.0032, 0.0037, 0.0042, 0.0053, 0.0079,	0.0106,
                0.0025, 0.0030, 0.0035, 0.0040, 0.0050, 0.0075,	0.0101)
  )
  
  ## Preparing input data (persondata) ###
  if(!('time' %in% colnames(persondata))) 
    stop('Must include time column in persondata')
  if(!('n' %in% colnames(persondata)))
    stop('Must include n column in persondata')
  if(!('CO2rate' %in% colnames(persondata))) persondata$CO2rate <- NA
  if(!('MET' %in% colnames(persondata))){
    persondata$minMET <- 1.4 # assume sitting
  } else {
    persondata$MET <- as.numeric(persondata$MET)
    # set input data METs to 1, 1.2, 1.4, 1.6, 2, 3 or 4
    persondata$minMET <- cut(persondata$MET, 
                             breaks=c(0, 1, 1.2, 1.4, 1.6, 2, 3, 4, Inf), 
                             labels=c(1, 1, 1.2, 1.4, 1.6, 2, 3, 4), 
                             right=FALSE, 
                             include.lowest=TRUE)
    persondata$minMET[is.na(persondata$minMET)] <- 1.4 # assume sitting
  }
  if(!('gender' %in% colnames(persondata))) persondata$gender <- NA # use M/F average
  if(!('age' %in% colnames(persondata))) {
    persondata$minage <- 30 # assume adult
  } else {
    persondata$age <- as.numeric(persondata$age)
    # set input data ages to 0, 1, 3, 6, 11, 16, 21, 30, 40, 50, 60, 70, 80
    persondata$minage <- cut(persondata$age, 
                             breaks=c(0, 1, 3, 6, 11, 16, 21, 30, 40, 50, 60, 70, 80, Inf), 
                             labels=c(0, 1, 3, 6, 11, 16, 21, 30, 40, 50, 60, 70, 80), 
                             right = FALSE, 
                             include.lowest=TRUE)
    persondata$minage[is.na(persondata$minage)] <- 30 # assume adult
  }
  
  # calculate averages for male/female CO2 rates
  CO2ratedata$avg_gender_CO2rate <- ave(CO2ratedata$CO2rate, CO2ratedata$minage, CO2ratedata$minMET, FUN=mean)
  
  # use CO2 rate data to fill in missing values in persondata$CO2rate
  # use gender if available; otherwise, average male and female CO2 rate
  for(i in 1:nrow(persondata)){
    if(persondata$n[i]>0 & is.na(persondata$CO2rate[i])) {
      if(is.na(persondata$gender[i])) {
        persondata$CO2rate[i] <- CO2ratedata$avg_gender_CO2rate[(CO2ratedata$minage %in% persondata$minage[i]) & 
                                                                 (CO2ratedata$minMET %in% persondata$minMET[i])][1]
      } else 
        persondata$CO2rate[i] <- CO2ratedata$CO2rate[(CO2ratedata$minage %in% persondata$minage[i]) & 
                                                      (CO2ratedata$minMET %in% persondata$minMET[i]) &
                                                      (CO2ratedata$gender %in% persondata$gender[i])]
    } 
  }
  
  # set CO2 rate to 0 whenever there are 0 people 
  persondata$CO2rate[persondata$n==0] <- 0
  
  
  # some unit conversions
  # CO2rate adjusted to room temperature and convert to mg/s
  # Gmole = 0.0446*G (L/s)
  # Gadj (L/s) = 8.314 (J/mol/K) * (273.15+temp) (K) / 101.325 (kPa) * Gmole (mol/s)
  # note: kPa = J/L
  # Gadj (g/s) = 1.965 * Gadj (L/s)
  # putting it all together: convert from L/s to ug/s at a particular temperature
  persondata$adj_CO2rate <- persondata$CO2rate * 1.965 * .0046 * (8.314*(273.15+temp)/101.325) * 1000
  
  # create time and n*CO2rate vectors
  times = seq(min(as.numeric(persondata$time)), 
              max(as.numeric(persondata$time)), 
              by=freq/60/60)
  nadj_CO2rate = rep(0, length(times))
  unique_times = sort(unique(as.numeric(persondata$time)))
  for(i in 1:(length(unique_times)-1)) {
    indtimes <- which(times >= unique_times[i] & times < unique_times[i+1])
    indinput <- which(as.numeric(persondata$time)==unique_times[i])
    nadj_CO2rate[indtimes] <- sum(persondata$n[indinput]*persondata$adj_CO2rate[indinput])
  }
  
  return(list(times=times, nadj_CO2rate=nadj_CO2rate))
}
