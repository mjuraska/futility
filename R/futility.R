#' @import graphics
NULL
#' @import stats
NULL
#' @import utils
NULL

globalVariables(c("trialObj", "pp", "eventTime"))

is.TRUE <- function(x){
  if ( !is.logical(x) ) stop("Argument to 'is.TRUE' must be of type Logical")
  return(x & !is.na(x))
}

is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


## This function takes as input a vector of times ( we are
## using time in weeks since enrollment), and returns the
## time of the next or previous scheduled visit.  The returned times
## are according to the visit schedule

getVisitWeek <- function( week, visitWeeks, whichVisit=c("next","previous")){
  whichVisit <- match.arg( whichVisit )

  ## check for NAs
  noNAs <- !is.na( week )
  wk <- week[ noNAs ]

  ## find the intervals that the values of 'wk' lie in
  interval <- cut( wk, breaks=visitWeeks, right=FALSE, labels=FALSE)

  ## If we want the "previous" visit, we return the lower bound of the interval,
  ## If we want the "next" we return the upper bound of the interval
  if (whichVisit == "previous")
  {
    week[ noNAs ] <- visitWeeks[ interval ]
  } else if (whichVisit == "next")
  {
    week[ noNAs ] <- visitWeeks[ interval + 1 ]
  }

  return(week)
}

FillinInterimdata.Pooled <- function(interimData, rates, visitSchedule, visitSchedule2 = NULL, Nppt, fuTime, ppAnalysis=FALSE, missVaccProb=NULL, ppAtRiskTimePoint=NULL, Seed = NULL){

    ## Finish enrollment and get enrollment times
    Nenroll <- Nppt - NROW(interimData)

    #Keep recruiting week by week untill Nenroll participants have been recruited
    ii<-0
    Nenrolled<-0
    enr<-NULL
    repeat{
      ii<-ii+1
      if ( !is.null(Seed) ) set.seed( Seed+ii )
      N.thisweek<-rpois(1, lambda = rates$enrollmentRate)
      enr<-c(enr, max(interimData$entry)+(ii-1)+sort(runif(N.thisweek, min=0, max=1)))
      Nenrolled<-Nenrolled+N.thisweek
      if ( Nenrolled >= Nenroll ) break
    }

    #restrict to the first nEnroll enrollees
    enr<-enr[1:Nenroll]

    ## generate dropout times for the Nenroll participants from an exponential distribution
    if ( !is.null(Seed) ) set.seed( Seed )
    dropout <- rexp(Nenroll, rates$dropRate)

    ## generate Time-to-Event for the Nenroll participants
    #here we can't use the same Seed as dropout because when rates$dropRate<rates$eventRate,
    #the same seed will always generate dropout>event.
    if ( !is.null(Seed) ) set.seed( Seed+1000 )
    event.time <- rexp(Nenroll, rates$eventRate)

    ## censor the Dropout Time  and Time-to-Event based on fuTime
    ##-----------------------------------------------------------
    ## we let the support of dropout being [0, fuTime)
    ## if dropout=0, it means dropout happens right after time=0
    ## if dropout=fuTime, it mean dropout happens right after fuTime, thus won't be observed
    ## similar for event
    ##-----------------------------------------------------------
    dropout[dropout >= fuTime] <- NA
    event.time[event.time >= fuTime] <- NA

    ## Censor events that occur after dropout, as they can't be observed
    event.time[ is.TRUE( event.time > dropout ) ] <- NA


    ## *** NOTE  ****************************************************
    ##  We still have records with both Dropout Time and Time-to-Event for those with dropout > event
    ##  We need to keep both times until we calculate the "DX time" for the events,
    ##  then we can compare dropout to the DX time and keep the earlier of the two.
    ## **************************************************************


    ## ---------------------------------------------------------------------- ##
    ## Now, we start to construct the "observed" data, by applying the        ##
    ## study visit Map (Schedule 1 and Schedule 4) to the simulated data.     ##
    ## This will allow us to figure out when: (a) the events are diagnosed##
    ## (b) the last visit occurred (during the simulated period).             ##
    ## ---------------------------------------------------------------------- ##

    ## create observed data object to fill in
    obsEDI <- data.frame( enrollTime = enr,
                          dropTime = dropout,
                          eventDXTime = NA
    )

    ## get event diagnosis dates for each non-NA Time-to-Event

    ## Do the following only if there's at least one event
    if ( !all( is.na(event.time) ) )
    {
      #for simulated newly enrolled subjects, everyone is in visitShcedule1
      eventDX <- getVisitWeek( week=event.time, visitWeeks=visitSchedule, whichVisit = "next")

      ## Compute the minimum of the eventDX time, dropout time
      ## Any events occurring *strictly* after this time are censored
      minTime <- pmin( eventDX, dropout, na.rm=TRUE)

      ## censor eventDX and dropout times
      eventDX[ is.TRUE(eventDX > minTime) ] <- NA
      dropout[ is.TRUE(dropout > minTime) ] <- NA

      ## store info in obsEDI
      obsEDI$eventDXTime <- eventDX
      obsEDI$dropTime <- dropout
    }

    ## we compute the amount of "follow-up time"  for each
    ## ppt.  This is equal to the 'fuTime' for ppt.s without events,
    ## and equal to the event time for participants with events.
    obsEDI$futime <- pmin( obsEDI$eventDXTime,
                           obsEDI$dropTime, fuTime, na.rm=TRUE)

    exit  <- obsEDI$enrollTime + obsEDI$futime

    ## create flag for event
    event <- !is.na(obsEDI$eventDXTime)
    droppedout<- !is.na(obsEDI$dropTime)

    out <- data.frame(
      entry = enr,
      exit  = exit,
      event = as.integer(event),
      dropout=as.integer(droppedout)
    )

    # generate the indicator of membership in the per-protocol cohort
    if (ppAnalysis){
      if (!is.null(Seed)){ set.seed(Seed) }
      out$missVacc <- rbinom(NROW(out), 1, prob=missVaccProb)
      out$pp <- as.numeric(out$missVacc==0 & out$exit - out$entry > ppAtRiskTimePoint)
    }

    ##Now we are done filling in data that are newly enrolled, which is "out"
    ##Next we need to fill in data for the already enrolled in the interimData

    ##the subjects that need to be filled in are those still being followed-up,
    ##i.e. followup=1, or exit=missing

    #we need to generate a dropout time and a Time-to-Event
    #because exponetial distribution is memoryless, the times are generated with same weekly rates
    Nfollowup<-sum(is.na(interimData$exit))

    tFU<-ifelse(interimData$followup==1,interimData$last_visit_dt-interimData$entry,interimData$exit-interimData$entry)
    #followup time for subjects in follow-up
    adj.fUP <- tFU[interimData$followup==1]
    fuTime2 <- fuTime-adj.fUP

    ####################################
    #Starting here, we changed reference time (time=0) to the end of observed follow-up time to each subject in follow-up.
    ####################################
    ## similar process as before
    if ( !is.null(Seed) ) set.seed( Seed+2000 )
    dropout <- rexp(Nfollowup, rate=rates$dropRate)
    if ( !is.null(Seed) ) set.seed( Seed+3000 )
    event.time <- rexp(Nfollowup, rate=rates$eventRate)

    ## fuTime needs to be adjusted to the amount of time the patients were already in the study

    ## censor the dropout time and Time-to-Event based on fuTime
    dropout[dropout >= fuTime2] <- NA
    event.time[event.time >= fuTime2] <- NA

    ## Censor events that occur after dropout, as they can't be observed
    event.time[ is.TRUE( event.time > dropout ) ] <- NA

    ## create observed data object to fill in
    if(is.null(visitSchedule2)){    obsEDI <- data.frame( enrollTime = 0, #interimData$entry[is.na(interimData$exit)],
                                                          schedule2 = FALSE,
                                                          dropTime = dropout,
                                                          eventDXTime = NA
    )}else{
      if(!"schedule2" %in% colnames(interimData)){
        stop("interimData does not have a variable named 'schedule2'")
      }else{
        obsEDI <- data.frame( enrollTime = 0,
                              schedule2 = interimData$schedule2[is.na(interimData$exit)],
                              dropTime = dropout,
                              eventDXTime = NA
        )
      }
    }

    ## get event diagnosis dates for each non-NA Time-to-Event
    event.time1<-event.time[obsEDI$schedule2==FALSE]
    event.time2<-event.time[obsEDI$schedule2==TRUE]
    dropout1<-dropout[obsEDI$schedule2==FALSE]
    dropout2<-dropout[obsEDI$schedule2==TRUE]

    if ( !all( is.na(event.time1) ) )
    {
      eventDX1 <- getVisitWeek( week=event.time1, visitWeeks=visitSchedule, whichVisit = "next")

      ## Compute the minimum of the eventDX time, dropout time
      ## Any events occurring *strictly* after this time are censored
      minTime1 <- pmin( eventDX1, dropout1, na.rm=TRUE)

      ## censor eventDX and dropout times
      eventDX1[ is.TRUE(eventDX1 > minTime1) ] <- NA
      dropout1[ is.TRUE(dropout1 > minTime1) ] <- NA

      ## store info in obsEDI
      obsEDI$eventDXTime[obsEDI$schedule2==FALSE] <- eventDX1
      obsEDI$dropTime[obsEDI$schedule2==FALSE]  <- dropout1
    }

    if ( !all( is.na(event.time2) ) )
    {
      eventDX2 <- getVisitWeek( week=event.time2, visitWeeks=visitSchedule2, whichVisit = "next")

      ## Compute the minimum of the eventDX time, dropout time
      ## Any events occurring *strictly* after this time are censored
      minTime2 <- pmin( eventDX2, dropout2, na.rm=TRUE)

      ## censor eventDX and dropout times
      eventDX2[ is.TRUE(eventDX2 > minTime2) ] <- NA
      dropout2[ is.TRUE(dropout2 > minTime2) ] <- NA

      ## store info in obsEDI
      obsEDI$eventDXTime[obsEDI$schedule2==TRUE] <- eventDX2
      obsEDI$dropTime[obsEDI$schedule==TRUE]  <- dropout2
    }

    obsEDI$futime <- pmin( obsEDI$eventDXTime,
                           obsEDI$dropTime, fuTime2, na.rm=TRUE)

    exit  <- obsEDI$enrollTime + obsEDI$futime

    ## create flag for event
    event <- !is.na(obsEDI$eventDXTime)
    droppedout<- !is.na(obsEDI$dropTime)


    interimData.filled <- data.frame(
      entry = interimData$entry,
      exit  = interimData$exit,
      event = as.integer(interimData$event),
      dropout=as.integer(interimData$dropout)
    )

    if (ppAnalysis){
      interimData.filled$missVacc <- interimData$missVacc
      interimData.filled$pp <- interimData$pp
    }

    ####################################
    #Change exit from reference time being the end of observed follow-up time to each subject in follow-up back to orignial timeline.
    ####################################
    #interimData.filled$exit[is.na(interimData$exit)]<-exit
    interimData.filled$exit[is.na(interimData$exit)]<-interimData.filled$entry[is.na(interimData$exit)]+adj.fUP+exit

    #sanity check
    summary(interimData.filled$exit[is.na(interimData$exit)]-interimData.filled$entry[is.na(interimData$exit)])

    interimData.filled$event[is.na(interimData$exit)]<-event
    interimData.filled$dropout[is.na(interimData$exit)]<-droppedout

    # generate the indicator of membership in the per-protocol cohort
    if (ppAnalysis){
      if (!is.null(Seed)){ set.seed(Seed+10000) }
      interimData.filled$missVacc <- ifelse(interimData$missVacc==0 & is.na(interimData$pp), rbinom(NROW(interimData), 1, prob=missVaccProb), interimData$missVacc)
      interimData.filled$pp <- ifelse(is.na(interimData$pp), as.numeric(interimData$missVacc==0 & interimData$exit - interimData$entry > ppAtRiskTimePoint), interimData$pp)
    }

    out <- rbind(interimData.filled, out)
    out$missVacc <- NULL

    return(out)
  }


FillinInterimdata.byArm <- function(interimData, rates, visitSchedule, visitSchedule2 = NULL, trtNames, N, fuTime, ppAnalysis=FALSE, missVaccProb=NULL, ppAtRiskTimePoint=NULL, Seed = NULL){
    nArms <- length(N)
    ## Finish enroll by arm and get enrollment times
    Nenroll<-NULL
    for(j in 1:nArms){
      Nenroll<-c(Nenroll,N[j]-sum(interimData$arm==trtNames[j]) )
    }

    #Keep recruiting week by week untill Nenroll participants have been recruited, do this for each arm
    enr<-NULL
    arm<-NULL

    for(j in 1:nArms){
      ii<-0
      Nenrolled<-0
      repeat{
        ii<-ii+1
        if ( !is.null(Seed) ) set.seed( Seed+ii )
        N.thisweek<-rpois(1, lambda = rates$enrollmentRate)
        enr<-c(enr, max(interimData$entry)+(ii-1)+sort(runif(N.thisweek, min=0, max=1)))
        arm<-c(arm, rep(trtNames[j],N.thisweek))
        Nenrolled<-Nenrolled+N.thisweek
        if ( Nenrolled >= Nenroll[j] ) break
      }

      #restrict to the first nEnroll enrollees
      Nenroll2 <-sum(Nenroll[1:j])
      enr<-enr[1:Nenroll2]
      arm<-arm[1:Nenroll2]
    }

    ## generate dropout times for the Nenroll participants from an exponential distribution
    if ( !is.null(Seed) ) set.seed( Seed )
    dropout <- rexp(length(enr), rates$dropRate)

    ## generate Time-to-Event by treatment arm
    #here we can't use the same Seed as dropout because when rates$dropRate<rates$eventRate,
    #the same seed will always generate dropout>event.
    event.time <- NULL
    for(j in 1:nArms){
      if ( !is.null(Seed) ) set.seed( Seed+1000*j )
      event.time<-c(event.time,  rexp(Nenroll[j], rates$eventRate[j]))
    }


    ## censor the dropout time and Time-to-Event based on fuTime
    ##-----------------------------------------------------------
    ## we let the support of dropout being [0, fuTime)
    ## if dropout=0, it means dropout happens right after time=0
    ## if dropout=fuTime, it mean dropout happens right after fuTime, thus won't be observed
    ## similar for event
    ##----------------------------------------------------------
    dropout[dropout >= fuTime] <- NA
    event.time[event.time >= fuTime] <- NA

    ## Censor events that occur after dropout, as they can't be observed
    event.time[ is.TRUE( event.time > dropout ) ] <- NA


    ## *** NOTE  ****************************************************
    ##  We still have records with both dropout time and Time-to-Event for those with dropout > event
    ##  We need to keep both times until we calculate the "DX time" for the evnets,
    ##  then we can compare dropout to the DX time and keep the earlier of the two.
    ## **************************************************************


    ## ---------------------------------------------------------------------- ##
    ## Now, we start to construct the "observed" data, by applying the        ##
    ## study visit Map (Schedule 1 and Schedule 4) to the simulated data.     ##
    ## This will allow us to figure out when: (a) the event are diagnosed##
    ## (b) the last visit occurred (during the simulated period).             ##
    ## ---------------------------------------------------------------------- ##

    ## create observed data object to fill in
    obsEDI <- data.frame( arm = arm,
                          enrollTime = enr,
                          dropTime = dropout,
                          eventDXTime = NA
    )

    ## get event diagnosis dates for each non-NA event time

    ## Do the following only if there's at least one event
    if ( !all( is.na(event.time) ) )
    {
      #for simulated newly enrolled subjects, everyone is in visitShcedule1
      eventDX <- getVisitWeek( week=event.time, visitWeeks=visitSchedule, whichVisit = "next")

      ## Compute the minimum of the eventDX time, dropout time
      ## Any events occurring *strictly* after this time are censored
      minTime <- pmin( eventDX, dropout, na.rm=TRUE)

      ## censor eventDX and dropout times
      eventDX[ is.TRUE(eventDX > minTime) ] <- NA
      dropout[ is.TRUE(dropout > minTime) ] <- NA

      ## store info in obsEDI
      obsEDI$eventDXTime <- eventDX
      obsEDI$dropTime <- dropout
    }

    ## we compute the amount of "follow-up time"  for each
    ## ppt.  This is equal to the 'fuTime' for ppt.s without events,
    ## and equal to the event time for participants with events.
    obsEDI$futime <- pmin( obsEDI$eventDXTime,
                           obsEDI$dropTime, fuTime, na.rm=TRUE)

    exit  <- obsEDI$enrollTime + obsEDI$futime

    ## create flag for event
    event <- !is.na(obsEDI$eventDXTime)
    droppedout<- !is.na(obsEDI$dropTime)

    out <- data.frame(
      arm = arm,
      entry = enr,
      exit  = exit,
      event = as.integer(event),
      dropout=as.integer(droppedout)
    )

    # generate the indicator of membership in the per-protocol cohort
    if (ppAnalysis){
      if (!is.null(Seed)){ set.seed(Seed) }
      out$missVacc <- rbinom(NROW(out), 1, prob=missVaccProb)
      out$pp <- as.numeric(out$missVacc==0 & out$exit - out$entry > ppAtRiskTimePoint)
    }

    ##Now we are done filling in data that are newly enrolled, which is "out"
    ##Next we need to fill in data for the already enrolled in the interimData

    ##the subjects that need to be filled in are those still being followed-up,
    ##i.e. followup=1, or exit=missing

    #we need to generate a dropout time and a Time-to-Event
    #because exponetial distribution is memoryless, the times are generated with same weekly rates
    arm<-interimData$arm[is.na(interimData$exit)]
    Nfollowup<-NULL
    for(j in 1:nArms){
      Nfollowup<-c(Nfollowup,sum(is.na(interimData$exit[interimData$arm==trtNames[j]])))
    }

    tFU<-ifelse(interimData$followup==1,interimData$last_visit_dt-interimData$entry,interimData$exit-interimData$entry)
    #followup time for subjects in follow-up
    adj.fUP <- tFU[interimData$followup==1]
    fuTime2 <- fuTime-adj.fUP

    ####################################
    #Starting here, we changed reference time (time=0) to the end of observed follow-up time to each subject in follow-up.
    ####################################
    ## similar process as before
    if ( !is.null(Seed) ) set.seed( Seed+4000 )
    dropout <- rexp(sum(Nfollowup), rate=rates$dropRate)

    event.time <- rep(NA, sum(Nfollowup))
    for(j in 1:nArms){
      if ( !is.null(Seed) ) set.seed( Seed+4000+1000*j )
      event.time[arm==trtNames[j]] <- rexp(Nfollowup[j], rates$eventRate[j] )
    }

    ## fuTime needs to be adjusted to the amount of time the patients were already in the study

    ## censor the dropout time and Time-to-Event based on fuTime
    dropout[dropout >= fuTime2] <- NA
    event.time[event.time >= fuTime2] <- NA

    ## Censor events that occur after dropout, as they can't be observed
    event.time[ is.TRUE( event.time > dropout ) ] <- NA

    ## create observed data object to fill in
    if(is.null(visitSchedule2)){     obsEDI <- data.frame( arm = arm,
                                                           enrollTime = 0, #interimData$entry[is.na(interimData$exit)],
                                                           schedule2 = FALSE,
                                                           dropTime = dropout,
                                                           eventDXTime = NA
    )}else{
      if(!"schedule2" %in% colnames(interimData)){
        stop("interimData does not have a variable named 'schedule2'")
      }else{
        obsEDI <- data.frame( arm = arm,
                              enrollTime = 0, #interimData$entry[is.na(interimData$exit)],
                              schedule2 = interimData$schedule2[is.na(interimData$exit)],
                              dropTime = dropout,
                              eventDXTime = NA
        )
      }
    }


    ## get event diagnosis dates for each non-NA Time-to-Event
    event.time1<-event.time[obsEDI$schedule2==FALSE]
    event.time2<-event.time[obsEDI$schedule2==TRUE]
    dropout1<-dropout[obsEDI$schedule2==FALSE]
    dropout2<-dropout[obsEDI$schedule2==TRUE]

    if ( !all( is.na(event.time1) ) )
    {
      eventDX1 <- getVisitWeek( week=event.time1, visitWeeks=visitSchedule, whichVisit = "next")

      ## Compute the minimum of the eventDX time, dropout time
      ## Any events occurring *strictly* after this time are censored
      minTime1 <- pmin( eventDX1, dropout1, na.rm=TRUE)

      ## censor eventDX and dropout times
      eventDX1[ is.TRUE(eventDX1 > minTime1) ] <- NA
      dropout1[ is.TRUE(dropout1 > minTime1) ] <- NA

      ## store info in obsEDI
      obsEDI$eventDXTime[obsEDI$schedule2==FALSE] <- eventDX1
      obsEDI$dropTime[obsEDI$schedule2==FALSE]  <- dropout1
    }

    if ( !all( is.na(event.time2) ) )
    {
      eventDX2 <- getVisitWeek( week=event.time2, visitWeeks=visitSchedule2, whichVisit = "next")

      ## Compute the minimum of the eventDX time, dropout time
      ## Any events occurring *strictly* after this time are censored
      minTime2 <- pmin( eventDX2, dropout2, na.rm=TRUE)

      ## censor eventDX and dropout times
      eventDX2[ is.TRUE(eventDX2 > minTime2) ] <- NA
      dropout2[ is.TRUE(dropout2 > minTime2) ] <- NA

      ## store info in obsEDI
      obsEDI$eventDXTime[obsEDI$schedule2==TRUE] <- eventDX2
      obsEDI$dropTime[obsEDI$schedule2==TRUE]  <- dropout2
    }

    obsEDI$futime <- pmin( obsEDI$eventDXTime,
                           obsEDI$dropTime, fuTime2, na.rm=TRUE)

    exit  <- obsEDI$enrollTime + obsEDI$futime

    ## create flag for event
    event <- !is.na(obsEDI$eventDXTime)
    droppedout<- !is.na(obsEDI$dropTime)


    interimData.filled <- data.frame(
      arm = interimData$arm,
      entry = interimData$entry,
      exit  = interimData$exit,
      event = as.integer(interimData$event),
      dropout=as.integer(interimData$dropout)
    )

    if (ppAnalysis){
      interimData.filled$missVacc <- interimData$missVacc
      interimData.filled$pp <- interimData$pp
    }

    ####################################
    #Change exit from reference time being the end of observed follow-up time to each subject in follow-up back to orignial timeline.
    ####################################
    #interimData.filled$exit[is.na(interimData$exit)]<-exit
    interimData.filled$exit[is.na(interimData$exit)]<-interimData.filled$entry[is.na(interimData$exit)]+adj.fUP+exit

    interimData.filled$event[is.na(interimData$exit)]<-event
    interimData.filled$dropout[is.na(interimData$exit)]<-droppedout

    # generate the indicator of membership in the per-protocol cohort
    if (ppAnalysis){
      if (!is.null(Seed)){ set.seed(Seed+10000) }
      interimData.filled$missVacc <- ifelse(interimData$missVacc==0 & is.na(interimData$pp), rbinom(NROW(interimData), 1, prob=missVaccProb), interimData$missVacc)
      interimData.filled$pp <- ifelse(is.na(interimData$pp), as.numeric(interimData$missVacc==0 & interimData$exit - interimData$entry > ppAtRiskTimePoint), interimData$pp)
    }

    out <- rbind(interimData.filled, out)
    out$missVacc <- NULL

    return(out)
  }

#' Treatment Arm-Pooled Simulation-Based Completion of a Randomized Efficacy Trial with a Time-to-Event Endpoint and Fixed Follow-up Using an Interim Data-set
#'
#' Considers MITT data collected through an interim timepoint and generates independent time-to-event data-sets, ignoring treatment assignments, to assess the distribution of the number of treatment arm-pooled endpoints
#' at the end of the follow-up period. A Bayesian model for the treatment arm-pooled endpoint rate, offering the option to specify a robust mixture prior distribution, is used for generating future data (see the vignette).
#'
#' @param interimData a data frame capturing observed MITT data at an interim timepoint that contains one row per enrolled participant in the MITT cohort and the following variables: \code{arm} (treatment arm), \code{schedule2} (an indicator that a participant follows the \code{visitSchedule2} schedule, e.g., participants who discontinue study product administration may remain in primary follow-up on a different schedule), \code{entry} (number of weeks since the reference date until the enrollment date), \code{exit} (number of weeks since the reference date until the trial exit date defined as the date of either infection diagnosis, dropout, or primary follow-up completion, whichever occurs first; \code{NA} for participants still in primary follow-up), \code{last_visit_dt} (number of weeks since the reference date until the last visit date), \code{event} (event indicator), \code{dropout} (dropout indicator), \code{complete} (indicator of completed follow-up), \code{followup} (indicator of being in primary follow-up). The reference date is defined as the enrollment date of the first participant. The variables \code{entry}, \code{exit}, and \code{last_visit_dt} use week as the unit of time. Month is defined as 52/12 weeks.
#' @param nTrials the number of trials to be simulated
#' @param N the total target number of enrolled participants
#' @param enrollRate a treatment arm-pooled weekly enrollment rate used for completing enrollment if fewer than \code{N} participants were enrolled in \code{interimData}. If \code{NULL} (default), the rate is calculated as the average over the last \code{enrollRatePeriod} weeks of enrollment in \code{interimData}. If equal to a numerical value, then \code{enrollRatePeriod} is ignored.
#' @param enrollRatePeriod the length (in weeks) of the time period preceding the time of the last enrolled participant in \code{interimData} that the average weekly enrollment rate will be based on and used for completing enrollment. If \code{NULL} (default), then \code{enrollRate} must be specified.
#' @param eventPriorWeight a numeric value in \eqn{[0,1]} representing a weight assigned to the prior gamma distribution of the treatment arm-pooled event rate at the time when 50\% of the estimated total person-time at risk has been accumulated (see the vignette)
#' @param eventPriorRate a numeric value of a treatment arm-pooled prior mean incidence rate for the endpoint, expressed as the number of events per person-year at risk. If \code{NULL} (default), then use the observed rate in \code{interimData}.
#' @param fixedDropOutRate the pre-trial assumed annual dropout rate. If \code{NULL} (default), then the observed treatment arm-pooled dropout rate is used.
#' @param ppAnalysis a logical value (\code{FALSE} by default) indicating whether an indicator of membership in the per-protocol cohort shall be generated based on complete MITT data. If \code{TRUE}, then \code{interimData} must include two additional variables: \code{missVacc} (an indicator of a missed vaccination) and \code{pp} (an indicator of membership in the per-protocol cohort; \code{NA} for participants with an indeterminate status).
#' @param missVaccProb a probability that a participant misses at least one vaccination. If \code{NULL} (default) and \code{ppAnalysis=TRUE}, then \code{missVaccProb} is calculated as the sample proportion of MITT participants in \code{interimData} with a missed vaccination using the \code{missVacc} variable. If \code{ppAnalysis=TRUE}, then the indicator of a missed vaccination for participants in \code{interimData} with \code{pp=NA} and future enrolled participants is sampled from the Bernoulli distribution with probability \code{missVaccProb}.
#' @param ppAtRiskTimePoint a minimal follow-up time (in weeks) for a participant to qualify for inclusion in the per-protocol cohort (\code{NULL} by default)
#' @param fuTime a follow-up time (in weeks) of each participant
#' @param mixture a logical value indicating whether to use the robust mixture approach (see the vignette). If equal to \code{FALSE} (default), then \code{mix.weights} and \code{eventPriorWeightRobust} are ignored.
#' @param mix.weights a numeric vector of length 2 representing prior weights (values in \eqn{[0,1]}) of the informative and the weakly informative component, respectively, of the prior gamma-mixture distribution of the treatment arm-pooled event rate. The two weights must sum up to 1. If \code{NULL} (default) and \code{mixture=TRUE}, then \code{c(0.8,0.2)} is used.
#' @param eventPriorWeightRobust a numeric value representing the weight \eqn{w} used to calculate the \eqn{\beta} parameter of the weakly informative gamma distribution in the mixture prior. If \code{NULL} (default) and \code{mixture=TRUE}, then \eqn{1/200} is used.
#' @param visitSchedule a numeric vector of visit weeks at which testing for the endpoint is conducted
#' @param visitSchedule2 a numeric vector of visit weeks at which testing for the endpoint is conducted in a subset of participants (e.g., those who discontinue administration of the study product but remain in follow-up). If \code{NULL} (default), everyone is assumed to follow \code{visitSchedule}.
#' @param saveFile a character string specifying an \code{.RData} file storing the output list. If \code{NULL} and \code{saveDir} is specified, the file name will be generated. If, in turn, \code{saveFile} is specified but \code{saveDir} equals \code{NULL}, then \code{saveFile} is ignored, and the output list will be returned.
#' @param saveDir a character string specifying a path for the output directory. If supplied, the output is saved as an \code{.RData} file in the directory; otherwise the output is returned as a list.
#' @param randomSeed seed of the random number generator for simulation reproducibility
#'
#' @return If \code{saveDir} is specified, the output list (named \code{trialObj}) is saved as an \code{.RData} file; otherwise it is returned. The output object is a list with the following components:
#' \itemize{
#' \item \code{trialData}: a list with \code{nTrials} components each of which is a \code{data.frame} with the variables \code{arm}, \code{entry}, \code{exit}, \code{event}, and \code{dropout} storing the treatment assignments, enrollment times, study exit times, event indicators, and dropout indicators respectively. The observed follow-up times can be recovered as \code{exit} - \code{entry}. If \code{ppAnalysis=TRUE}, then the indicators of belonging to the per-protocol cohort (named \code{pp}) are included.
#' \item \code{nTrials}: the number of simulated trials
#' \item \code{N}: the total number of enrolled trial participants
#' \item \code{rates}: a list with three components:
#' \itemize{
#' \item \code{enrollRate}: the treatment arm-pooled \emph{weekly} enrollment rate
#' \item \code{dropRate}: \code{fixedDropOutRate}, or, if \code{NULL}, the \emph{annual} treatment arm-pooled dropout rate in \code{interimData}
#' \item \code{eventPostRate}: a numeric vector of length \code{nTrials} of the treatment arm-pooled \emph{annual} event rates sampled from the posterior distribution
#' }
#' \item \code{BetaOverBetaPlusTk}: the weight placed on the prior mean event rate
#' \item \code{TkOverTstar}: the ratio of the observed person-time at risk to the estimated total person-time at risk, with the event rate set equal to \code{eventPriorRate} in the estimator for the total person-time at risk
#' \item \code{randomSeed}: seed of the random number generator for simulation reproducibility
#' \item \code{w.post}: the weights, summing up to 1, of the gamma components of the posterior mixture distribution of the treatment arm-pooled event rate. If \code{mixture=FALSE}, then \code{w.post=NA}.
#' }
#'
#' @examples
#' arm <- rep(c("C3","T1","T2"), each=250)
#' schedule <- rbinom(length(arm), 1, 0.01)
#' entry <- rpois(length(arm), lambda=60)
#' entry <- entry - min(entry)
#' last_visit_dt <- entry + runif(length(arm), min=0, max=80)
#' event <- rbinom(length(arm), 1, 0.01)
#' dropout <- rbinom(length(arm), 1, 0.02)
#' dropout[event==1] <- 0
#' exit <- rep(NA, length(arm))
#' exit[event==1] <- last_visit_dt[event==1] + 5
#' exit[dropout==1] <- last_visit_dt[dropout==1] + 5
#' followup <- ifelse(event==1 | dropout==1, 0, 1)
#' interimData <- data.frame(arm=arm, schedule2=schedule, entry=entry, exit=exit,
#' last_visit_dt=last_visit_dt, event=event, dropout=dropout, complete=0,
#' followup=followup)
#'
#' completeData <- completeTrial.pooledArms(interimData=interimData, nTrials=5, N=1500,
#' enrollRatePeriod=24, eventPriorWeight=0.5, eventPriorRate=0.001, fuTime=80,
#' visitSchedule=seq(0, 80, by=4),
#' visitSchedule2=c(0,seq(from=8,to=80,by=12)), randomSeed=9)
#' ### alternatively, to save the .RData output file (no '<-' needed):
#' completeTrial.pooledArms(interimData=interimData, nTrials=5, N=1500,
#' enrollRatePeriod=24, eventPriorWeight=0.5, eventPriorRate=0.001, fuTime=80,
#' visitSchedule=seq(0, 80, by=4),
#' visitSchedule2=c(0,seq(from=8,to=80,by=12)), saveDir="./", randomSeed=9)
#'
#' @seealso \code{\link{completeTrial.byArm}}
#' @export
completeTrial.pooledArms <-
  function(
    #in interimData contains variables: "arm","schedule","entry","exit","last_visit_dt","event","dropout","complete","followup"
    #entry, exit and last_visit_dt are in unit of weeks
    #event, dropout, complete, followup are indicators, one of these four indicators equals 1
    interimData,

    # number of trials to simulate
    nTrials,

    # the targeted total number of enrolled participants
    N,

    # a user-specified enrollment rate. If enrollRate equals a numerical value, then parameter enrollRatePeriod is ignored
    enrollRate = NULL,

    # the period that the enrollment rate will be calculated based on, in weeks, when enrollRate equals null. Parameters 'enrollRate' and 'enrollRatePeriod' can not be both null.
    enrollRatePeriod = NULL,

    # the weight for event rate Prior
    eventPriorWeight,

    # pre-trial assumed event rate, needs to be in unit of # per person-year
    # if eventPriorRate=NULL, then use observed event rate
    eventPriorRate = NULL,

    # pre-trial assumptions on annual drop-out rate. If NULL, the observed drop-out rate is used.
    fixedDropOutRate=NULL,

    # should the PP indicator be generated?
    ppAnalysis = FALSE,

    #used to create "per-protocol" indicators
    #If specified, indicator for belonging to a per-protocol cohort is created
    missVaccProb = NULL,

    ppAtRiskTimePoint = NULL,

    fuTime,

    # this is a dummy variable to call the robust mixture approach
    mixture=FALSE,

    # this is a vector of lenght 2 to indicate the weights of the
    # informative part and of the uniformative part (default is 0.80/0.20).
    # the elements should sum to 1.
    mix.weights=NULL,

    # the robust / non-informative parts is specified via a weight
    eventPriorWeightRobust=NULL,

    #Schedule 1-Procedures at CRS for event-free participants
    visitSchedule,

    #Schedule 4-Procedures at CRS for participants who discontinue infusions for reasons other than event
    # if NULL, everyone is assumed to follow visitSchedule
    visitSchedule2 = NULL,

    saveFile= NULL,
    saveDir = NULL,
    randomSeed = NULL ){
    ## define 'eps' (epsilon) a fudge-factor used where we're concerned with the limits
    ## of floating point accuracy
    eps <- sqrt( .Machine$double.eps )

    ## total number of trial participants
    Nppt <- N

    if(!is.null(enrollRate)){
      enrollmentRate <- enrollRate
    } else {
      if(is.null(enrollRatePeriod)){
        stop("Arguments 'enrollRate' and 'enrollRatePeriod' cannot be both NULL.")
      }

      # calculate weekly enrollment rate based on interimData and enrollRatePeriod, rate=enrolled participants per week
      enrollmentRate <- sum(interimData$entry >= (max(interimData$entry) - enrollRatePeriod)) / enrollRatePeriod
    }

    # a vector of individual follow-up times in weeks
    tFU <- ifelse(interimData$followup==1, interimData$last_visit_dt - interimData$entry, interimData$exit - interimData$entry)
    # the total observed person-weeks at risk in the interim data
    totFU <- sum(tFU, na.rm=TRUE)

    n_k <- sum(interimData$event)
    T_k <- totFU

    if (is.null(eventPriorRate)) {
      # calculate weekly event rate based on observed event rate
      # T_k is in weeks
      eventPriorRate <- n_k / T_k #events/person-week
    } else {
      #change unit of per person-year to per person-week
      eventPriorRate <- eventPriorRate / 52
    }

    #Estimation of Total Person-Weeks at Risk (T_star)
    #pre-trial assumed dropout rate d_star=0.1 per person-year=0.1/52 per person-week
    if(!is.null(fixedDropOutRate)){
      dropRate <- fixedDropOutRate / 52
    } else {
      #calculate weekly dropout rate based on interimData, rate=observed dropout per person-week
      # 'totFU' is in weeks
      dropRate <- sum(interimData$dropout) / totFU
    }

    d_star <- dropRate
    # estimate the total person-weeks at risk in the completed data using the user-specified event rate, if available, or else the observed event rate in the interim data
    T_star <- N * (1 - exp(-(d_star + eventPriorRate) * fuTime)) / (d_star + eventPriorRate)

    ## if mixture==TRUE then alpha and beta should become vectors
    ## furthermore the weights need to be updated as well.

    if(mixture){
        if(is.null(eventPriorWeightRobust)){
            warning("'eventPriorWeightRobust' unspecified and thus set to 1/200.")
            eventPriorWeight <- c(eventPriorWeight, 1/200)
        }else{
         eventPriorWeight <- c(eventPriorWeight,eventPriorWeightRobust)
        }
    }

    ## if mixture=FALSE than these remain scalars
    ## because the two gamma have the same expected value, eventPriorRate can be used
    beta <- T_star * eventPriorWeight / (2 * (1 - eventPriorWeight))
    alpha <- beta * eventPriorRate

    ## create lists for storage of trial data
    trialList  <- vector("list", nTrials)

    #record the weekly event rates
    eventPostRate <- rep(NA, nTrials)

    ## core of the mixture part
    if(mixture){
        # verify that two weights are specified
        if(is.null(mix.weights)) {
            # first element is assigned a weight of 0.8
            warning("mixing weights not specified, set to c(0.8,0.2)")
            mix.weights=c(0.8,0.2)
        }
        if(length(mix.weights)!=2){
            warning("mixing weights not specified, set to c(0.8,0.2)")
            mix.weights=c(0.8,0.2)
        }

        # determine how likely the informative and uninformative part are
        # via the marginal likelihoods
        marg_lik <- numeric(2)
        marg_lik[1] <- dgamma(n_k/T_k, alpha[1]+n_k, beta[1]+T_k)
        marg_lik[2] <- dgamma(n_k/T_k, alpha[2]+n_k, beta[2]+T_k)
        # posterior weights
        w_post <- marg_lik * mix.weights / sum(marg_lik *mix.weights)
    }

    for (i in 1:nTrials){
      if (!is.null(randomSeed)){ set.seed(randomSeed+i) }
      #weekly event rate
      #event rate=rgamma( 1,  alpha + n_Events,  beta + totFU_for_event)
      if(mixture){
           ## simple sampling from the mixture distribution
          ## sample runif from 0-1; if smaller than max-weight -> sample from dist.
          ## with max weight otherwise sample from the other distribution.
          tmp <- runif(1, min=0, max=1)
          IDX <- which(w_post==max(w_post))
          if(tmp < max(w_post)) {
            eventRate <- rgamma(1, alpha[IDX]+n_k, beta[IDX]+T_k)
          }else{
              #IDX <- which(w_post!=max(w_post))
            eventRate <- rgamma(1, alpha[-IDX]+n_k, beta[-IDX]+T_k)
          }
      } else {
        # sample the weekly event rate from the posterior distribution
        eventRate <- rgamma(1, alpha + n_k, beta + T_k)
      }

      ## 'parSet' contains weekly rates
      rates <- list(enrollmentRate=enrollmentRate, eventRate=eventRate, dropRate=dropRate)

      if (ppAnalysis){
        if (is.null(missVaccProb)){
          missVaccProb <- sum(interimData$missvac, na.rm=TRUE) / NROW(interimData)
        }
      }

      ## generate data
      trialList[[i]] <- FillinInterimdata.Pooled(interimData=interimData,
                                      rates = rates,
                                      visitSchedule = visitSchedule,
                                      visitSchedule2 = visitSchedule2,
                                      Nppt = Nppt,
                                      fuTime = fuTime,
                                      ppAnalysis = ppAnalysis,
                                      missVaccProb = missVaccProb,
                                      ppAtRiskTimePoint = ppAtRiskTimePoint,
                                      Seed = randomSeed+i)

      ## store sampled posterior event rate
      eventPostRate[i] <- eventRate
    }

    ## Put everything into a "trial Object"
    trialObj <- list(trialData = trialList, nTrials = nTrials, N = Nppt,
                     rates = list(enrollRate=enrollmentRate, dropRate=dropRate * 52, eventPostRate=eventPostRate * 52),
                     BetaOverBetaPlusTk = beta / (beta + T_k),
                     TkOverTstar = T_k / T_star,
                     randomSeed = randomSeed,
                     w.post = ifelse(mixture, w_post, NA))

    # save trial output and information on used rates
    if (!is.null(saveDir)){
      if (is.null(saveFile)){
        # in the file name, eventPriorRate shows the per person-YEAR rate
        saveFile <- paste0("completeTrial_pooled_eventPriorRate=", round(eventPriorRate * 52, 4), "_eventPriorWt=", round(eventPriorWeight, 3), ".RData")
      }

      save(trialObj, file=file.path(saveDir, saveFile))
    }

    return(invisible(trialObj))

  }   ######################## End of completeTrial.pooledArms function #######################


#' Treatment Arm-Specific Simulation-Based Completion of a Randomized Efficacy Trial with a Time-to-Event Endpoint and Fixed Follow-up Using an Interim Data-set
#'
#' Considers MITT data collected through an interim timepoint and generates independent time-to-event data-sets, by treatment arm, to assess the distribution of the number of treatment arm-specific endpoints
#' at the end of the follow-up period. A Bayesian model for treatment arm-specific endpoint rates is used for generating future data (see the vignette).
#'
#' @param interimData a data frame capturing observed MITT data at an interim timepoint that contains one row per enrolled participant in the MITT cohort and the following variables: \code{arm} (treatment arm), \code{schedule2} (an indicator that a participant follows the \code{visitSchedule2} schedule, e.g., participants who discontinue study product administration may remain in primary follow-up on a different schedule), \code{entry} (number of weeks since the reference date until the enrollment date), \code{exit} (number of weeks since the reference date until the trial exit date defined as the date of either infection diagnosis, dropout, or primary follow-up completion, whichever occurs first; \code{NA} for participants still in primary follow-up), \code{last_visit_dt} (number of weeks since the reference date until the last visit date), \code{event} (event indicator), \code{dropout} (dropout indicator), \code{complete} (indicator of completed follow-up), \code{followup} (indicator of being in primary follow-up). The reference date is defined as the enrollment date of the first participant. The variables \code{entry}, \code{exit}, and \code{last_visit_dt} use week as the unit of time. Month is defined as 52/12 weeks.
#' @param nTrials the number of trials to be simulated
#' @param trtNames a character vector of treatment labels as specified in \code{interimData$arm} determining the order of treatment arms in other input arguments
#' @param N a numeric vector specifying the target number of enrolled participants in each treatment arm, with the arms in the same order as in \code{trtNames}
#' @param enrollRate a treatment arm-pooled weekly enrollment rate used for completing enrollment if \code{interimData}'s enrollment is incomplete. If \code{NULL} (default), the rate is calculated as the average over the last \code{enrollRatePeriod} weeks of enrollment in \code{interimData}. If equal to a numeric value, then \code{enrollRatePeriod} is ignored.
#' @param enrollRatePeriod the length (in weeks) of the time period preceding the time of the last enrolled participant in \code{interimData} that the average weekly enrollment rate will be based on and used for completing enrollment. If \code{NULL} (default), then \code{enrollRate} must be specified.
#' @param eventPriorWeight a numeric value in \eqn{[0,1]} representing a weight assigned to the prior gamma distribution of the treatment arm-specific event rates at the time when 50\% of the estimated person-time at risk in each arm has been accumulated (see the vignette)
#' @param eventPriorRate a numeric vector of treatment arm-specific prior mean incidence rates for the endpoint, expressed as numbers of events per person-year at risk, with the arms in the same order as in \code{trtNames}
#' @param fixedDropOutRate the pre-trial assumed annual treatment arm-pooled dropout rate. If \code{NULL} (default), then the observed treatment arm-pooled dropout rate is used.
#' @param ppAnalysis a logical value (\code{FALSE} by default) indicating whether an indicator of membership in the per-protocol cohort shall be generated based on complete MITT data. If \code{TRUE}, then \code{interimData} must include two additional variables: \code{missVacc} (an indicator of a missed vaccination) and \code{pp} (an indicator of membership in the per-protocol cohort; \code{NA} for participants with an indeterminate status).
#' @param missVaccProb a probability that a participant misses at least one vaccination. If \code{NULL} (default) and \code{ppAnalysis=TRUE}, then \code{missVaccProb} is calculated as the sample proportion of MITT participants in \code{interimData} with a missed vaccination using the \code{missVacc} variable. If \code{ppAnalysis=TRUE}, then the indicator of a missed vaccination for participants in \code{interimData} with \code{pp=NA} and future enrolled participants is sampled from the Bernoulli distribution with probability \code{missVaccProb}.
#' @param ppAtRiskTimePoint a minimal follow-up time (in weeks) for a participant to qualify for inclusion in the per-protocol cohort (\code{NULL} by default)
#' @param fuTime a follow-up time (in weeks) of each participant
#' @param visitSchedule a numeric vector of visit weeks at which testing for the endpoint is conducted
#' @param visitSchedule2 a numeric vector of visit weeks at which testing for the endpoint is conducted in a subset of participants (e.g., those who discontinue administration of the study product but remain in follow-up). If \code{NULL} (default), everyone is assumed to follow \code{visitSchedule}.
#' @param saveFile a character string specifying an \code{.RData} file storing the output list. If \code{NULL} and \code{saveDir} is specified, the file name will be generated. If, in turn, \code{saveFile} is specified but \code{saveDir} equals \code{NULL}, then \code{saveFile} is ignored, and the output list will be returned.
#' @param saveDir a character string specifying a path for the output directory. If supplied, the output is saved as an \code{.RData} file in the directory; otherwise the output is returned as a list.
#' @param randomSeed seed of the random number generator for simulation reproducibility
#'
#' @return If \code{saveDir} is specified, the output list (named \code{trialObj}) is saved as an \code{.RData} file; otherwise it is returned. The output object is a list with the following components:
#' \itemize{
#' \item \code{trialData}: a list with \code{nTrials} components each of which is a \code{data.frame} with the variables \code{arm}, \code{entry}, \code{exit}, \code{event}, and \code{dropout} storing the treatment assignments, enrollment times, study exit times, event indicators, and dropout indicators, respectively. The observed follow-up times can be recovered as \code{exit} - \code{entry}. If \code{ppAnalysis=TRUE}, then the indicators of belonging to the per-protocol cohort (named \code{pp}) are included.
#' \item \code{nTrials}: the number of simulated trials
#' \item \code{N}: the total number of enrolled trial participants
#' \item \code{rates}: a list with three components:
#' \itemize{
#' \item \code{enrollRate}: the treatment arm-pooled \emph{weekly} enrollment rate
#' \item \code{dropRate}: \code{fixedDropOutRate}, or, if \code{NULL}, the \emph{annual} treatment arm-pooled dropout rate in \code{interimData}
#' \item \code{eventPostRate}: a list with \code{length(trtNames)} components (labeled by the levels of the \code{arm} variable in \code{interimData}) each of which is a numeric vector of length \code{nTrials} of the sampled treatment arm-specific posterior \emph{annual} event rates
#' }
#' \item \code{BetaOverBetaPlusTk}: a list with \code{length(trtNames)} components (labeled by the levels of the \code{arm} variable in \code{interimData}) each of which is the arm-specific weight placed on the prior mean event rate
#' \item \code{TkOverTstar}: a list with \code{length(trtNames)} components (labeled by the levels of the \code{arm} variable in \code{interimData}) each of which is the ratio of the observed arm-specific person-time at risk to the estimated total arm-specific person-time at risk, with the arm-specific event rates set equal to the components of \code{eventPriorRate} in the estimator for the total arm-specific person-time at risk
#' \item \code{randomSeed}: seed of the random number generator for simulation reproducibility
#' }
#'
#' @examples
#' arm <- rep(c("C3","T1","T2"), each=250)
#' schedule <- rbinom(length(arm), 1, 0.01)
#' entry <- rpois(length(arm), lambda=60)
#' entry <- entry - min(entry)
#' last_visit_dt <- entry + runif(length(arm), min=0, max=80)
#' event <- rbinom(length(arm), 1, 0.01)
#' dropout <- rbinom(length(arm), 1, 0.02)
#' dropout[event==1] <- 0
#' exit <- rep(NA, length(arm))
#' exit[event==1] <- last_visit_dt[event==1] + 5
#' exit[dropout==1] <- last_visit_dt[dropout==1] + 5
#' followup <- ifelse(event==1 | dropout==1, 0, 1)
#' interimData <- data.frame(arm=arm, schedule2=schedule, entry=entry, exit=exit,
#' last_visit_dt=last_visit_dt, event=event, dropout=dropout, complete=0, followup=followup)
#'
#' completeData <- completeTrial.byArm(interimData=interimData, nTrials=5,
#' trtNames=c("C3","T1","T2"), N=c(500,500,500), enrollRatePeriod=24, eventPriorWeight=0.5,
#' eventPriorRate=c(0.001,0.0004,0.0004), fuTime=80, visitSchedule=seq(0, 80, by=4),
#' visitSchedule2=c(0,seq(from=8,to=80,by=12)), randomSeed=9)
#' ### alternatively, to save the .RData output file (no '<-' needed):
#' completeTrial.byArm(interimData=interimData, nTrials=5, trtNames=c("C3","T1","T2"),
#' N=c(500,500,500), enrollRatePeriod=24, eventPriorWeight=0.5,
#' eventPriorRate=c(0.001,0.0004,0.0004), fuTime=80, visitSchedule=seq(0, 80, by=4),
#' visitSchedule2=c(0,seq(from=8,to=80,by=12)), saveDir="./", randomSeed=9)
#'
#' @seealso \code{\link{completeTrial.pooledArms}}
#' @export
completeTrial.byArm <-
  function(
    #in interimData contains variables: "arm","schedule","entry","exit","last_visit_dt","event","dropout","complete","followup"
    #entry, exit and last_visit_dt are in unit of weeks
    #event, dropout, complete, followup are indicators, one of these four indicators equals 1
    interimData,

    # number of trials to simulate
    nTrials,

    ## treatment arm names: a vector of character strings, one per treatment arm
    trtNames,

    # number of subjects in in each treatment arm, same order as trtNames
    N,

    #a user-specified arm-pooled enrollment rate vector. If enrollRate equals a numerical value, then parameter enrollRatePeriod is ignored
    enrollRate = NULL,

    # the period that the enrollment rate will be calculated based on, in weeks
    enrollRatePeriod,

    eventPriorWeight,
    # pre-trial assumed event rate, needs to be in unit of # per person-year
    # a vector for each treatment arm, same order as trtNames
    # eventPriorRate can not be null
    eventPriorRate,

    # pre-trial assumptions on annual drop-out rate. If NULL, the observed drop-out rate is used.
    fixedDropOutRate=NULL,

    # should the PP indicator be generated?
    ppAnalysis = FALSE,

    #used to create "per-protocol" indicators
    #If specified, indicator for belonging to a per-protocol cohort is created
    missVaccProb = NULL,

    ppAtRiskTimePoint = NULL,

    fuTime,

    #Schedule 1-Procedures at CRS for event-free participants
    visitSchedule,

    #Schedule 4-Procedures at CRS for participants who discontinue infusions for reasons other than event
    #if Null, everyone is assumed to follow visitSchedule
    visitSchedule2 = NULL,

    saveFile= NULL,
    saveDir = NULL,
    randomSeed = NULL ){
    ## define 'eps' (epsilon) a fudge-factor used where we're concerned with the limits
    ## of floating point accuracy
    eps <- sqrt( .Machine$double.eps )

    if (!is.null(enrollRate)){
      enrollmentRate <- enrollRate
    } else {
      if (is.null(enrollRatePeriod)){
        stop("Arguments 'enrollRate' and 'enrollRatePeriod' cannot be both NULL.")
      }

      #calculate arm-pooled weekly enrollment rate based on interimData and enrollRatePeriod, rate=enrolled participants per week
      enrollmentRate <- sum(interimData$entry >= (max(interimData$entry) - enrollRatePeriod)) / enrollRatePeriod
    }

    # a vector of individual follow-up times in weeks
    tFU <- ifelse(interimData$followup==1, interimData$last_visit_dt - interimData$entry, interimData$exit - interimData$entry)
    # the total observed person-weeks at risk in the interim data
    totFU <- sum(tFU, na.rm=TRUE)

    ## check 'trtNames' contains the same elements as interimData$arm
    if (!setequal(trtNames, interimData$arm)){
      stop("The input argument 'trtNames' and the variable 'arm' in 'interimData' must contain the same unique elements.")
    }

    nArms <- length(N)

    # a vector of arm-specifc event counts
    n_k <- NULL
    for(i in 1:nArms){
      n_k <- c(n_k, sum(interimData$event[interimData$arm==trtNames[i]]))
    }

    # a vector of arm-specific person-weeks at risk
    T_k <- NULL
    for(i in 1:nArms){
      T_k <- c(T_k, sum(tFU[interimData$arm==trtNames[i]], na.rm=TRUE))
    }

    if (is.null(eventPriorRate)) {
      # calculate weekly arm-specific event rates based on observed arm-specific event rates
      # T_k is in weeks
      eventPriorRate <- n_k / T_k #events/person-week
    } else {
      #change unit of per person-year to per person-week
      eventPriorRate <- eventPriorRate / 52
    }

    #Estimation of Arm Specific Total Person-Weeks at Risk (T_star)
    #pre-trial assumed dropout rate d_star=0.1 per person-year=0.1/52 per person-week
    if (!is.null(fixedDropOutRate)){
      dropRate <- fixedDropOutRate / 52
    } else {
      #calculate weekly dropout rate based on interimData, rate=observed dropout per person-week
      # 'totFU' is in weeks
      dropRate <- sum(interimData$dropout) / totFU
    }

    d_star <- dropRate
    # estimate the arm-specific person-weeks at risk in the completed data using the user-specified arm-specific event rates, if available, or else the observed event rates in the interim data
    T_star <- N * (1 - exp(-(d_star + eventPriorRate) * fuTime)) / (d_star + eventPriorRate)

    beta <- T_star * eventPriorWeight / (2 * (1 - eventPriorWeight))
    alpha <- beta * eventPriorRate

    ## create lists for storage of trial data
    trialList  <- vector("list", nTrials)

    #record arm-specific weekly event rates
    eventPostRate    <- vector("list", nArms)
    for (i in 1:nArms){
      eventPostRate[[i]] <- rep(NA, nTrials)
    }

    for (i in 1:nTrials){
      #arm-specific weekly event rates
      eventRate <- NULL
      for (j in 1:nArms){
        if (!is.null(randomSeed)){ set.seed(randomSeed+i) }
        # sample the weekly event rate from the posterior distribution
        eventRate <- c(eventRate, rgamma(1, alpha[j] + n_k[j],  beta[j] + T_k[j]))
      }

      ## 'parSet' contains all arm-specific weekly rates
      rates <- list(enrollmentRate=enrollmentRate, dropRate=dropRate, eventRate=eventRate)

      if (ppAnalysis){
        if (is.null(missVaccProb)){
          # the sample proportion of participants with at least 1 missed vaccination pooling over all arms
          missVaccProb <- sum(interimData$missvac, na.rm=TRUE) / NROW(interimData)
        }
      }

      ## generate data
      trialList[[i]] <- FillinInterimdata.byArm(interimData=interimData,
                                                rates = rates,
                                                visitSchedule = visitSchedule,
                                                visitSchedule2 = visitSchedule2,
                                                trtNames=trtNames,
                                                N = N,
                                                fuTime = fuTime,
                                                ppAnalysis = ppAnalysis,
                                                missVaccProb = missVaccProb,
                                                ppAtRiskTimePoint = ppAtRiskTimePoint,
                                                Seed = randomSeed+i)

      ## store sampled posterior arm-specific event rates
      for (j in 1:nArms){
        eventPostRate[[j]][i] <- eventRate[j]
      }
    }

    ## Put everything into a "trial Object"
    trialObj <- list(trialData = trialList, nTrials = nTrials, N = N,
                     rates = list(enrollRate=enrollmentRate, dropRate=dropRate * 52, eventPostRate=lapply(eventPostRate, function(r){ r * 52 })),
                     BetaOverBetaPlusTk = beta / (beta+T_k),
                     TkOverTstar = T_k / T_star,
                     randomSeed = randomSeed)

    # save trial output and information on used rates
    if (!is.null(saveDir)){
      if (is.null(saveFile)){

        # in the file name, eventPriorRate shows the per person-YEAR rate
        saveFile <- paste0("completeTrial_byArm_", paste(paste0("eventPriorRate", trtNames, "=", round(eventPriorRate * 52, 3)), collapse="_"), "_eventPriorWt=", round(eventPriorWeight, 3), ".RData")
      }

      save(trialObj, file=file.path(saveDir, saveFile))
    }

    return(invisible(trialObj))
  }   ######################## End of completeTrial.byArm function #######################


to_rcdf <- function(s, min){
  var.floored <- pmax(min,s) # Apply floor function at minimum of x
  var.sorted <- sort(unique(var.floored)) # Get all unique values of x, sorted (after applying minimum)
  var.rcdf <- 1 - ecdf(var.floored)(var.sorted) # This gives the proportion of observations > x
  var.rcdf <- c(1,var.rcdf[-length(var.rcdf)]) # This gives the proportion of observations >= x
  out <- data.frame(var.sorted,var.rcdf)
  names(out) <- c("var","RCDF")
  return(out)
} # to_rcd


#' Plot Characteristics of the Estimated Distribution of the Treatment Arm-Pooled Number of Endpoints
#'
#' Takes the output from the \code{\link{completeTrial.pooledArms}} function and generates a plot describing characteristics of the estimated distribution of the treatment arm-pooled number of endpoints.
#'
#' @param eventTimeFrame a time frame within which endpoints are counted, specified in weeks as \code{c(start, end)}. If \code{NULL} (default), then all endpoints are counted.
#' @param eventPPcohort a logical value. If \code{TRUE}, only endpoints in the per-protocol cohort are counted. The default value is \code{FALSE}.
#' @param target a vector of target numbers of endpoints for reporting of the estimated probability that the total number of endpoints will be \eqn{\ge} \code{target}, with a 95\% credible interval
#' @param power.axis a logical value. If \code{TRUE} (default), then a top axis is added to the plot, showing power to reject \eqn{H_0}: TE \eqn{\le} 0\% using a 1-sided 0.025-level Wald test if TE = \code{power.TE} throughout the trial.
#' @param power.TE a numeric value of treatment efficacy for which power is shown on the top axis. If \code{power.axis} is \code{FALSE}, then \code{power.TE} is ignored.
#' @param eventPriorRate a numeric value of the treatment arm-pooled prior mean incidence rate for the endpoint, expressed as the number of events per person-year at risk
#' @param eventPriorWeight a numeric vector in which each value represents a weight (i.e., a separate scenario) assigned to the prior gamma distribution of the treatment arm-pooled event rate at the time when 50\% of the estimated total person-time at risk has been accumulated
#' @param xlim a numeric vector of the form \code{c(xmin, xmax)} for the user-specified x-axis limits. If \code{NULL} (default), then the computed range of x-axis values will be used.
#' @param xlab a character string for the user-specified x-axis label. If \code{NULL} (default), then the label "Total Number of Infections (n)" will be used.
#' @param ylab a character string for the user-specified y-axis label. If \code{NULL} (default), then the label "P( Total Number of Infections >= n ) x 100" will be used.
#' @param power.lab a character string for the user-specified power-axis label. If \code{NULL} (default), then the label "Power for TE = \code{power.TE} (x 100)" will be used.
#' @param xPosLegend a numeric value in \eqn{[0,1]} (0.67 by default) specifying the x-coordinate for the position of the legend
#' @param fileDir a character string specifying a path for the input directory
#'
#' @return None. The function is called solely for plot generation.
#'
#' @examples
#' arm <- rep(c("C3","T1","T2"), each=250)
#' schedule <- rbinom(length(arm), 1, 0.01)
#' entry <- rpois(length(arm), lambda=60)
#' entry <- entry - min(entry)
#' last_visit_dt <- entry + runif(length(arm), min=0, max=80)
#' event <- rbinom(length(arm), 1, 0.01)
#' dropout <- rbinom(length(arm), 1, 0.02)
#' dropout[event==1] <- 0
#' exit <- rep(NA, length(arm))
#' exit[event==1] <- last_visit_dt[event==1] + 5
#' exit[dropout==1] <- last_visit_dt[dropout==1] + 5
#' followup <- ifelse(event==1 | dropout==1, 0, 1)
#' interimData <- data.frame(arm=arm, schedule2=schedule, entry=entry, exit=exit,
#' last_visit_dt=last_visit_dt, event=event, dropout=dropout, complete=0,
#' followup=followup)
#'
#' weights <- c(0.2, 0.4, 0.6)
#' for (j in 1:length(weights)){
#'   completeTrial.pooledArms(interimData=interimData, nTrials=50, N=1500, enrollRatePeriod=24,
#'   eventPriorWeight=weights[j], eventPriorRate=0.06, fuTime=80, visitSchedule=seq(0, 80, by=4),
#'   visitSchedule2=c(0,seq(from=8,to=80,by=12)), saveDir="./", randomSeed=9)
#' }
#'
#' pdf(file=paste0("./","rcdf_pooled_eventPriorRate=",0.06,".pdf"), width=6, height=5)
#' plotRCDF.pooledArms(target=c(60,30), power.axis=FALSE, eventPriorRate=0.06,
#' eventPriorWeight=weights, fileDir="./")
#' dev.off()
#'
#' @seealso \code{\link{completeTrial.pooledArms}} and \code{\link{plotRCDF.byArm}}
#' @export
plotRCDF.pooledArms <- function(eventTimeFrame=NULL, #the time frame to count events, in a format of c(start, end), in weeks. If null, then count all events.
                                eventPPcohort=FALSE,
                                target,
                                power.axis=TRUE,
                                power.TE=NULL,
                                eventPriorRate,
                                eventPriorWeight,
                                xlim=NULL,
                                xlab=NULL,
                                ylab=NULL,
                                power.lab=NULL,
                                xPosLegend=0.67,
                                fileDir){

  dat  <- vector("list", length(eventPriorWeight))

  for (j in 1:length(eventPriorWeight)){
    wt<-eventPriorWeight[j]
    # a list named 'trialObj'
    load(file.path(fileDir, paste0("completeTrial_pooled_eventPriorRate=", eventPriorRate, "_eventPriorWt=", wt, ".RData")))
    legend.Prior.weight<-trialObj$BetaOverBetaPlusTk
    TNI<-rep(NA,length(trialObj$trialData))

    if (is.null(eventTimeFrame)){
      if (!eventPPcohort){
        for (i in 1:length(TNI)){
          TNI[i]<-sum(trialObj$trialData[[i]]$event)
        }
      } else {
        if (!"pp" %in% colnames(trialObj$trialData[[1]])){
          stop("trialData does not have a variable named 'pp'.")
        }
        for (i in 1:length(TNI)){
          TNI[i]<-sum(subset(trialObj$trialData[[i]], pp==1)$event)
        }
      }
    } else {
      if (!eventPPcohort){
        for (i in 1:length(TNI)){
          trialObj$trialData[[i]]$eventTime <- trialObj$trialData[[i]]$exit - trialObj$trialData[[i]]$entry
          TNI[i]<-sum(subset(trialObj$trialData[[i]], eventTime >= eventTimeFrame[1] & eventTime <= eventTimeFrame[2])$event)
        }
      } else {
        if (!"pp" %in% colnames(trialObj$trialData[[1]])){
          stop("trialData does not have a variable named 'pp'.")
        }
        for (i in 1:length(TNI)){
          trialObj$trialData[[i]]$eventTime <- trialObj$trialData[[i]]$exit - trialObj$trialData[[i]]$entry
          TNI[i] <- sum(subset(trialObj$trialData[[i]], eventTime >= eventTimeFrame[1] & eventTime <= eventTimeFrame[2] & pp==1)$event)
        }
      }
    }

    rm(trialObj)
    #mean and 95% confidence interval
    mu.TNI<-round(mean(TNI),0)
    CI.TNI<-round(c(mean(TNI) - qnorm(0.975)*sd(TNI)/sqrt(length(TNI)),
                    mean(TNI) + qnorm(0.975)*sd(TNI)/sqrt(length(TNI))),0)

    #estimated prob that the total number of events>=target with 95% credible intervals
    target.p  <- vector("list", length(target))
    for(k in 1:length(target)){
      p<-sum(TNI>=target[k])/length(TNI)
      p.CI<- c(p-qnorm(0.975)*sqrt(p*(1-p)/length(TNI)),
                      p+qnorm(0.975)*sqrt(p*(1-p)/length(TNI)))
      target.p[[k]]<-list(p=p,p.CI=p.CI)
    }

    # get a data frame of x- and y-coordinates of jumps in the RCDF step function
    TNI_rcdf<-to_rcdf(TNI,min=min(TNI))

    dat[[j]]<-list(TNI_rcdf=TNI_rcdf, mu.TNI=mu.TNI, CI.TNI=CI.TNI, target.p=target.p, legend.Prior.weight=legend.Prior.weight)
  }

  var.all<-dat[[1]]$TNI_rcdf$var
  if(length(dat)>1){
    for(i in 2:length(dat)){
      var.all<-c(var.all, dat[[i]]$TNI_rcdf$var)
    }
  }

  if(is.null(xlim)){
    myxlim<-c(min(var.all, target), max(var.all, target))
  }else{
    myxlim<-xlim
  }

  tmp<-seq(10,200,10)
  tmp<-tmp[tmp>=(myxlim[1]+5) & tmp<=(myxlim[2]-6)]

  x.label<-c(myxlim[1],tmp,myxlim[2])
  colors<-c("blue","red","seagreen")
  pchar <- 15:17

  #par(fig=c(0.1,1,0,1))
  par(mar=c(2.5,2.8,2.5,0.5),oma=c(1,1,1,1),las=1)
  mycex<-0.75
  mycex2<-1
  plot(myxlim, c(0,1), type="n", main="", axes=FALSE,xlab="", ylab="")

  # mean of x
  polygon(x=c(rep(min(x.label),2),rep(max(x.label),2)), y=c(0,-2,-2,0), border=NA, col="gray90")
  axis(side=1, at=min(x.label), #+0.01*(max(x.label)-min(x.label)),
       labels="Mean:", tick=FALSE, line=-2, cex.axis=0.6, hadj=0)
  tmp<-NULL
  for(i in 1:length(dat)){
    tmp<-c(tmp, dat[[i]]$mu.TNI)
  }
  widthchange<-length(tmp)!=length(unique(tmp))
  widthloc<-NULL
  for(i in 2:length(tmp)){
    if(tmp[i] %in% tmp[1:i-1]) widthloc<-i
  }

  for(j in 1:length(dat)){
    mylwd<-1.5
    if(!is.null(widthloc)){if(j==widthloc) mylwd<-0.8}
    segments(x0=dat[[j]]$mu.TNI, y0=0, y1=1, col=colors[[j]], lwd=mylwd)
  }

  tmp <- unique(tmp)

  # suppressWarnings() is used because if 'tmp' is of length 1, then diff() returns an empty vector, and 'min' issues a warning
  if(suppressWarnings(min(diff(sort(tmp))))>=3 | max(tmp)-min(tmp)==0){
    axis(side=1, at=tmp[c(1,3)], labels=tmp[c(1,3)], tick=FALSE, line=-2, cex.axis=0.5)
    if (length(tmp)>=2){ axis(side=1, at=tmp[2], labels=tmp[2], tick=FALSE, line=-2, cex.axis=0.5) }
  } else {
    # at least one of the differences is <3, i.e., 'idx' is not empty
    idx <- which(diff(sort(tmp))<3)
    if (length(idx)==1){
      if (idx==1){ if (length(tmp)==3){ atValue <- sort(tmp) + c(-1,0.5,0) } else { atValue <- sort(tmp) + c(-0.7,0.7) } }
      if (idx==2){ atValue <- sort(tmp) + c(0,-0.5,1) }
    } else {
      atValue <- sort(tmp) + c(-1,0,1)
    }
    axis(side=1, at= atValue[c(1,3)], labels=sort(tmp)[c(1,3)], tick=FALSE, line=-2, cex.axis=0.5)
    if (length(tmp)>=2){ axis(side=1, at= atValue[2], labels=sort(tmp)[2], tick=FALSE, line=-2, cex.axis=0.5) }
  }

  for(i in 1:length(dat)){
    # ensure that the step function starts at the minimum of the x-axis range
    lines(c(min(x.label),dat[[i]]$TNI_rcdf$var),c(1,dat[[i]]$TNI_rcdf$RCDF), type="S",lty=1,col=colors[i],lwd=1.2)
    points(dat[[i]]$TNI_rcdf$var,dat[[i]]$TNI_rcdf$RCDF,lty=1, pch=pchar[i], col=colors[i], cex=ifelse(i==1,0.4,0.5))
  }

  axis(side=1, at=x.label, labels=FALSE, cex.axis=mycex)
  axis(side=1, at=x.label[seq(1,length(x.label),by=2)], labels=x.label[seq(1,length(x.label),by=2)], tick=FALSE, line=-0.6, cex.axis=mycex)
  axis(side=1, at=x.label[seq(2,length(x.label),by=2)], labels=x.label[seq(2,length(x.label),by=2)], tick=FALSE, line=-0.6, cex.axis=mycex)

  yTicks <- seq(0,1,by=0.2)
  axis(side=2, at=yTicks,labels=F, cex.axis=mycex)
  axis(side=2, at=yTicks,labels=100*yTicks, tick=F, line=0, cex.axis=mycex)

  #calculate the 3 y values for target

  for(i in 1:length(target)){
    mytarget<-target[i]
    if(mytarget>=myxlim[1] & mytarget<=myxlim[2]){
      targetY<-NULL
      if (!(mytarget %in% x.label)){
        offset <- 0.019*(max(x.label)-min(x.label))
        axis(side=1, at=mytarget, labels=F,cex.axis=mycex)
        axis(side=1, at=mytarget + ifelse(mytarget-x.label[which.min(abs(mytarget-x.label))]>0, offset, -offset), labels=mytarget, tick=FALSE, line=-1.1, cex.axis=0.7)
      }

      for(j in 1:length(dat)){
        if(length(unique(dat[[j]]$TNI_rcdf$RCDF[dat[[j]]$TNI_rcdf$var==mytarget]))==1){
          targetY<-c(targetY,unique(dat[[j]]$TNI_rcdf$RCDF[dat[[j]]$TNI_rcdf$var==mytarget]))
        }else{
          if(min(dat[[j]]$TNI_rcdf$var)<mytarget){
            target2<-max(dat[[j]]$TNI_rcdf$var[dat[[j]]$TNI_rcdf$var<mytarget])
            targetY<-c(targetY,unique(dat[[j]]$TNI_rcdf$RCDF[dat[[j]]$TNI_rcdf$var==target2]))
          }else{
            targetY<-c(targetY,1)
          }
        }
      }

      segments(x0=mytarget, y0=0, y1=1.2, col="gray60", lty="dotted", lwd=1.5)
      atValue <- rep(NA,3)
      for(j in 1:length(dat)){
        lines(x=c(-10,mytarget),y=rep(targetY[j],2), col=colors[j])
        axis(side=2, at=targetY[j], labels=F,cex.axis=mycex, col=colors[j])
        # if targetY equals one of the plotted tick labels, do not duplicate the tick label
        if (!(targetY[j] %in% yTicks)){
          # if 'targetY[j]' is too close to one of the plotted tick marks
          if (min(abs(targetY[j]-yTicks))<=0.01){
            atValue[j] <- targetY[j] + ifelse(targetY[j]-yTicks[which.min(abs(targetY[j]-yTicks))]<0, -0.022, 0.022)
          } else if (min(abs(targetY[j]-yTicks))>0.01 & min(abs(targetY[j]-yTicks))<0.04){
            atValue[j] <- targetY[j] + ifelse(targetY[j]-yTicks[which.min(abs(targetY[j]-yTicks))]<0, -0.01, 0.01)
          } else {
            atValue[j] <- targetY[j]
          }
        }
      }

      idx <- which(!is.na(atValue))
      atValue <- atValue[idx]
      targetY <- targetY[idx]

      idx <- which(as.numeric(format(100*targetY, digits=1, nsmall=1))==100)
      if (length(idx)>0){
        atValue <- atValue[-idx]
        targetY <- targetY[-idx]
      }

      # if there are any additional tickmarks to plot
      if (length(atValue)>0){
        if( (max(targetY)-min(targetY)) > 0.05 | max(targetY)-min(targetY)==0){
          axis(side=2, at=atValue, labels=format(100*targetY, digits=1, nsmall=1), tick=FALSE, line=-0.4, cex.axis=0.5)
        }else{
          atValue <- atValue[order(targetY)]
          targetY <- sort(targetY)
          if (length(atValue)==2){
            if (format(100*targetY[2], digits=1, nsmall=1) != format(100*targetY[1], digits=1, nsmall=1)){
              axis(side=2, at=atValue[1]-0.01, labels=format(100*targetY[1], digits=1, nsmall=1), tick=FALSE, line=-0.4, cex.axis=0.5)
              axis(side=2, at=atValue[2]+0.01, labels=format(100*targetY[2], digits=1, nsmall=1), tick=FALSE, line=-0.4, cex.axis=0.5)
            } else {
              axis(side=2, at=atValue[1], labels=format(100*targetY[1], digits=1, nsmall=1), tick=FALSE, line=-0.4, cex.axis=0.5)
            }
          } else {
            axis(side=2, at=atValue[c(1,3)]+c(-0.01,0.01), labels=format(100*targetY[c(1,3)], digits=1, nsmall=1), tick=FALSE, line=-0.4, cex.axis=0.5)
            axis(side=2, at=atValue[2], labels=format(100*targetY[2], digits=1, nsmall=1), tick=FALSE, line=-0.4, cex.axis=0.5)
          }
        }
      }
    }
  }



  if(is.null(xlab)){mtext("Total Number of Infections (n)", side=1, las=0, line=2, cex=mycex2)
    }else{mtext(xlab, side=1, las=0, line=2, cex=mycex2)}
  if(is.null(ylab)){mtext("P( Total Number of Infections >= n ) x 100", side=2, las=0, line=2.5, cex=mycex2)
  }else{mtext(ylab, side=2, las=0, line=2.5, cex=mycex2)}

  legend((xPosLegend+0.1)*max(x.label),0.9,legend=round(c(dat[[1]]$legend.Prior.weight,dat[[2]]$legend.Prior.weight,dat[[3]]$legend.Prior.weight),2), cex=0.7, col=colors, pch=pchar,
         lty=1, bty = "n", title="Prior weight")

  legend.text <- NULL
  for (j in 1:length(dat)){
    legend.text <- c(legend.text, paste0("P(>=",target[1],") = ",format(dat[[j]]$target.p[[1]]$p*100, digits=1, nsmall=1-(dat[[j]]$target.p[[1]]$p>=0.9995)),"% (95% CI, ",format(dat[[j]]$target.p[[1]]$p.CI[1]*100, digits=1, nsmall=1-(dat[[j]]$target.p[[1]]$p.CI[1]>=0.9995))," to ",format(dat[[j]]$target.p[[1]]$p.CI[2]*100, digits=1, nsmall=1-(dat[[j]]$target.p[[1]]$p.CI[2]>=0.9995)),")","\n",
                                         "P(>=",target[2],") = ",format(dat[[j]]$target.p[[2]]$p*100, digits=1, nsmall=1-(dat[[j]]$target.p[[2]]$p>=0.9995)),"% (95% CI, ",format(dat[[j]]$target.p[[2]]$p.CI[1]*100, digits=1, nsmall=1-(dat[[j]]$target.p[[2]]$p.CI[1]>=0.9995))," to ",format(dat[[j]]$target.p[[2]]$p.CI[2]*100, digits=1, nsmall=1-(dat[[j]]$target.p[[2]]$p.CI[2]>=0.9995)),")"))
  }
  legend(xPosLegend*max(x.label), 0.7, lty=1, col=colors, pch=pchar, bty="n", cex=0.5, legend=legend.text, x.intersp=0.5, y.intersp=1.6, pt.cex=0.7)

  if(power.axis==TRUE){
    if(!is.numeric(power.TE)){
      stop("'power.TE' must be a numeric value when power.axis=TRUE.")
    }
    #Calculate power each x.label value and target1 and target2
    x.label2<-target
    Zbeta1<-sqrt(x.label*(1/3)*(2/3)*(log(1-power.TE)^2))-qnorm(1-0.025)
    x.label.power1<-100*pnorm(Zbeta1)
    Zbeta2<-sqrt(x.label2*(1/3)*(2/3)*(log(1-power.TE)^2))-qnorm(1-0.025)
    x.label.power2<-100*pnorm(Zbeta2)
    axis(side=3, at=x.label, labels=FALSE, cex.axis=mycex)
    axis(side=3, at=x.label[seq(1,length(x.label),by=2)], labels=round(x.label.power1[seq(1,length(x.label.power1),by=2)],0), tick=FALSE, line=-0.5, cex.axis=mycex)
    axis(side=3, at=x.label[seq(2,length(x.label),by=2)], labels=round(x.label.power1[seq(2,length(x.label.power1),by=2)],0), tick=FALSE, line=-0.5, cex.axis=mycex)

    for (k in 1:length(x.label2)){
      if(x.label2[k]>myxlim[1] & x.label2[k]<myxlim[2]){
        offset <- 0.019*(max(x.label)-min(x.label))
        axis(side=3, at=x.label2[k], labels=FALSE, cex.axis=mycex)
        axis(side=3, at=x.label2[k] + ifelse(x.label2[k]-x.label[which.min(abs(x.label2[k]-x.label))]>0, offset, -offset), labels=round(x.label.power2[k],0), tick=FALSE, line=-0.95, cex.axis=0.7)
      }
    }
    if(is.null(power.lab)){
      mtext(paste0("Power for TE = ",100*power.TE,"% (x 100)"), side=3, las=0, line=1.8, cex=mycex2)
    } else {
      mtext(power.lab, side=3, las=0, line=1.8, cex=mycex2)
    }
  }
}



#' Plot Characteristics of the Estimated Distribution of the Treatment Arm-Specific Number of Endpoints
#'
#' Takes the output from the \code{\link{completeTrial.byArm}} function and generates a plot describing characteristics of the estimated distribution of the treatment arm-specific number of endpoints.
#'
#' @param arm a character string matching a treatment label in the \code{arm} variable in \code{interimData} that indicates the treatment arm for which the plot will be generated
#' @param trtNames a character vector of all treatment labels listed in the same order as in \code{trtNames} in \code{\link{completeTrial.byArm}}
#' @param eventTimeFrame a time frame within which endpoints are counted, specified in weeks as \code{c(start, end)}. If \code{NULL} (default), then all endpoints are counted.
#' @param eventPPcohort a logical value. If \code{TRUE}, only endpoints in the per-protocol cohort are counted. The default value is \code{FALSE}.
#' @param eventPriorRate a numeric vector of treatment arm-specific prior mean incidence rates for the endpoint, expressed as numbers of events per person-year at risk, matching the order of treatment arms in \code{trtNames}
#' @param eventPriorWeight a numeric vector in which each value represents a weight (i.e., a separate scenario) assigned to the prior gamma distribution of the treatment arm-specific event rate at the time when 50\% of the estimated person-time at risk in the given \code{arm} has been accumulated
#' @param xlim a numeric vector of the form \code{c(xmin, xmax)} for the user-specified x-axis limits. If \code{NULL} (default), then the computed range of x-axis values will be used.
#' @param xlab a character string for the user-specified x-axis label. If \code{NULL} (default), then the label "Number of Infections in Group \code{arm} (n)" will be used.
#' @param ylab a character string for the user-specified y-axis label. If \code{NULL} (default), then the label "P(Number of Infections in Group \code{arm} >= n ) x 100" will be used.
#' @param fileDir a character string specifying a path for the input directory
#'
#' @return None. The function is called solely for plot generation.
#'
#' @examples
#' arm <- rep(c("C3","T1","T2"), each=250)
#' schedule <- rbinom(length(arm), 1, 0.01)
#' entry <- rpois(length(arm), lambda=60)
#' entry <- entry - min(entry)
#' last_visit_dt <- entry + runif(length(arm), min=0, max=80)
#' event <- rbinom(length(arm), 1, 0.01)
#' dropout <- rbinom(length(arm), 1, 0.02)
#' dropout[event==1] <- 0
#' exit <- rep(NA, length(arm))
#' exit[event==1] <- last_visit_dt[event==1] + 5
#' exit[dropout==1] <- last_visit_dt[dropout==1] + 5
#' followup <- ifelse(event==1 | dropout==1, 0, 1)
#' interimData <- data.frame(arm=arm, schedule2=schedule, entry=entry, exit=exit,
#' last_visit_dt=last_visit_dt, event=event, dropout=dropout, complete=0, followup=followup)
#'
#' weights <- c(0.2, 0.4, 0.6)
#' for (j in 1:length(weights)){
#'   completeTrial.byArm(interimData=interimData, nTrials=50,
#'   trtNames=c("C3","T1","T2"),N=c(500,500,500),
#'   enrollRatePeriod=24, eventPriorWeight=weights[j], eventPriorRate=c(0.06,0.03,0.03),
#'   fuTime=80, visitSchedule=seq(0, 80, by=4), visitSchedule2=c(0,seq(from=8,to=80,by=12)),
#'   saveDir="./", randomSeed=9)
#' }
#'
#' pdf(file=paste0("./","rcdf_byArm_arm=T1_",
#' "eventPriorRateC3=0.06_eventPriorRateT1=0.03_eventPriorRateT2=0.03.pdf"), width=6,
#' height=5)
#' plotRCDF.byArm(arm="T1", trtNames=c("C3","T1","T2"), eventPriorRate=c(0.06,0.03,0.03),
#' eventPriorWeight=weights, fileDir="./")
#' dev.off()
#'
#' @seealso \code{\link{completeTrial.byArm}} and \code{\link{plotRCDF.pooledArms}}
#' @export
plotRCDF.byArm<-function(arm,
                         trtNames,
                         eventTimeFrame=NULL,
                         eventPPcohort=FALSE,
                         eventPriorRate,
                         eventPriorWeight,
                         xlim=NULL,
                         xlab=NULL,
                         ylab=NULL,
                         fileDir){
  if (!(arm %in% trtNames)){ stop("'arm' is not included in 'trtNames'.") }

  #plot arm
  dat  <- vector("list", length(eventPriorWeight))

  for (j in 1:length(eventPriorWeight)){
    wt<-eventPriorWeight[j]

    fileName <- paste0("completeTrial_byArm_", paste(paste0("eventPriorRate", trtNames, "=", eventPriorRate), collapse="_"), "_eventPriorWt=", wt, ".RData")
    if (file.exists(file.path(fileDir, fileName))){ load(file.path(fileDir, fileName)) }

    legend.Prior.weight<-trialObj$BetaOverBetaPlusTk[which(trtNames==arm)]
    TNI<-rep(NA,length(trialObj$trialData))

    if (is.null(eventTimeFrame)){
      if (!eventPPcohort){
        for (i in 1:length(TNI)){
          TNI[i] <- sum(subset(trialObj$trialData[[i]], arm==arm)$event)
        }
      } else {
        if (!"pp" %in% colnames(trialObj$trialData[[1]])){
          stop("trialData does not have a variable named 'pp'")
        }
        for (i in 1:length(TNI)){
          TNI[i] <- sum(subset(trialObj$trialData[[i]], arm==arm & pp==1)$event)
        }
      }
    } else {
      if (!eventPPcohort){
        for (i in 1:length(TNI)){
          trialObj$trialData[[i]]$eventTime <- trialObj$trialData[[i]]$exit - trialObj$trialData[[i]]$entry
          TNI[i] <- sum(subset(trialObj$trialData[[i]], eventTime >= eventTimeFrame[1] & eventTime <= eventTimeFrame[2] & arm==arm)$event)
        }
      } else {
        if (!"pp" %in% colnames(trialObj$trialData[[1]])){
          stop("trialData does not have a variable named 'pp'")
        }
        for (i in 1:length(TNI)){
          trialObj$trialData[[i]]$eventTime <- trialObj$trialData[[i]]$exit - trialObj$trialData[[i]]$entry
          TNI[i] <- sum(subset(trialObj$trialData[[i]], eventTime >= eventTimeFrame[1] & eventTime <= eventTimeFrame[2] & arm==arm & pp==1)$event)
        }
      }
    }

    rm(trialObj)
    mu.TNI<-round(mean(TNI),0)
    CI.TNI<-round(c(mean(TNI) - qnorm(0.975)*sd(TNI)/sqrt(length(TNI)),
                    mean(TNI) + qnorm(0.975)*sd(TNI)/sqrt(length(TNI))),0)


    # get a data frame of x- and y-coordinates of jumps in the RCDF step function
    TNI_rcdf<-to_rcdf(TNI,min=min(TNI))

    dat[[j]]<-list(TNI_rcdf=TNI_rcdf, mu.TNI=mu.TNI, CI.TNI=CI.TNI,legend.Prior.weight=legend.Prior.weight)
  }

  var.all<-dat[[1]]$TNI_rcdf$var
  if(length(dat)>1){
    for(i in 2:length(dat)){
      var.all<-c(var.all, dat[[i]]$TNI_rcdf$var)
    }
  }

  if(is.null(xlim)){
    myxlim<-c(min(var.all),max(var.all))
  }else{
    myxlim<-xlim
  }

  #if(myxlim[2]>=70){
  #  tmp<-seq(20,200,20)
  #}else{
  tmp<-seq(10,200,10)
  #}
  tmp<-tmp[tmp>=(myxlim[1]+5) & tmp<=(myxlim[2]-5)]

  x.label<-c(myxlim[1],tmp,myxlim[2])
  colors<-c("blue","red","seagreen")
  pchar <- 15:17

  par(mar=c(2.5,2.8,2.5,0.5),oma=c(1,1,1,1),las=1)
  mycex<-0.75
  mycex2<-1
  plot(myxlim, c(0,1), type="n", main="", axes=FALSE,xlab="", ylab="")

  #mean of n
  polygon(x=c(rep(min(x.label),2),rep(max(x.label),2)), y=c(0,-2,-2,0), border=NA, col="gray90")
  axis(side=1, at=min(x.label), #+0.01*(max(x.label)-min(x.label))
       labels="Mean:", tick=FALSE, line=-2, cex.axis=0.6, hadj=0)

  tmp<-c(dat[[1]]$mu.TNI,dat[[2]]$mu.TNI,dat[[3]]$mu.TNI)

  widthchange<-length(tmp)!=length(unique(tmp))
  widthloc<-NULL
  for(i in 2:length(tmp)){
    if(tmp[i] %in% tmp[1:i-1]) widthloc<-i
  }

  for(j in 1:length(dat)){
    mylwd<-1.5
    if(!is.null(widthloc)){if(j==widthloc) mylwd<-0.7}
    segments(x0=dat[[j]]$mu.TNI, y0=0, y1=1, col=colors[[j]], lwd=mylwd)
  }

  tmp <- unique(tmp)

  # suppressWarnings() is used because if 'tmp' is of length 1, then diff() returns an empty vector, and 'min' issues a warning
  if (suppressWarnings(min(diff(sort(tmp))))>=3 | max(tmp)-min(tmp)==0){
    axis(side=1, at=tmp, labels=tmp, tick=FALSE, line=-2, cex.axis=0.5)
  } else {
    if (length(tmp)==3){
      # at least one of the differences is <3, i.e., 'idx' is not empty
      idx <- which(diff(sort(tmp))<3)
      if (length(idx)==1){
        if (idx==1){
          atValue <- sort(tmp) + c(-0.5,0,0)
        }
        if (idx==2){
          atValue <- sort(tmp) + c(0,0,0.5)
        }
      } else {
        atValue <- sort(tmp) + c(-0.4,0,0.4)
      }
      for (k in 1:length(tmp)){
        axis(side=1, at= atValue[k], labels=sort(tmp)[k], tick=FALSE, line=-2, cex.axis=0.5)
      }
    } else {
      atValue <- sort(tmp) + c(-0.25,0.25)
      for (k in 1:length(tmp)){
        axis(side=1, at= atValue[k], labels=sort(tmp)[k], tick=FALSE, line=-2, cex.axis=0.5)
      }
    }
  }

  for(i in 1:length(dat)){
    lines(c(min(x.label),dat[[i]]$TNI_rcdf$var),c(1,dat[[i]]$TNI_rcdf$RCDF), type="S",lty=1,col=colors[i],lwd=1.2)
    points(dat[[i]]$TNI_rcdf$var,dat[[i]]$TNI_rcdf$RCDF,lty=1, pch=pchar[i], col=colors[i], cex=ifelse(i==1,0.4,0.5))
  }

  axis(side=1, at=x.label, labels=F,cex.axis=mycex)
  axis(side=1, at=x.label, labels=x.label,tick=F, line=-0.6,cex.axis=mycex)
  yTicks <- seq(0,1,by=0.2)
  axis(side=2, at=yTicks,labels=F, cex.axis=mycex)
  axis(side=2, at=yTicks,labels=100*yTicks, tick=F, cex.axis=mycex)

  if(is.null(xlab)){mtext(paste0("Number of Infections in Group ",arm," (n)"), side=1, las=0, line=2, cex=mycex2)
  }else{mtext(xlab, side=1, las=0, line=2, cex=mycex2)}
  if(is.null(ylab)){mtext(paste0("P( Number of Infections in Group ",arm," >= n ) x 100"), side=2, las=0, line=2.5, cex=mycex2)
  }else{mtext(ylab, side=2, las=0, line=2.5, cex=mycex2)}
  legend(0.75*max(x.label),0.9,legend=round(sapply(dat, "[[", "legend.Prior.weight"), 2), cex=0.7, col=colors, pch=pchar, lty=1, bty = "n",
         title="Prior weight")
}
