
####################
# 1.  Load libraries
####################

library(tidyverse)
library(lubridate)

library(png)
library(ggalluvial)
library(RColorBrewer)
library(cowplot)
library(flextable)


library(survival)
library(survminer)

library(grid)
library(gridExtra)
library(lattice)

# library(rlang)
# library(magrittr)



####################
# 2.  Load functions
####################


# function takes kaplan-meier fit objects (output from survfit) and
# gives coverage for given time
# example for fit below:
# fit <- survfit(Surv(V1.SA.time,V1.SA.censor) ~ 1,
#                type = "kaplan-meier",
#                error = "greenwood",
#                conf.type = "log-log",
#                data = analysis.DTP)

getCoverage <- function(KMfit, t){
        
        #objects without strata
        if(is.null(KMfit$strata)){
                
                covDf <- data.frame(
                        time = round(KMfit$time,1),
                        cov = round(1- KMfit$surv,2),
                        lb =round(1-KMfit$upper,2),
                        ub =round(1-KMfit$lower,2)
                )  %>%
                        mutate(diff = abs(time - t ),
                               strata = "Overall") %>%
                        filter(
                                diff==min(diff)
                        ) %>% select(-diff) %>%
                        filter(!duplicated(strata)) #if two points are equidistant
                
                #objects with strata
        } else{
                s <- NULL
                
                for (i in 1:length(KMfit$strata)){
                        n <- names(KMfit$strata)[i]
                        s <- c(s,
                               rep(n,KMfit$strata[i]) )
                        
                }
                
                
                covDf <- data.frame(
                        time = round(KMfit$time,1),
                        cov = round(1- KMfit$surv,2),
                        lb =round(1-KMfit$upper,2),
                        ub =round(1-KMfit$lower,2),
                        strata = s
                )  %>%
                        mutate(diff = abs(time - t )) %>%
                        group_by(strata) %>%
                        filter(
                                diff==min(diff)
                        ) %>% select(-diff)%>%
                        filter(!duplicated(strata))#if two points are equidistant
                
                
        }
        
        return(covDf)
        
}


getSurvivalTime <- function(KMfit, quant){
        
        #objects without strata
        if(is.null(KMfit$strata)){
                
                quantDf <- data.frame(
                        time = round(KMfit$time,1),
                        cov = round(1- KMfit$surv,2)
                )  %>%
                        mutate(diff = abs(cov - quant ),
                               strata = "Overall") %>%
                        filter(
                                diff==min(diff)
                        ) %>% select(-diff) %>%
                        filter(!duplicated(strata)) #if two points are equidistant
                
                #objects with strata
        } else{
                s <- NULL
                
                for (i in 1:length(KMfit$strata)){
                        n <- names(KMfit$strata)[i]
                        s <- c(s,
                               rep(n,KMfit$strata[i]) )
                        
                }
                
                
                quantDf <- data.frame(
                        time = round(KMfit$time,1),
                        cov = round(1- KMfit$surv,2),
                        strata = s
                )  %>%
                        mutate(diff = abs(cov - quant)) %>%
                        group_by(strata) %>%
                        filter(
                                diff==min(diff)
                        ) %>% select(-diff)%>%
                        filter(!duplicated(strata))#if two points are equidistant
                
                
        }
        
        return(quantDf)
        
}



#calculates 1-year coverage and the proportion that are given early or ontime

KMestim <- function(data){
        
        #coverage
        
        d1 <- data %>%
                survfit(data = .,
                        Surv(V1.SA.time,V1.SA.censor) ~ 1,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "First")
        
        
        
        d2 <- data %>%
                survfit(data = .,
                        Surv(V2.SA.time2,V2.SA.censor2) ~ 1,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Second")
        
        d3 <- data %>%
                survfit(data = .,
                        Surv(V3.SA.time2,V3.SA.censor2) ~ 1,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Third")
        
        d4<- data %>%
                survfit(data = .,
                        Surv(V1.SA.time,V1.SA.censor) ~ District,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "First")
        
        
        
        d5<- data %>%
                survfit(data = .,
                        Surv(V2.SA.time2,V2.SA.censor2) ~ District,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Second")
        
        d6 <- data %>%
                survfit(data = .,
                        Surv(V3.SA.time2,V3.SA.censor2) ~ District,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Third")
        
        d7<- data %>%
                survfit(data = .,
                        Surv(V1.SA.time,V1.SA.censor) ~ birth.month,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "First")
        
        
        
        d8<- data %>%
                survfit(data = .,
                        Surv(V2.SA.time2,V2.SA.censor2) ~ birth.month,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Second")
        
        d9 <- data %>%
                survfit(data = .,
                        Surv(V3.SA.time2,V3.SA.censor2) ~ birth.month,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Third")
        
        
        
        # Lateness part1
        l1 <- data %>%
                survfit(data = .,
                        Surv(V1.SA.time,V1.SA.censor) ~ 1,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(10) %>% mutate(dose = "First")
        
        l2 <- data %>%
                survfit(data = .,
                        Surv(V2.SA.time,V2.SA.censor) ~ 1,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(8) %>% mutate(dose = "Second")
        
        l3 <- data %>%
                survfit(data = .,
                        Surv(V3.SA.time,V3.SA.censor) ~ 1,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(8) %>% mutate(dose = "Third")
        
        
        l4 <- data %>%
                survfit(data = .,
                        Surv(V1.SA.time,V1.SA.censor) ~ District,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(10) %>% mutate(dose = "First")
        
        l5 <- data %>%
                survfit(data = .,
                        Surv(V2.SA.time,V2.SA.censor) ~ District,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(8)%>% mutate(dose = "Second")
        
        l6 <- data %>%
                survfit(data = .,
                        Surv(V3.SA.time,V3.SA.censor) ~ District,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(8) %>% mutate(dose = "Third")
        
        l7 <- data %>%
                survfit(data = .,
                        Surv(V1.SA.time,V1.SA.censor) ~ birth.month,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(10) %>% mutate(dose = "First")
        
        l8 <- data %>%
                survfit(data = .,
                        Surv(V2.SA.time,V2.SA.censor) ~ birth.month,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(8)%>% mutate(dose = "Second")
        
        l9 <- data %>%
                survfit(data = .,
                        Surv(V3.SA.time,V3.SA.censor) ~ birth.month,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(8) %>% mutate(dose = "Third")
        
        # Lateness part2
        q1 <- data %>%
                survfit(data = .,
                        Surv(V1.SA.time,V1.SA.censor) ~ 1,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "First")
        
        q2 <- data %>%
                survfit(data = .,
                        Surv(V2.SA.time,V2.SA.censor) ~ 1,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Second")
        
        q3 <- data %>%
                survfit(data = .,
                        Surv(V3.SA.time,V3.SA.censor) ~ 1,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Third")
        
        
        q4 <- data %>%
                survfit(data = .,
                        Surv(V1.SA.time,V1.SA.censor) ~ District,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "First")
        
        q5 <- data %>%
                survfit(data = .,
                        Surv(V2.SA.time,V2.SA.censor) ~ District,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52)%>% mutate(dose = "Second")
        
        q6 <- data %>%
                survfit(data = .,
                        Surv(V3.SA.time,V3.SA.censor) ~ District,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Third")
        
        q7 <- data %>%
                survfit(data = .,
                        Surv(V1.SA.time,V1.SA.censor) ~ birth.month,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "First")
        
        q8 <- data %>%
                survfit(data = .,
                        Surv(V2.SA.time,V2.SA.censor) ~ birth.month,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52)%>% mutate(dose = "Second")
        
        q9 <- data %>%
                survfit(data = .,
                        Surv(V3.SA.time,V3.SA.censor) ~ birth.month,
                        type = "kaplan-meier",
                        error = "greenwood",
                        conf.type = "log-log") %>%
                getCoverage(52) %>% mutate(dose = "Third")
        
        
        
        d <- bind_rows(d1,d2,d3,d4,d5,d6,d7,d8,d9) %>%
                mutate(type="coverage")
        l <- bind_rows(l1,l2,l3,l4,l5,l6,l7,l8,l9) %>%
                mutate(type="notLate")
        q <- bind_rows(q1,q2,q3,q4,q5,q6,q7,q8,q9) %>%
                mutate(type="cov|previous")
        
        return( bind_rows(d,l,q))
        
}


#get value of survival from KM curve
getSurv <- function(fit,t=52){
        
        df <-data.frame(time=fit$time,
                        surv=fit$surv,
                        lb=fit$lower,
                        ub=fit$upper) %>%
                mutate(diff=abs(t-time))%>%
                filter(diff==min(diff) )%>%
                filter(time==max(time) )%>%
                select(-diff,-time)
        
        df
        
}

#get weights from last KM curve
getWeights <- function(data, #bring in the survival data
                       fit, # last survival fit
                       dose # 1 2 or 3 for the last fit
){
        
        # get the pieces of the fit out that are important
        l <- length(fit$surv)
        times <- fit$time
        events <- fit$n.event
        survival <- fit$surv
        jumps <- c(1,survival[1:(l-1)])-survival
        atrisk <- fit$n.risk[1]
        
        #Figure out which dose this is
        Tstring <- paste0("T",dose)
        Dstring <- paste0("D",dose)
        
        # set up the dataframe to include new weights
        data$newweight <- NA
        
        
        #loop through each time
        for(i in 1:l){
                if(events[i]>0){
                        
                        t <- round(times[i],1)
                        w <-jumps[i]*atrisk/events[i]
                        
                        data$newweight[data[Tstring]==t & data[Dstring]==1]<-w
                        
                        
                }
                
                
        }
        
        #multiply this output by all previous weights
        weightedVector <- data$newweight
        
        
        return(weightedVector)
}


#input: data from individuals of particular population

getCovLat<- function(data){
        #get dataset ready
        analysis <- data %>%
                mutate(
                        B1=B1.SA.time,
                        B2=B2.SA.time,
                        B3=B3.SA.time
                ) %>%
                select(ID,
                       T1=T1.SA.time,
                       D1=V1.SA.censor,
                       T2=T2.SA.time,
                       D2=V2.SA.censor,
                       T3=T3.SA.time,
                       D3=V3.SA.censor,
                       B1,
                       B2,
                       B3
                ) %>%
                mutate(
                        T1=round(T1,1),
                        T2=round(T2,1),
                        T3=round(T3,1),
                        
                        B1=round(B1,1),
                        B2=round(B2,1),
                        B3=round(B3,1)
                )%>%
                mutate(W1=1)
        
        #fit the models and weights
        #first dose
        fit1 <- survfit(Surv(T1,D1) ~ 1, data = analysis)
        #second dose
        analysis$W2 <- analysis$W1 * getWeights(analysis,fit1,1)
        fit2 <- survfit(Surv(T2, D2) ~ 1, data = analysis,weights=W2)
        Lfit2 <- survfit(Surv(B2, D2) ~ 1, data = analysis,weights=W2)
        # third dose
        analysis$W3 <- analysis$W2 * getWeights(analysis,fit2,2)
        fit3 <- survfit(Surv(T3, D3) ~ 1, data = analysis,weights=W3)
        Lfit3 <- survfit(Surv(B3, D3) ~ 1, data = analysis,weights=W3)
        
        
        
        
        
        #prepare output
        covdf <- bind_rows(
                (1-getSurv(fit1,52)) %>%
                        mutate(dose="First"),
                ((1-getSurv(fit2,52))*(1-getSurv(fit1,52))$surv) %>%
                        mutate(dose="Second"),
                ((1-getSurv(fit3,52))*(1-getSurv(fit2,52))$surv*(1-getSurv(fit1,52))$surv) %>%
                        mutate(dose="Third")
        ) %>%
                rename(cov=surv) %>%
                #switch upperbound and lowerbound
                rename(tmp=lb,
                       lb=ub) %>%
                rename(ub=tmp)
        
        
        latdf <- data.frame(
                lateness=c(
                        getSurv(fit1,10)$surv-getSurv(fit1,52)$surv,
                        (getSurv(Lfit2,8)$surv-getSurv(fit2,52)$surv)*(1-getSurv(fit1,52)$surv),
                        (getSurv(Lfit3,8)$surv-getSurv(fit3,52)$surv)*(1-getSurv(fit2,52)$surv)*(1-getSurv(fit1,52)$surv)
                ),
                dose= c("First","Second","Third")
        )
        
        
        #coverage and lateness
        return(list(covdf,latdf))
}


stratifyCovLat <- function(data,
                           column #column of interest to stratify
) {
        
        
        lvlies <- unlist(unique(data[column]))
        
        covdf <- NULL
        latdf <- NULL
        
        for (i in lvlies){
                
                # print(i)
                filt <- as.vector(data[column] == i)
                fdata <- data[filt,]
                
                tmp <- getCovLat(fdata)
                
                covdf <- bind_rows(
                        covdf,
                        tmp[[1]]%>% mutate(stratum=i)
                        
                )
                
                latdf <- bind_rows(
                        latdf,
                        tmp[[2]] %>% mutate(stratum=i)
                )
        }
        
        
        return(list(covdf,latdf))
        
}

weightedKMcurves<- function(data){
        analysis <- data %>% 
                mutate(
                        B1=B1.SA.time,
                        B2=B2.SA.time,
                        B3=B3.SA.time
                ) %>%
                select(ID,
                       T1=T1.SA.time,
                       D1=V1.SA.censor,
                       T2=T2.SA.time,
                       D2=V2.SA.censor,
                       T3=T3.SA.time,
                       D3=V3.SA.censor,
                       B1,
                       B2,
                       B3
                ) %>%
                mutate(
                        T1=round(T1,1),
                        T2=round(T2,1),
                        T3=round(T3,1),
                        
                        B1=round(B1,1),
                        B2=round(B2,1),
                        B3=round(B3,1)
                )%>%
                mutate(W1=1)
        
        
        
        
        #first dose
        fit1 <- survfit(Surv(T1,D1) ~ 1, data = analysis)
        #second dose
        analysis$W2 <- analysis$W1 * getWeights(analysis,fit1,1)
        fit2 <- survfit(Surv(T2, D2) ~ 1, data = analysis,weights=W2)
        #third dose
        analysis$W3 <- analysis$W2 * getWeights(analysis,fit2,2)
        fit3 <- survfit(Surv(T3, D3) ~ 1, data = analysis,weights=W3)
        
        plot_df <-bind_rows(
                data.frame(
                        time=fit1$time,
                        vaccinated=(1-fit1$surv)*sum(analysis$W1,na.rm = TRUE),
                        upper=(1-fit1$lower)*sum(analysis$W1,na.rm = TRUE),
                        lower=(1-fit1$upper)*sum(analysis$W1,na.rm = TRUE),
                        dose="First"
                ),
                data.frame(
                        time=fit2$time,
                        vaccinated=(1-fit2$surv)*sum(analysis$W2,na.rm = TRUE),
                        upper=(1-fit2$lower)*sum(analysis$W2,na.rm = TRUE),
                        lower=(1-fit2$upper)*sum(analysis$W2,na.rm = TRUE),
                        dose="Second"
                ),
                data.frame(
                        time=fit3$time,
                        vaccinated=(1-fit3$surv)*sum(analysis$W3,na.rm = TRUE),
                        upper=(1-fit3$lower)*sum(analysis$W3,na.rm = TRUE),
                        lower=(1-fit3$upper)*sum(analysis$W3,na.rm = TRUE),
                        dose="Third"
                )
        )
        
        return(plot_df)
}




#used for figure 4 and for getting the median time to vacination
createWeightedDF<- function(data){
        
        analysis <- data %>%
                mutate(
                        B1=B1.SA.time,
                        B2=B2.SA.time,
                        B3=B3.SA.time,
                        
                        T1=T1.SA.time,
                        D1=V1.SA.censor,
                        T2=T2.SA.time,
                        D2=V2.SA.censor,
                        T3=T3.SA.time,
                        D3=V3.SA.censor
                ) %>%
                mutate(
                        T1=round(T1,1),
                        T2=round(T2,1),
                        T3=round(T3,1),
                        
                        B1=round(B1,1),
                        B2=round(B2,1),
                        B3=round(B3,1)
                )%>%
                mutate(W1=1)
        
        
        #first dose
        fit1 <- survfit(Surv(T1,D1) ~ 1, data = analysis)
        #second dose
        analysis$W2 <- analysis$W1 * getWeights(analysis,fit1,1)
        fit2 <- survfit(Surv(T2, D2) ~ 1, data = analysis,weights=W2)
        # third dose
        analysis$W3 <- analysis$W2 * getWeights(analysis,fit2,2)
        
        return(analysis)
        
}


#Creates river plot using ggalluvial for one vaccine
makeRiverPlot <- function(data , antigen) {
        
        #filter data to only children who were observed for one year
        data <- data %>% filter(full_observe)
        # data <- data %>% filter(DDN < (last.obs-365) )
        
        
        # Limit down to columns used
        ifelse(antigen != "Rota",
               allu <- data %>% select(ID, V1.timeliness, V2.timeliness ,
                                       V3.timeliness),
               allu <- data %>% select(ID, V1.timeliness, V2.timeliness))
        
        
        
        
        # switch to long dataset and rename categories
        allu2<- gather(allu, vacc, Timeliness, -ID) %>%
                mutate(Timeliness = ifelse(is.na(Timeliness),"Never",Timeliness),
                       vacc = as.factor(vacc)) %>%
                mutate(Timeliness = factor(Timeliness,
                                           levels= c("Early", "OnTime", "Late", "Never")))
        levels(allu2$vacc)<- c("First","Second" , "Third")
        
        
        
        #allu2$Timeliness[is.na(allu2$Timeliness)==TRUE] <-"Never"
        #allu2$Timeliness[allu2$Timeliness=="Early"|allu2$Timeliness=="OnTime"] <-"Early/On Time"
        #allu2$vacc <- as.factor(allu2$vacc)
        
        # add counting variable
        allu2$freq <- 1 
        
        #create label data frame
        obs <- dim(data)[1]
        
        ldf <- allu2 %>% group_by(vacc, Timeliness) %>%
                summarize(freq=sum(freq))
        
        lbs <-                ldf %>%
                spread(Timeliness, freq) %>%
                mutate(Early =Never+ Late + OnTime + Early/2,
                       OnTime = Never + Late + OnTime/2,
                       Late = Never + Late/2,
                       Never = Never/2)%>%
                # mutate(`Early/On Time` = `Early/On Time` + Late/2,
                #        Late = Never + Late/2,
                #        Never = Never/2) %>%
                gather(Timeliness, ypos, -vacc) %>%
                left_join(ldf, by=c("Timeliness", "vacc")) %>%
                mutate(percent = round(100*freq/obs)) %>%
                filter(percent >= 7) %>%
                mutate(percent = paste(percent,"%",sep=""))
        
        
        # print(as.data.frame(ldf))
        # print(as.data.frame(lbs))
        # #create plot
        rp <- ggplot(allu2,
                     aes(x = vacc, stratum = Timeliness, alluvium = ID,
                         y = freq,
                         fill = Timeliness ,
                         #label = freq
                     )
        )+
                scale_x_discrete(expand = c(.1, .1)) +
                scale_fill_manual(values=c("blue","#4daf4a","#ff7f00","#e41a1c")) +
                geom_flow() +
                geom_stratum(alpha = .5)+
                #geom_text(stat = "stratum") +
                geom_text(inherit.aes = FALSE,data=lbs, size = 4,
                          aes(x=vacc,y=ypos,
                              label=percent))+
                ylab("Vaccination Cards") +
                xlab("Dose") + ggtitle(antigen)
        
        return(rp)
        
        
        
}



# Overall Goal
### Break dataset up by calendar week and timeliness
### Calculate the hazard at the weekly level

splitter <- function(dat) {
        
        cutPoints <- floor_date(seq.Date(as.Date("2015-01-01"),as.Date("2017-02-10"),by=7),"weeks")+1
        
        dat %>% survSplit(data=.,Surv(as.numeric(entry), as.numeric(exit), censor)~.,
                          cut = as.numeric(cutPoints),
                          end = "exit") %>%
                mutate(PT=exit-tstart,
                       tstart=as.Date(tstart, origin="1970-01-01"))%>%
                mutate(week=cut(tstart,cutPoints))
        
}




fitter <-function(data){
        
        #fit poisson model and predict values
        model <- data %>% glm(data=., censor~ -1 + factor(week)+District,
                              offset=log(PT),
                              family=poisson(link="log")
        ) %>% predict.glm(data.frame(week=sort(unique(data$week)),
                                     District="TULEAR I",PT=7),
                          se.fit=TRUE)
        
        conFint <- data.frame(fit=model$fit,
                              lwr=model$fit-1.96*model$se.fit,
                              upr=model$fit+1.96*model$se.fit
        )
        
        # interval='confidence', level=0.95) 
        # output predictions
        data.frame(week=sort(unique(data$week)),
                   Dose=unique(data$Dose),
                   timeliness=unique(data$timeliness),
                   exp(conFint)*100)
        
        
}


#runs simulations and tracks doses/coverage
VWsimulation <- function(catchup=1){
        
        #load rates of vaccination
        rates <- read_csv("vaccRateMatrix.csv") %>% select(-X1)
        RegRates <- as.matrix(rates)[,2:4]
        
        #define rates of those who are "caught up"
        vwRates <- matrix(0,nrow=52,ncol=2)
        vwRates[1:47,] <- RegRates[c(1:3,9:52),c(2,3)]
        colnames(vwRates) <- c("VWrate2", "VWrate3")
        
        #define which birth weeks
        birthweeks<- seq.Date(as.Date("2014-01-01"),as.Date("2017-03-01"),by=7)
        VWs <- as.Date(c("2014-05-14","2014-10-22",
                         "2015-05-13","2015-10-21",
                         "2016-05-11","2016-10-19"
        ))
        # VWs <- as.Date(c("2020-05-14"))
        
        #tracking
        covTracker  <- NULL
        doseTracker <- NULL
        
        
        #go through each cohort
        for (i in birthweeks){
                #create compartments to fill in
                #Start with no vaccination
                #1000 people in a cohort 
                C.matrix <- matrix(NA,nrow=53,ncol=4) #followup in weeks and doses
                C.matrix[1,] <-0
                C.matrix[1,1] <- 1000
                
                #Compartments for those caught up by VW
                VW.matrix <- matrix(NA,nrow=53,ncol=3) #followup in weeks and doses
                VW.matrix[1,] <-0
                
                #go through each week of the year
                for (j in 1:52){
                        date <- as.Date( i + j*7,origin="1970-01-01")
                        
                        if(date %in% VWs ){ # VACCINATION WEEK
                                
                                #vaccinations
                                new.matrix <- matrix(NA, nrow=52, ncol = 3)
                                
                                #keep rates for those getting vaccinated on-time or early
                                new.matrix[1:10,1] <- C.matrix[1:10,1] * RegRates[1:10,1]
                                new.matrix[1:8,2] <- C.matrix[1:8,2] * RegRates[1:8,2]
                                new.matrix[1:8,3] <- C.matrix[1:8,3] * RegRates[1:8,3]
                                
                                #catch up everyone else at specified rate
                                new.matrix[11:52,1] <- C.matrix[11:52,1] * (RegRates[11:52,1] + catchup*(1-RegRates[11:52,1]))
                                new.matrix[9:52,2] <- C.matrix[9:52,2]   * (RegRates[9:52,2] + catchup*(1-RegRates[9:52,2]))
                                new.matrix[9:52,3] <- C.matrix[9:52,3]   * (RegRates[9:52,3] + catchup*(1-RegRates[9:52,3]))
                                # everyone left in VW.matrix will be late and eligible for VW
                                VW.new.matrix <- VW.matrix[-53,1:2] * (catchup + vwRates*(1-catchup))
                                
                                #follow doses
                                normaldoses <- colSums(C.matrix[-53,1:3] * RegRates,na.rm = TRUE)
                                doses[2:3] <- doses[2:3] + colSums(VW.matrix[-53,1:2] * vwRates, na.rm = TRUE)
                                
                                VWdoses <- c(    sum(C.matrix[11:52,1] * catchup*(1-RegRates[11:52,1]),na.rm = TRUE),
                                                 sum(C.matrix[9:52,2]   * RegRates[9:52,2] + catchup*(1-RegRates[9:52,2]),na.rm = TRUE),
                                                 sum(C.matrix[9:52,3]   * RegRates[9:52,3] + catchup*(1-RegRates[9:52,3]),na.rm = TRUE)) +
                                        c(0,colSums(VW.new.matrix,na.rm = TRUE))
                                
                                # #follow doses
                                # 
                                # 
                                # 
                                # normaldoses <- c(    sum(new.matrix[1:10,1],na.rm = TRUE),
                                #                      sum(new.matrix[1:8,2],na.rm = TRUE),
                                #                      sum(new.matrix[1:8,3],na.rm = TRUE))
                                # 
                                # VWdoses <- c(    sum(new.matrix[11:52,1],na.rm = TRUE),
                                #                  sum(new.matrix[9:52,2],na.rm = TRUE),
                                #                  sum(new.matrix[9:52,3],na.rm = TRUE)) +
                                #         c(0,colSums(VW.new.matrix,na.rm = TRUE))
                                # 
                                
                                
                                #temporary matrices to fill in
                                tempC.matrix <- matrix(NA,nrow=53,ncol=4) #followup in weeks and doses
                                tempVW.matrix <- matrix(NA,nrow=53,ncol=3) #followup in weeks and doses
                                
                                #first row
                                tempC.matrix[1,] <- c(0,normaldoses)
                                tempVW.matrix[1,] <- VWdoses
                                
                                #the rest (let people age)
                                tempC.matrix[2:53,1:3] <- C.matrix[1:52,1:3]-new.matrix
                                tempC.matrix[2:53,4]   <- C.matrix[1:52,4]
                                
                                tempVW.matrix[2:53,1:2] <- VW.matrix[1:52,1:2]-VW.new.matrix
                                tempVW.matrix[2:53,3]   <- VW.matrix[1:52,3]
                                
                                # bring it back to the compartment matrix
                                C.matrix <- tempC.matrix
                                VW.matrix <- tempVW.matrix
                                
                                
                                
                                # tracking
                                doses <-  normaldoses + VWdoses
                                
                                
                                
                        } else{ # NOT VACCINATION WEEKS
                                
                                #vaccinations
                                new.matrix <- C.matrix[-53,1:3] * RegRates
                                VW.new.matrix <- VW.matrix[-53,1:2] * vwRates
                                
                                #follow doses
                                doses <- colSums(new.matrix,na.rm = TRUE)
                                doses[2:3] <- doses[2:3] + colSums(VW.new.matrix,na.rm = TRUE)
                                
                                #temporary matrices to fill in
                                tempC.matrix <- matrix(NA,nrow=53,ncol=4) #followup in weeks and doses
                                tempVW.matrix <- matrix(NA,nrow=53,ncol=3) #followup in weeks and doses
                                
                                #first row
                                tempC.matrix[1,] <- c(0,doses)
                                tempVW.matrix[1,] <- c(0,0,0)
                                
                                #the rest (let people age)
                                tempC.matrix[2:53,1:3] <- C.matrix[1:52,1:3]-new.matrix
                                tempC.matrix[2:53,4]   <- C.matrix[1:52,4]
                                
                                tempVW.matrix[2:53,1:2] <- VW.matrix[1:52,1:2]-VW.new.matrix
                                tempVW.matrix[2:53,3]   <- VW.matrix[1:52,3]
                                
                                # bring it back to the compartment matrix
                                C.matrix <- tempC.matrix
                                VW.matrix <- tempVW.matrix
                        }
                        
                        doseage <-data.frame(Catchup=catchup,
                                             birthdate=as.Date( i ,origin="1970-01-01"),
                                             date,
                                             First=doses[1],
                                             Second=doses[2],
                                             Third=doses[3])
                        
                        doseTracker<- rbind(doseTracker,doseage)
                        
                        
                }
                
                coverage <- cumsum(colSums(C.matrix) + c(0,colSums(VW.matrix)))
                coverage <- round((1000 - coverage)/10,2)
                coverage <- data.frame(Catchup=catchup,
                                       birthdate=as.Date( i ,origin="1970-01-01"),
                                       date,
                                       First=coverage[1],
                                       Second=coverage[2],
                                       Third=coverage[3])
                
                covTracker <- rbind(covTracker,coverage)
        }
        
        return(list(covTracker,doseTracker))
        
}



#for sensitivity analysis for reviwer 1

VWsimulationSens <- function(catchup=1){
        
        #load rates of vaccination
        rates <- read_csv("vaccRateMatrix.csv") %>% select(-X1)
        RegRates <- as.matrix(rates)[,2:4]
        
        #define rates of those who are "caught up"
        vwRates <-  RegRates[c(1:52),c(2,3)]
        colnames(vwRates) <- c("VWrate2", "VWrate3")
        
        #define which birth weeks
        birthweeks<- seq.Date(as.Date("2014-01-01"),as.Date("2017-03-01"),by=7)
        VWs <- as.Date(c("2014-05-14","2014-10-22",
                         "2015-05-13","2015-10-21",
                         "2016-05-11","2016-10-19"
        ))
        # VWs <- as.Date(c("2020-05-14"))
        
        #tracking
        covTracker  <- NULL
        doseTracker <- NULL
        
        
        #go through each cohort
        for (i in birthweeks){
                #create compartments to fill in
                #Start with no vaccination
                #1000 people in a cohort 
                C.matrix <- matrix(NA,nrow=53,ncol=4) #followup in weeks and doses
                C.matrix[1,] <-0
                C.matrix[1,1] <- 1000
                
                #Compartments for those caught up by VW
                VW.matrix <- matrix(NA,nrow=53,ncol=3) #followup in weeks and doses
                VW.matrix[1,] <-0
                
                #go through each week of the year
                for (j in 1:52){
                        date <- as.Date( i + j*7,origin="1970-01-01")
                        
                        if(date %in% VWs ){ # VACCINATION WEEK
                                
                                #vaccinations
                                new.matrix <- matrix(NA, nrow=52, ncol = 3)
                                
                                #keep rates for those getting vaccinated on-time or early
                                new.matrix[1:10,1] <- C.matrix[1:10,1] * RegRates[1:10,1]
                                new.matrix[1:8,2] <- C.matrix[1:8,2] * RegRates[1:8,2]
                                new.matrix[1:8,3] <- C.matrix[1:8,3] * RegRates[1:8,3]
                                
                                #catch up everyone else at specified rate
                                new.matrix[11:52,1] <- C.matrix[11:52,1] * (RegRates[11:52,1] + catchup*(1-RegRates[11:52,1]))
                                new.matrix[9:52,2] <- C.matrix[9:52,2]   * (RegRates[9:52,2] + catchup*(1-RegRates[9:52,2]))
                                new.matrix[9:52,3] <- C.matrix[9:52,3]   * (RegRates[9:52,3] + catchup*(1-RegRates[9:52,3]))
                                # everyone left in VW.matrix will be late and eligible for VW
                                VW.new.matrix <- VW.matrix[-53,1:2] * (catchup + vwRates*(1-catchup))
                                
                                #follow doses
                                normaldoses <- colSums(C.matrix[-53,1:3] * RegRates,na.rm = TRUE)
                                doses[2:3] <- doses[2:3] + colSums(VW.matrix[-53,1:2] * vwRates, na.rm = TRUE)
                                
                                VWdoses <- c(    sum(C.matrix[11:52,1] * catchup*(1-RegRates[11:52,1]),na.rm = TRUE),
                                                 sum(C.matrix[9:52,2]   * RegRates[9:52,2] + catchup*(1-RegRates[9:52,2]),na.rm = TRUE),
                                                 sum(C.matrix[9:52,3]   * RegRates[9:52,3] + catchup*(1-RegRates[9:52,3]),na.rm = TRUE)) +
                                        c(0,colSums(VW.new.matrix,na.rm = TRUE))
                                
                                # #follow doses
                                # 
                                # 
                                # 
                                # normaldoses <- c(    sum(new.matrix[1:10,1],na.rm = TRUE),
                                #                      sum(new.matrix[1:8,2],na.rm = TRUE),
                                #                      sum(new.matrix[1:8,3],na.rm = TRUE))
                                # 
                                # VWdoses <- c(    sum(new.matrix[11:52,1],na.rm = TRUE),
                                #                  sum(new.matrix[9:52,2],na.rm = TRUE),
                                #                  sum(new.matrix[9:52,3],na.rm = TRUE)) +
                                #         c(0,colSums(VW.new.matrix,na.rm = TRUE))
                                # 
                                
                                
                                #temporary matrices to fill in
                                tempC.matrix <- matrix(NA,nrow=53,ncol=4) #followup in weeks and doses
                                tempVW.matrix <- matrix(NA,nrow=53,ncol=3) #followup in weeks and doses
                                
                                #first row
                                tempC.matrix[1,] <- c(0,normaldoses)
                                tempVW.matrix[1,] <- VWdoses
                                
                                #the rest (let people age)
                                tempC.matrix[2:53,1:3] <- C.matrix[1:52,1:3]-new.matrix
                                tempC.matrix[2:53,4]   <- C.matrix[1:52,4]
                                
                                tempVW.matrix[2:53,1:2] <- VW.matrix[1:52,1:2]-VW.new.matrix
                                tempVW.matrix[2:53,3]   <- VW.matrix[1:52,3]
                                
                                # bring it back to the compartment matrix
                                C.matrix <- tempC.matrix
                                VW.matrix <- tempVW.matrix
                                
                                
                                
                                # tracking
                                doses <-  normaldoses + VWdoses
                                
                                
                                
                        } else{ # NOT VACCINATION WEEKS
                                
                                #vaccinations
                                new.matrix <- C.matrix[-53,1:3] * RegRates
                                VW.new.matrix <- VW.matrix[-53,1:2] * vwRates
                                
                                #follow doses
                                doses <- colSums(new.matrix,na.rm = TRUE)
                                doses[2:3] <- doses[2:3] + colSums(VW.new.matrix,na.rm = TRUE)
                                
                                #temporary matrices to fill in
                                tempC.matrix <- matrix(NA,nrow=53,ncol=4) #followup in weeks and doses
                                tempVW.matrix <- matrix(NA,nrow=53,ncol=3) #followup in weeks and doses
                                
                                #first row
                                tempC.matrix[1,] <- c(0,doses)
                                tempVW.matrix[1,] <- c(0,0,0)
                                
                                #the rest (let people age)
                                tempC.matrix[2:53,1:3] <- C.matrix[1:52,1:3]-new.matrix
                                tempC.matrix[2:53,4]   <- C.matrix[1:52,4]
                                
                                tempVW.matrix[2:53,1:2] <- VW.matrix[1:52,1:2]-VW.new.matrix
                                tempVW.matrix[2:53,3]   <- VW.matrix[1:52,3]
                                
                                # bring it back to the compartment matrix
                                C.matrix <- tempC.matrix
                                VW.matrix <- tempVW.matrix
                        }
                        
                        doseage <-data.frame(Catchup=catchup,
                                             birthdate=as.Date( i ,origin="1970-01-01"),
                                             date,
                                             First=doses[1],
                                             Second=doses[2],
                                             Third=doses[3])
                        
                        doseTracker<- rbind(doseTracker,doseage)
                        
                        
                }
                
                coverage <- cumsum(colSums(C.matrix) + c(0,colSums(VW.matrix)))
                coverage <- round((1000 - coverage)/10,2)
                coverage <- data.frame(Catchup=catchup,
                                       birthdate=as.Date( i ,origin="1970-01-01"),
                                       date,
                                       First=coverage[1],
                                       Second=coverage[2],
                                       Third=coverage[3])
                
                covTracker <- rbind(covTracker,coverage)
        }
        
        return(list(covTracker,doseTracker))
        
}








