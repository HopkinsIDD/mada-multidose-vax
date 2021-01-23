vcard <- read_csv("data/raw_data/vaccard_final_KM.csv")

#limit dataset to those born in 2015 and 2016
analysis1 <- vcard %>% 
        mutate(DDN = as.Date(DDN,format="%m/%d/%Y")) %>%
        filter(year(DDN) %in% c(2015,2016))



DTPissues1 <- analysis1 %>% 
                filter( (DTP1==0 & DTP2==1) |
                        (DTP1==0 & DTP3==1) |
                        (DTP2==0 & DTP3==1)
                       )


PCV10issues1 <- analysis1 %>% 
        filter( (PCV10_1==0 & PCV10_2==1) |
                        (PCV10_1==0 & PCV10_3==1) |
                        (PCV10_2==0 & PCV10_3==1)
        )

rotaissues1 <- analysis1 %>% 
        filter(rota1==0 & rota2==1)

length(unique(analysis1$ID))
length(unique(c(DTPissues1$ID, PCV10issues1$ID,rotaissues1$ID)))

#remove the IDs with non-sensical vaccination status and change the dates according to the 1's
analysis2 <- filter(analysis1, !ID %in% c(DTPissues1$ID, PCV10issues1$ID,rotaissues1$ID)) %>%
                        mutate(
                                DTPHBHiB1_date=ifelse(DTP1==1,DTPHBHiB1_date,NA),
                                PCV10_1_date = ifelse(PCV10_1==1,PCV10_1_date,NA),
                                rota_1_date = ifelse(rota1==1,rota_1_date,NA),
                                
                                DTPHBHiB2_date = ifelse(DTP2==1,DTPHBHiB2_date,NA),
                                PCV10_2_date =  ifelse(PCV10_2==1,PCV10_2_date,NA),
                                rota_2_date = ifelse(rota2==1,rota_2_date,NA),
                                
                                DTPHBHiB3_date = ifelse(DTP3==1,DTPHBHiB3_date,NA),
                                PCV10_3_date = ifelse(PCV10_3==1,PCV10_3_date,NA)
                                ) %>%
        #change characters to dates
        mutate(DTPHBHiB1_date=as.Date(DTPHBHiB1_date,format="%m/%d/%Y"),
                PCV10_1_date = as.Date(PCV10_1_date,format="%m/%d/%Y"),
                rota_1_date = as.Date(rota_1_date,format="%m/%d/%Y"),
                
                DTPHBHiB2_date = as.Date(DTPHBHiB2_date,format="%m/%d/%Y"),
                PCV10_2_date = as.Date(PCV10_2_date,format="%m/%d/%Y"),
                rota_2_date = as.Date(rota_2_date,format="%m/%d/%Y"),
                
                DTPHBHiB3_date = as.Date(DTPHBHiB3_date,format="%m/%d/%Y"),
                PCV10_3_date = as.Date(PCV10_3_date,format="%m/%d/%Y")
                
        ) 

#find those with dates that are out of order
DTPissues2 <- analysis2 %>% 
        filter(          DTPHBHiB1_date <= DDN |
                         DTPHBHiB2_date <= DTPHBHiB1_date |
                         DTPHBHiB3_date <= DTPHBHiB2_date )
PCV10issues2 <- analysis2 %>% 
        filter(         PCV10_1_date <= DDN |
                        PCV10_2_date <= PCV10_1_date |
                        PCV10_3_date <= PCV10_2_date )

rotaissues2 <- analysis2 %>% 
        filter(         rota_1_date <= DDN |
                        rota_2_date <= rota_1_date )


length(unique(analysis2$ID))
length(unique(c(DTPissues2$ID, PCV10issues2$ID,rotaissues2$ID)))

#remove those with dates out of order
analysis3 <- filter(analysis2, !ID %in% c(DTPissues2$ID, PCV10issues2$ID,rotaissues2$ID)) 
        

table(analysis3$DTP1,is.na(analysis3$DTPHBHiB1_date))
table(analysis3$DTP2,is.na(analysis3$DTPHBHiB2_date))
table(analysis3$DTP3,is.na(analysis3$DTPHBHiB3_date))


table(analysis3$PCV10_1,is.na(analysis3$PCV10_1_date))
table(analysis3$PCV10_2,is.na(analysis3$PCV10_2_date))
table(analysis3$PCV10_3,is.na(analysis3$PCV10_3_date))

table(analysis3$rota1,is.na(analysis3$rota_1_date))
table(analysis3$rota2,is.na(analysis3$rota_2_date))

#remove those with NA's for dates but have evidence of vaccination
analysis4 <- filter(analysis3,!(
                    (DTP1==1 & is.na(DTPHBHiB1_date))|
                    (DTP2==1 & is.na(DTPHBHiB2_date))|
                    (DTP3==1 & is.na(DTPHBHiB3_date))|
                            
                    (PCV10_1==1 & is.na(PCV10_1_date))|
                    (PCV10_2==1 & is.na(PCV10_2_date))|
                    (PCV10_3==1 & is.na(PCV10_3_date))|
                          
                    (rota1==1 & is.na(rota_1_date))|
                    (rota2==1 & is.na(rota_2_date)))
                    ) %>%
        mutate(rota3=0, #add third dose of rota to simplify code
               rota_3_date=as.Date(NA)
               )


tidyVariables<-function(data, antigen){
        #select the antigen we are looking at
        if (antigen == "DTP"){
                subdata <- data %>% 
                        select(ID,Region,District,HCC,DDN,
                               #DTP Vaccine
                               DTP1,DTP2,DTP3,DTPHBHiB1_date,DTPHBHiB2_date,DTPHBHiB3_date)
        }
        
        if (antigen == "PCV10"){
                subdata <- data %>% 
                        select(ID,Region,District,HCC,DDN,
                               #PCV Vaccine
                               PCV10_1,PCV10_2,PCV10_3, PCV10_1_date, PCV10_2_date, PCV10_3_date)
        }
        
        if (antigen == "Rota"){
                subdata <- data %>% 
                        select(ID,Region,District,HCC, DDN,
                               #Rota Vaccine
                               rota1,rota2, rota3, rota_1_date, rota_2_date,rota_3_date)
        }
        
        #change the variable names
        colnames(subdata)[6:11] <- c("V1","V2","V3","V1.date","V2.date","V3.date")
        
        
        #Fix geographic variables and birth
        subdata <- subdata %>% mutate(District=factor(District,
                                                      levels = c("MAHAJANGA I",
                                                                 "MAROVOAY" ,
                                                                 "TULEAR I" ,
                                                                 "TULEAR II",
                                                                 "MANAKARA",
                                                                 "VOHIPENO")),
                                      urban = ifelse(District %in% 
                                                              c("MAHAJANGA I",
                                                                "MANAKARA",
                                                                "TULEAR I"),
                                                      1,0)) %>%
                mutate(
                        birth.yearmonth=format(DDN, '%Y-%m'),
                        birth.year= year(DDN),
                        birth.month= month(DDN))
        
        #administratively censoring at data collection
        last.obs <- as.Date("2017-02-01") 
        subdata <- subdata %>%
                        mutate(
                                #ensure that all observations are before the administrative censoring date
                                V1 =ifelse(V1.date>last.obs,0,1),
                                V2 =ifelse(V2.date>last.obs,0,1),
                                V3 =ifelse(V3.date>last.obs,0,1),
                                V1.date=as.Date(ifelse(V1.date>last.obs,NA,as.character(V1.date))),
                                V2.date=as.Date(ifelse(V2.date>last.obs,NA,as.character(V2.date))),
                                V3.date=as.Date(ifelse(V3.date>last.obs,NA,as.character(V3.date)))
                        ) 
        
        
        # #create survival analysis variables
        SAdata <- subdata %>% mutate(
                        #Make cards who miss one dose to be listed as NA for future doses
                        #this doesnt match the equations
                        # V3 = ifelse(V2==0,NA,V3),
                        # V2 = ifelse(V1==0,NA,V2),

                        # Calculate censor date to determine if observation was censored or not
                        full_observe = ifelse(DDN + 365 < last.obs,TRUE,FALSE),
                        censor_date = as.Date(ifelse(full_observe,
                                                     as.character(DDN + 365),
                                                     as.character(last.obs))),
                        
                        # Determine if observation was censored
                        V1.SA.censor = ifelse(censor_date < V1.date | is.na(V1.date) ,
                                              0, V1),
                        V2.SA.censor = ifelse(censor_date < V2.date| is.na(V2.date), 0, V2),
                        V3.SA.censor = ifelse(censor_date < V3.date| is.na(V3.date), 0, V3),
                        
                        # V3.SA.censor = ifelse(V2.SA.censor==0, NA, V3.SA.censor),
                        # V2.SA.censor = ifelse(V1.SA.censor==0, NA, V2.SA.censor),
                        
                        #Define T1, T2 and T3
                        T1.SA.time =   ifelse(V1.SA.censor==1,
                                                         V1.date - DDN, censor_date - DDN),
                        T2.SA.time  =  ifelse(V2.SA.censor==1,
                                                         V2.date - DDN, censor_date - DDN),
                        T2.SA.time =   ifelse(V1.SA.censor==0,
                                              T1.SA.time,T2.SA.time),
                        T3.SA.time  =  ifelse(V3.SA.censor==1,
                                                         V3.date - DDN, censor_date - DDN),
                        T3.SA.time =   ifelse(V2.SA.censor==0,
                                              T2.SA.time,T3.SA.time),
                        
                        #Define B1, B2 and B3
                        B1.SA.time =  T1.SA.time,
                        B2.SA.time =  ifelse(V1.SA.censor==0, 0,
                                             T2.SA.time- T1.SA.time),
                        B3.SA.time =  ifelse(V2.SA.censor==0, 0,
                                             T3.SA.time- T2.SA.time)
                        
)

               

        # Define timeliness
        final <- SAdata %>%
                mutate(
                        V1.timeliness = ifelse(V1.SA.censor==0,
                                               ifelse(full_observe==TRUE,"Never","Censored"),
                                               as.character(cut(B1.SA.time,c(0,6,10,Inf),
                                                                right = FALSE,
                                                                labels=c("Early","OnTime","Late")))
                        ),
                        V2.timeliness = ifelse(V2.SA.censor==0,
                                               ifelse(full_observe==TRUE,"Never","Censored"),
                                               as.character(cut(B2.SA.time,c(0,4,8,Inf),
                                                                right = FALSE,
                                                                labels=c("Early","OnTime","Late")))
                        ),
                        V3.timeliness = ifelse(V3.SA.censor==0,
                                               ifelse(full_observe==TRUE,"Never","Censored"),
                                               as.character(cut(B3.SA.time,c(0,4,8,Inf),
                                                                right = FALSE,
                                                                labels=c("Early","OnTime","Late")))
                        )
                )%>%
                mutate(
                        V2.timeliness = ifelse(is.na(V2.timeliness),V1.timeliness,V2.timeliness),
                        V3.timeliness = ifelse(is.na(V3.timeliness),V2.timeliness,V3.timeliness)

                )

        return(SAdata)
        
} 



analysis.DTP <- tidyVariables(analysis4,"DTP")
tidyVariables(analysis4,"PCV10")
tidyVariables(analysis4,"Rota")



survfit(data = analysis.DTP,
        Surv(T1.SA.time,V1.SA.censor) ~ 1,
        type = "kaplan-meier",
        error = "greenwood",
        conf.type = "log-log") %>% plot()

survfit(data = analysis.DTP,
        Surv(T2.SA.time,V2.SA.censor) ~ 1,
        type = "kaplan-meier",
        error = "greenwood",
        conf.type = "log-log") %>% lines()

survfit(data = analysis.DTP,
        Surv(T3.SA.time,V3.SA.censor) ~ 1,
        type = "kaplan-meier",
        error = "greenwood",
        conf.type = "log-log") %>% lines()



                    

