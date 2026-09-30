#Data analysis for "Global relationships between non-crop habitat, crop pollinator diversity and crop pollination"
#Last modified 2026-09-22

#Load packages
library(dplyr)
library(nlme)
library(ggplot2)
library(mapdata) 
library(metafor) 
library(patchwork)
library(caret) 
library(vegan) 
library(mgcv) 
library(tidyr) 
library(ggeffects)

#Graph theme
mytheme <- theme_bw()+ theme(panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(),
                             legend.background = element_rect(fill = "transparent"),
                             panel.background=element_rect(fill = "white"))

# 1. Read data and sorting
    Dat <- read.csv("Data.csv") 
    Dat$NonCrop <- with(Dat,(Herbaceous+Trees+Shrubland+Grassland+Mangroves)/(Total))
    SelectYield  <-  read.csv ("Yield_matrix.csv",na = "NA") %>%
                      filter((Type1+Type2)!=0 )
    
# Calculate the Service matrix
    Dat$Measure1 <- with(Dat,log(yield/yield_t)) #Services according to yield 1 - log ratio
    Dat$Measure2 <- with(Dat,log(yield2/yield2_t)) #Services according to yield 2 - log ratio
    
    Dat_long <- pivot_longer(Dat,
      cols = c(Measure1, Measure2),  # columns to pivot
      names_to = "MeasureType",       # new column name for type names
      values_to = "Service"  )     # new column name for values

    Sub <- Dat_long[complete.cases(Dat_long$Service*0),,drop=FALSE] #remove NA and inf values
    Sub$ServiceType <- ifelse(Sub$MeasureType == "Measure1",
                        SelectYield$Type1[match(Sub$study_id, SelectYield$study_id)],
                        SelectYield$Type2[match(Sub$study_id, SelectYield$study_id)])
    #Data subset
    Sub1 <- droplevels(subset(Sub,ServiceType == 1)) # pollination services1
    Sub2 <- droplevels(subset(Sub,ServiceType == 2)) # pollination services2

# 2. Plot the location of each study
    Study_Location_1 <- aggregate(cbind(latitude,longitude)~study_id,data=Sub1,mean)
    Study_Location_2 <- aggregate(cbind(latitude,longitude)~study_id,data=Sub2,mean)

    p <- ggplot() + coord_fixed() +
        xlab("") + ylab("")
    base_world <- map_data("world")
   p_site <-  p + geom_polygon(data=base_world, 
                    aes(x=long, y=lat, group=group), 
                    colour="light grey", fill="grey90")+
        geom_point(data=Study_Location_1, 
                   aes(x=longitude, y=latitude),col="black",alpha=0.8)+
        geom_point(data=Study_Location_2, 
                   aes(x=longitude, y=latitude),col="black",alpha=0.8)+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_rect(fill = 'white', colour = 'white'), 
              axis.line = element_line(colour = "white"), 
              axis.ticks = element_blank(), 
              axis.text.x = element_blank(),
              axis.text.y = element_blank())
   p_site


# 3. Overall pollination services - Figure 1
  # Modeling
    mod_reproductive_ratio <- lme(Service~1,random=~1|as.factor(study_climate),data=Sub1)
    LR1 <- coef(summary(mod_reproductive_ratio))
    LR1.study <- data.frame(value = coef(mod_reproductive_ratio),type ='reproductive') #Log ratio 1 at study level
    mod_weight_ratio <- lme(Service~1,random=~1|as.factor(study_climate),data=Sub2)
    summary(mod_weight_ratio)
    LR2 <- coef(summary(mod_weight_ratio))
    LR2.study <- data.frame(value = coef(mod_weight_ratio),type ='weight')    #Log ratio 2 at study level
    #Residual check
    plot(mod_reproductive_ratio)
    plot(mod_weight_ratio)
    
    #New data frame
    Ldf <- rbind(LR1.study,LR2.study)
    
    #Plot the results
    Ldf$type <- as.factor(Ldf$type)
    levels(Ldf$type) <- c("Reproductive-based","Weight-based")
    
    p_effect <- ggplot(data=NULL,aes(x,y))+
      geom_jitter(data=Ldf,aes(x=type,y=X.Intercept.),width = 0.1,alpha=0.1)+
      geom_point(aes(x=levels(Ldf$type)[1],y=LR1[1]),col="black",cex=3)+
      geom_errorbar(aes(x=levels(Ldf$type)[1],
                        y=LR1[1],
                        ymax=LR1[1]+1.96*LR1[1,2],
                        ymin=LR1[1]-1.96*LR1[1,2]),width=0.05,col="black")+
      geom_point(aes(x=levels(Ldf$type)[2],y=LR2[1]),col="black",cex=3)+
      geom_errorbar(aes(x=levels(Ldf$type)[2],
                        y=LR2[1],
                        ymax=LR2[1]+1.96*LR1[1,2],
                        ymin=LR2[1]-1.96*LR1[1,2]),width=0.05,col="black")+
      geom_hline(yintercept = 0,linetype = "dashed")+
      ylab("Effect size")+
      xlab("Type of pollination services")+
      mytheme
    p_effect

    
# 4 Pollination services for different crop - Figure 2
    Result_crop_S1 <- data.frame()
    for (i in unique(Sub1$crop)){
      Sub1.i <- subset(Sub1, crop == i)
      Mcrop_S1 <- lme(Service~1,random=~1|as.factor(study_climate),data=Sub1.i) 
      Coef_S1 <- coefficients(summary(Mcrop_S1))
      Out_S1 <- data.frame( 
        Type = "Reproductive",
        Crop=i,
        n = length(unique(Sub1.i$study_climate)),
        beta = Coef_S1[1,1],
        SE = Coef_S1[1,2],
        p = Coef_S1[1,5])
      Result_crop_S1 <- rbind(Result_crop_S1, Out_S1)
      print(i)
    }
    
    Result_crop_S2 <- data.frame()
    for (i in unique(Sub2$crop)){
      Sub2.i <- subset(Sub2, crop == i)
      Mcrop_S2 <- lme(Service~1,random=~1|as.factor(study_climate),data=Sub2.i) 
      Coef_S2 <- coefficients(summary(Mcrop_S2))
      Out_S2 <- data.frame( 
        Type = "Weight",
        Crop=i,
        n = length(unique(Sub2.i$study_climate)),
        beta = Coef_S2[1,1],
        SE = Coef_S2[1,2],
        p = Coef_S2[1,5]
      )
      Result_crop_S2 <- rbind(Result_crop_S2, Out_S2)
      print(i)
    }
      
    Overall_S1 <- data.frame(Type = "Reproductive",
                             Crop = "Overall", 
                           n = length(unique(Sub1$study_climate)),
                           beta = LR1[1,1],
                           SE = LR1[1,2],
                           p = LR1[1,5])
    
    Overall_S2 <- data.frame(Type = "Weight",
                           Crop="Overall", 
                           n = length(table(Sub2$study_climate)>=5),
                           beta = LR2[1,1],
                           SE = LR2[1,2],
                           p = LR2[1,5])
    
    R_Crop <- rbind(Result_crop_S1,Result_crop_S2,Overall_S1,Overall_S2)
    R_Crop$CI <- 1.96*R_Crop$SE
    R_Crop$Sig <- with(R_Crop,ifelse(p<=0.05&p>0.01,"*",
                                            ifelse(p<=0.01&p>0.001,"**",
                                                   ifelse(p<=0.001,"***","NS"))))
    R_Crop$Crop_n <- paste0(R_Crop$Sig,"(",R_Crop$n,")")
    
    
    #Data sorting
    R_Crop$Crop <- as.factor(R_Crop$Crop)
    R_Crop$Crop <- factor(R_Crop$Crop, levels = sort(levels(R_Crop$Crop),decreasing = F))
    R_Crop$Crop <- relevel(R_Crop$Crop,ref = "Overall")
    R_Crop$Type <- as.factor(R_Crop$Type)
    R_Crop$Benefit <- 1/exp(R_Crop$beta)-1
    R_Crop$Increase <- exp(R_Crop$beta)-1
    R_Crop$Rank <- rank(R_Crop$beta)
    
    #Plot
    p_service_crop <- ggplot(data= subset(R_Crop,Crop != "Overall"),
                             aes(x=beta,y= Crop))+
      geom_vline(xintercept=0,lty=3,col="black")+
       geom_errorbarh(data= subset(R_Crop,Crop == "Overall"),
                        aes(xmin=beta-CI, xmax=beta+CI),col="red",
                      height=0.1, position=position_dodge(width=0.5))+
       geom_point(data= subset(R_Crop,Crop == "Overall"), col = "red",
                  position=position_dodge(width=0.5))+
      geom_text(data= subset(R_Crop,Crop == "Overall"), col = "red",
        aes(label = Crop_n,group=Type),size = 2,vjust = -0.5,position=position_dodge(width=0.5)) +
      geom_errorbarh(aes(xmin=beta-CI, xmax=beta+CI),
                     height=0.1, position=position_dodge(width=0.5),col = "black")+
      geom_point(position=position_dodge(width=0.5),col = "black")+
      geom_text(aes(label = Crop_n,group=Type),size = 2,vjust = -0.5,position=position_dodge(width=0.5)) +
      ylab("")+
      xlab("Effect size")+
      labs(shape="")+
      guides(colour = FALSE)+
      facet_grid(~Type,drop = TRUE)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(face = "italic"),
            legend.position=c(0.8,0.1),
            legend.background = element_rect(fill = "transparent"),
            panel.background=element_rect(fill = "white"))
    p.service.crop
    
    
    
#  5. Pollinator diversity, data sorting and merging
    Dat_pol <- read.csv("CropPol_sampling_data.csv")
    Dt <- with(Dat_pol,aggregate(list(abundance=abundance),
                                 list(study_id=study_id,site_id=site_id,pollinator=pollinator),
                                 sum))
    #Calculate the diversity
    Study_site_diversity <- aggregate(abundance~site_id+study_id,
                                      FUN=diversity,
                                      data=Dat_pol) 
    names(Study_site_diversity)[3] <- "Shannon"
    Mgdat <- merge(Dat,Study_site_diversity,by=c("study_id","site_id"), all=T)
    Mgdat$Shannon.d <- exp(Mgdat$Shannon)
    
    
# 6. Check the relationship between Shannon diversity and pollination services
    #Overall
    Diversity_service <- merge(Sub,Study_site_diversity,by=c("study_id","site_id"), all=T)
    Diversity_service$Shannon.d <- exp(Diversity_service$Shannon)
    
    Diversity_service_dt1 <- Diversity_service %>%
      filter(is.finite(Shannon.d)&(ServiceType==1))
    Diversity_service_dt2 <- Diversity_service %>%
      filter(is.finite(Shannon.d)&(ServiceType==2))
    
    Diversity_service_m1 <- lme(Service~Shannon.d,
              random =~1|factor(study_climate), #use study_climate as random factor
              data=Diversity_service_dt1,
              method="REML",
              na.action = na.omit)
    Diversity_service_m2 <- lme(Service~Shannon.d,
               random =~1|factor(study_climate),
               data=Diversity_service_dt2,
               method="REML",
               na.action = na.omit)
    
    summary(Diversity_service_m1) # 
    summary(Diversity_service_m2) # Both model shows positive but not significant

    
  # Plot - Appendix 4
    #Service 1
    #Get the 95%CI for the model
    Diversity_service_m1coef <- coef(summary(Diversity_service_m1))
    x <- Diversity_service_dt1$Shannon.d
    xv<-seq(min(x),max(x),(max(x)- min(x))/100)
    yv<- Diversity_service_m1coef[1,1]+Diversity_service_m1coef[2,1]*xv 
    xm<-mean(x)
    n<-length(x)
    ssx<- sum(x^2)- sum(x)^2/n
    s.t<- qt(0.975,(n-2)) 
    se<-sqrt(summary(Diversity_service_m1)[[6]]^2*(1/n+(xv-xm)^2/ssx)) 
    ci<-s.t*se
    uyv<-yv+ci
    lyv<-yv-ci
    
    #Create a data.frame
    Pred <- data.frame(xv,yv,uyv,lyv) 
    
    #Plot the results 
    Diversity_service_p1 <-     ggplot(Diversity_service_dt1,aes(x=Shannon.d,y=Service))+
      xlab("Shannon")+
      ylab("Service")+
      geom_hline(yintercept=0,lty=2,col="red")+
      geom_point(size=2, col = "grey50", alpha = 0.4)+
      geom_line(data = Pred, aes(x= xv, y = yv), col  = "blue", size = 1)+ 
      geom_ribbon(data = Pred,aes(x= xv, y = yv, ymax=uyv,ymin=lyv),alpha=0.3,color="grey")  +
      xlab("Pollinator diversity")+
      ylab("Pollination services")+
      ggtitle("(a) Reproductive-based")+
      theme_bw()
    
    Diversity_service_p1
    
    #Service 2    
    #Get the 95%CI for the model
    Diversity_service_m2coef <- coef(summary(Diversity_service_m2))
    
    x <- Diversity_service_dt2$Shannon.d
    xv<-seq(min(x),max(x),(max(x)- min(x))/100)
    yv<- Diversity_service_m2coef[1,1]+Diversity_service_m2coef[2,1]*xv 
    
    xm<-mean(x)
    n<-length(x)
    ssx<- sum(x^2)- sum(x)^2/n
    s.t<- qt(0.975,(n-2)) 
    
    se<-sqrt(summary(Diversity_service_m2)[[6]]^2*(1/n+(xv-xm)^2/ssx)) 
    ci<-s.t*se
    uyv<-yv+ci
    lyv<-yv-ci
    
    #Build up a data.frame
    Pred <- data.frame(xv,yv,uyv,lyv) #Note: this is very similar to the build-in geom_smooth() function
    
    #Plot the results by ggplot2
    Diversity_service_p2 <-     ggplot(Diversity_service_dt2,aes(x=Shannon.d,y=Service))+
      xlab("Shannon")+
      ylab("Service")+
      geom_hline(yintercept=0,lty=2,col="red")+
      geom_point(size=2, col = "grey50", alpha = 0.4)+
      geom_line(data = Pred, aes(x= xv, y = yv), col  = "blue", size = 1)+ 
      geom_ribbon(data = Pred,aes(x= xv, y = yv, ymax=uyv,ymin=lyv),alpha=0.3,color="grey")  +
      xlab("Pollinator diversity")+
      ylab("")+
      ggtitle("(b) Weight-based ")+
      theme_bw()    
    
    
# 7. Check the relationship between landscape and pollinator diversity (Shannon index)
    #Overall
    subdat_diversity <- Mgdat[complete.cases(Mgdat$Shannon*0),,drop=FALSE] 
    d_noncrop_m <- lme(Shannon~scale(NonCrop), 
              random =~1|factor(study_climate),
              data=subdat_diversity,
              method="REML",
              na.action = na.omit)
    summary(d_noncrop_m) 
    
    # Landscape effect in different climate zones
    Climate_all_shannon <- unique(subdat_diversity$Climate)
    Result_climate_shannon <- data.frame()
    for (i in 1:length(Climate_all_shannon))
    {
      climate_i <- Climate_all_shannon[i] 
      Sub_climate <- droplevels(subset(subdat_diversity,Climate==climate_i)) 
      if (nrow(Sub_climate)<5){next}
      
      mod <- lme(Shannon~scale(NonCrop),
                 random =~1|factor(study_id), 
                 data=Sub_climate,
                 method="REML",
                 na.action = na.omit)
      Coef <- coefficients(summary(mod))
      
      Out <- data.frame( 
        Climate=climate_i,
        n = length(unique(Sub_climate$study_id)),
        beta = Coef[2,1],
        SE = Coef[2,2],
        p = Coef[2,5]
      )
      Result_climate_shannon <- rbind(Result_climate_shannon,Out)
    }
    
    Overall_shannon <- data.frame( Climate="Overall", 
                                   n = length(table(subdat_diversity$study_climate)>=5),
                                   beta = coef(summary(d_noncrop_m))[2,1],
                                   SE = coef(summary(d_noncrop_m))[2,2],
                                   p = coef(summary(d_noncrop_m))[2,5])
    
    Rshannon_climate<- rbind(Result_climate_shannon,Overall_shannon)
    Rshannon_climate$CI <- 1.96*Rshannon_climate$SE
    Rshannon_climate$Sig <- with(Rshannon_climate,
                                 ifelse(p<=0.05&p>0.01,"*",
                                        ifelse(p<=0.01&p>0.001,"**",
                                               ifelse(p<=0.001,"***","NS"))))
    Rshannon_climate$Climate_n <- paste0(Rshannon_climate$Sig,"(",Rshannon_climate$n,")")
    
    
    Rshannon_climate$Climate <- factor(Rshannon_climate$Climate,levels = c("Overall","D","C","B","A"))
    levels( Rshannon_climate$Climate) <- c("Overall","Cold", "Temperate", "Arid", "Tropical")
    
    #Plot the results
    p_shannon_climate <- 
      ggplot(data=Rshannon_climate,aes(x=beta,y=Climate,colour=Climate))+ 
      geom_vline(xintercept=0,lty=3,col="blue")+
      geom_errorbarh(aes(xmin=beta-CI, xmax=beta+CI),height=0.1,position=position_dodge(width=0.5))+
      geom_point(position=position_dodge(width=0.5))+
      scale_color_manual(values = c( "red","black","black","black","black")) +
      geom_text(aes(label = Climate_n),size = 2,vjust = -0.5,position=position_dodge(width=0.5)) +
      ylab("Climate")+
      xlab("Coefficient")+
      labs(shape="")+
      ggtitle("(a) Pollinator diversity")+
      guides(colour = FALSE)+
      theme_bw()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background=element_rect(fill = "white"))+
      theme(axis.text.x = element_text(colour="black",angle = 0, hjust = 0.5,vjust = 0.5))
    p_shannon_climate
    
    #For individual study
    study_all <- unique(subdat_diversity$study_climate)
    Result_div_study <- data.frame() 
    for (i in 1:length(study_all))
    {
      study_i <- study_all[i] 
      Sub_study <- subset(subdat_diversity,study_climate==study_i) 
      if (nrow(Sub_study)<5){next}
      if(length(nearZeroVar(Sub_study$NonCrop)!=0)){next} 
      mod <- glm(Shannon~scale(NonCrop),data=Sub_study) 
      Coef <- coefficients(summary(mod))
      Out <- data.frame( 
        study = study_i,
        Climate = unique(Sub_study$Climate),
        n = nrow(Sub_study),
        beta = Coef[2,1],
        SE = Coef[2,2],
        MedCrop = median(Sub_study$NonCrop,na.rm=T)
      )
      Result_div_study <- rbind(Result_div_study,Out)
      print(i)
    }

    #Plot at different climate zone - overall model 
    Beta_all <- coef(summary(d_noncrop_m))[2,1]
    SE_all <- coef(summary(d_noncrop_m))[2,2]
    
    #Plot along different climate zones - individual model
    Result_div_study$study_name <- as.factor(paste0("Study_shannon",1:length(Result_div_study$study)))
    Result_div_study$study_name <- factor (Result_div_study$study_name,
                                  levels=Result_div_study$study_name[order(Result_div_study$beta)])
    Result_div_study$Climate <- as.factor(substr(Result_div_study$Climate,1,1)) 
    
    Result_div_study <- Result_div_study %>%
      filter(is.finite(beta),is.finite(SE),SE > 0,is.finite(MedCrop)) %>%
      mutate(variance = SE^2, precision = 1 / variance,weight_iv = precision / mean(precision) )

    # GAM model    
    # Unweighted
    GAM_crop1 <- gam(beta ~ s(MedCrop, k = 10),data = Result_div_study,  method = "REML")
    # Weighted
    GAM_crop2 <- gam( beta ~ s(MedCrop, k = 10), data = Result_div_study, weights = weight_iv, method = "REML")
    summary(GAM_crop1)
    summary(GAM_crop2)
    gam.check(GAM_crop2)
    
    # Plot unweighted
    p_diversity_unweighted <- ggplot(subset(Result_div_study),aes(x=MedCrop,y=beta))+
      xlab("Non-crop habitat coverage")+
      ylab("Within-study slope")+
      ggtitle("(a) Pollinator diversity")+
      geom_hline(yintercept=0,lty=2,col="red")+
      geom_point(size=2,col = "grey50")+
      geom_smooth(method = "gam")+
      theme_bw()
    p_diversity_unweighted  
    
    # Plot weighted
    new_crop <- data.frame(MedCrop = seq(min(Result_div_study$MedCrop, na.rm = TRUE),max(Result_div_study$MedCrop, na.rm = TRUE), length.out = 200))
    pred_crop <- predict(GAM_crop2, newdata = new_crop,type = "response", se.fit = TRUE)
    new_crop$fit <- as.numeric(pred_crop$fit)
    new_crop$lower <- new_crop$fit - 1.96 * pred_crop$se.fit
    new_crop$upper <- new_crop$fit + 1.96 * pred_crop$se.fit
    p_diversity_weighted <- ggplot() +
      geom_hline(yintercept = 0, linetype = 2, colour = "red") +
      geom_ribbon(data = new_crop,   aes(x = MedCrop, ymin = lower, ymax = upper), fill = "grey50",  alpha = 0.35) +
      geom_line( data = new_crop,  aes(x = MedCrop, y = fit), colour = "blue" ) +
      geom_errorbar(  data = Result_div_study,  aes(x = MedCrop,ymin = beta -  SE, ymax = beta +  SE),
                      width = 0, colour = "grey70",  alpha = 0.5 ) +
      geom_point( data = Result_div_study, aes(x = MedCrop, y = beta),size = 2, colour = "grey40" ) +
      ggtitle("(a) Pollinator diversity")+
      labs( x = "Non-crop habitat coverage",
            y = "Within-study slope",
            shape = NULL) +
      mytheme
    
    p_diversity_weighted

# 8. Check the relationship between landscape and pollination services
#  Overall
    M1 <- lme(Service~scale(NonCrop),
              random =~1|factor(study_climate),
              data=Sub1,
              method="REML",
              na.action = na.omit)
    M2 <- lme(Service~NonCrop,
              random =~1|factor(study_climate),
              data=Sub2,
              method="REML",
              na.action = na.omit)
    summary(M1) 
    summary(M2) 
    
#Different climate zone  -  service 1
    Climate_all_s1 <- unique(Sub1$Climate)
    Result1_climate <- data.frame()
    for (i in 1:length(Climate_all_s1))
    {
      climate_i <- Climate_all_s1[i] 
      Sub_climate <- droplevels(subset(Sub1,Climate==climate_i)) 
      if (nrow(Sub_climate)<5){next}
      mod <- lme(Service~scale(NonCrop),
                random =~1|factor(study_id), 
                data=Sub_climate,
                method="REML",
                na.action = na.omit)
      Coef <- coefficients(summary(mod))
      Out <- data.frame( 
        Type = "Reproductive",
        Climate=climate_i,
        n = length(unique(Sub_climate$study_id)),
        beta = Coef[2,1],
        SE = Coef[2,2],
        p = Coef[2,5]
      )
      Result1_climate <- rbind(Result1_climate,Out)
      print(i)}

    #Climate - service 2
    Climate_all_s2 <- unique(Sub2$Climate)
    Result2_climate <- data.frame()
    for (i in 1:length(Climate_all_s2))
    {
      climate_i <- Climate_all_s2[i] 
      Sub_climate <- droplevels(subset(Sub2,Climate==climate_i)) # ***
      if (nrow(Sub_climate)<5){next}
      
      mod <- lme(Service~scale(NonCrop),
                 random =~1|factor(study_id), 
                 data=Sub_climate,
                 method="REML",
                 na.action = na.omit)
      coefficients(mod)
      Coef <- coefficients(summary(mod))
      
      Out <- data.frame( 
        Type = "Weight",
        Climate=climate_i,
        n = length(unique(Sub_climate$study_id)),
        beta = Coef[2,1],
        SE = Coef[2,2],
        p = Coef[2,5]
      )
      Result2_climate <- rbind(Result2_climate,Out)
      print(i)
    }
    
   
   Overall1 <- data.frame(Type = "Reproductive",
                          Climate="Overall", 
                          n = length(table(Sub1$study_climate)>=5),
                         beta = coef(summary(M1))[2,1],
                         SE = coef(summary(M1))[2,2],
                         p = coef(summary(M1))[2,5])
   
   Overall2 <- data.frame(Type = "Weight",
                          Climate="Overall", 
                          n = length(table(Sub2$study_climate)>=5),
                          beta = coef(summary(M2))[2,1],
                          SE = coef(summary(M2))[2,2],
                          p = coef(summary(M2))[2,5])

   R12_climate <- rbind(Result1_climate,Result2_climate,Overall1,Overall2)
   R12_climate$CI <- 1.96*R12_climate$SE
   R12_climate$Sig <- with(R12_climate, ifelse(p<=0.05&p>0.01,"*",
                                                     ifelse(p<=0.01&p>0.001,"**",
                                                            ifelse(p<=0.001,"***","NS"))))
   R12_climate$Climate_n <- paste0(R12_climate$Sig,"(",R12_climate$n,")")
   
   R12_climate$Climate <- factor(R12_climate$Climate,levels = c("Overall","D","C","B","A"))
   levels(R12_climate$Climate) <- c("Overall","Cold", "Temperate", "Arid", "Tropical")
   
   #Plot the results
   R12_climate$Type <- as.factor(R12_climate$Type)
   
   p_service_climate <- ggplot(data=R12_climate,aes(x=beta,y=Climate,shape=Type,colour=Climate))+ #Two studies were removed during the ploting but should be in appendix
     geom_vline(xintercept=0,lty=3,col="blue")+
     geom_errorbarh(aes(xmin=beta-CI, xmax=beta+CI),height=0.1,position=position_dodge(width=0.5))+
     geom_point(position=position_dodge(width=0.5))+
     scale_color_manual(values = c( "red","black","black","black","black")) +
     geom_text(aes(label = Climate_n,group=Type),size = 2,vjust = -0.5,position=position_dodge(width=0.5)) +
     ylab("")+
     xlab("Coefficient")+
     ggtitle("(b) Pollination services")+
     theme(panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background=element_rect(fill = "white"))+
     labs(shape="")+
     guides(colour = FALSE)+
       theme_bw()+
     theme(panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background=element_rect(fill = "white"),
           legend.position=c(0.8,0.1),
           legend.background = element_rect(fill = "transparent"))+
     theme(axis.text.x = element_text(colour="black",angle = 0, hjust = 0.5,vjust = 0.5))
   
   p_service_climate
   

# For individual study, i.e. beta_i
    #Service1
    study_all <- unique(Sub1$study_climate)
    

    Result_study_service1 <- data.frame() 
    for (i in 1:length(study_all))
    {
        study_i <- study_all[i] 
        Sub_study <- subset(Sub1,study_climate==study_i) # ***
        if (nrow(Sub_study)<5){next}
        if(length(nearZeroVar(Sub_study$Service)!=0)){next} #Remove studies without variance
        if(length(nearZeroVar(Sub_study$NonCrop)!=0)){next} #Remove studies without variance
        mod <- glm(Service~scale(NonCrop),data=Sub_study) # Standardized within each study
        Coef <- coefficients(summary(mod))
        Out <- data.frame( 
            study = study_i,
            Climate=unique(Sub_study$Climate),
            n = nrow(Sub_study),
            beta = Coef[2,1],
            SE = Coef[2,2],
            MedCrop = median(Sub_study$NonCrop,na.rm=T))
        Result_study_service1 <- rbind(Result_study_service1,Out)
        print(i)
    }

    
    #Plot along different climate zones
    Result_study_service1$study_name <- as.factor(paste0("Study_reproductive",1:length(Result_study_service1$study)))
   Result_study_service1$study_name <- factor (Result_study_service1$study_name,
                                levels=Result_study_service1$study_name[order(Result_study_service1$beta)])
    Result_study_service1$Climate <- as.factor(substr(Result_study_service1$Climate,1,1)) 
    
    
    #Service2
    study_all2 <- unique(Sub2$study_climate)
    
    Result_study_service2 <- data.frame() 
    for (i in 1:length(study_all2))
    {
        study_i <- study_all2[i] 
        Sub_study <- subset(Sub2,study_climate==study_i)
        if (nrow(Sub_study)<5){next} 
        if(length(nearZeroVar(Sub_study$Service)!=0)){next} 
        if(length(nearZeroVar(Sub_study$NonCrop)!=0)){next} 
        mod <- glm(Service~scale(NonCrop),data=Sub_study)
        Coef <- coefficients(summary(mod))
        Out <- data.frame( 
                    study = study_i,
                    Climate=unique(Sub_study$Climate),
                    n = nrow(Sub_study),
                    beta = Coef[2,1],
                    SE = Coef[2,2],
                    MedCrop = median(Sub_study$NonCrop,na.rm=T)
                    )
        Result_study_service2 <- rbind(Result_study_service2,Out)
        print(i)
    }

    Result_study_service2$study_name <- as.factor(paste0("Study_weight",1:length(Result_study_service2$study)))
    Result_study_service2$study_name <- factor (Result_study_service2$study_name,
                                  levels=Result_study_service2$study_name[order(Result_study_service2$beta)])
    Result_study_service2$Climate <- as.factor(substr(Result_study_service2$Climate,1,1)) 

    #Combine Service1 and Service2 together
    Result_study_service1$Type <- "Reproductive"
    Result_study_service2$Type <- "Weight"
    Result <- rbind(Result_study_service1, Result_study_service2) %>%
      filter(is.finite(beta),is.finite(SE),SE > 0,is.finite(MedCrop)) %>%
      mutate(variance = SE^2, precision = 1 / variance,weight_iv = precision / mean(precision)) 
        
  # GAM 
  # Unweighted 
    GAM1 <- gam(beta ~ s(MedCrop, k = 10),data = Result,  method = "REML")
  # Weighted 
    GAM2 <- gam( beta ~ s(MedCrop, k = 10), data = Result, weights = weight_iv, method = "REML")
    summary(GAM1)
    summary(GAM2)
    gam.check(GAM2)
    
    # Plot unweighted model
    p_service_unweighted <- ggplot(subset(Result),aes(x=MedCrop,y=beta))+
      geom_hline(yintercept=0,lty=2,col="red")+
      geom_point(size=2,aes(shape=Type),col = "grey50")+
      geom_smooth(method = "gam")+
      labs(shape="")+
      ggtitle("(b) Pollination services")+
      xlab("Non-crop habitat coverage")+
      ylab("")+
      theme_bw()+ 
      theme(legend.position=c(0.16,0.12),
            legend.background = element_rect(fill = "transparent"),
            panel.background=element_rect(fill = "white"))
    p_service_unweighted
    
    #Plot weighted model 
    new_service <- data.frame(MedCrop = seq(min(Result$MedCrop, na.rm = TRUE),
                                            max(Result$MedCrop, na.rm = TRUE), length.out = 200))
    pred_service <- predict(GAM2, newdata = new_service,type = "response", se.fit = TRUE)
    
    new_service$fit <- as.numeric(pred_service$fit)
    new_service$lower <- new_service$fit - 1.96 * pred_service$se.fit
    new_service$upper <- new_service$fit + 1.96 * pred_service$se.fit
    
    p_service_weighted <- ggplot() +
      geom_hline(yintercept = 0, linetype = 2, colour = "red") +
      geom_ribbon(data = new_service,   aes(x = MedCrop, ymin = lower, ymax = upper), fill = "grey50",  alpha = 0.35) +
      geom_line( data = new_service,  aes(x = MedCrop, y = fit), colour = "blue" ) +
      geom_errorbar(  data = Result,  aes(x = MedCrop,ymin = beta -  SE, ymax = beta +  SE),
                      width = 0, colour = "grey70",  alpha = 0.5 ) +
      geom_point( data = Result, aes(x = MedCrop, y = beta, shape = Type),size = 2, colour = "grey40" ) +
      ggtitle("(b) Pollination services")+
            labs( x = "Non-crop habitat coverage",
            y = "",
            shape = NULL) +
      mytheme +
      theme(legend.position=c(0.16,0.12))
    
    p_service_weighted
    

# 9. Plot the theoretical effect directly from the model - Figure 3
    #Pollinator diversity
    d_noncrop_p <-   plot(ggeffect(d_noncrop_m, "NonCrop")) +
      labs( x = "Non-corp habitat coverage", y = "Pollinator divesity", title = "(a)")  +
      mytheme
    
    # Pollination services
    p1 <- ggeffect(M1, terms = "NonCrop") %>% mutate(model = "Reproductive")
    p2 <- ggeffect(M2, terms = "NonCrop") %>% mutate(model = "Weight")
    pred <- bind_rows(p1, p2)
    
    # With CI ribbons + lines
    s_noncrop_p <-   ggplot(pred, aes(x = x, y = predicted, color = model, fill = model)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 1,lty = "dashed") +
      labs(x = "Non-corp habitat coverage", 
           y = "Pollination Services",
           title = "(b)",
           colour = "",
           fill = "") +
      mytheme +
      theme(legend.position=c(0.8,0.92))

# 10. Plots in the manuscript
    pdf("Figure 1.pdf",width=3.5,height=3)
    p_effect
    dev.off()
    
    pdf("Figure 2.pdf",width=5.5,height=4.0) 
    p_service_crop 
    dev.off()
    
    pdf("Figure 3.pdf",width=9,height=4.5)
    print(d_noncrop_p|s_noncrop_p) 
    dev.off()
    
    pdf("Figure 4.pdf",width=7.5,height=4.5)
    print(p_shannon_climate|p_service_climate)
    dev.off()
    
    pdf("Figure 5.pdf",width=9,height=4.5)
    print(p_diversity_unweighted|p_service_unweighted) 
    dev.off()
    
    pdf("Appendix 3.pdf",heigh=6,width=8)
    print(p_site) 
    dev.off()
    
    pdf("Appendix 4.pdf",width=9,height=4.5)
    print(Diversity_service_p1|Diversity_service_p2)
    dev.off()
    
    pdf("Appendix 5.pdf",width=9,height=4.5)
    print(p_diversity_weighted|p_service_weighted) 
    dev.off()
    