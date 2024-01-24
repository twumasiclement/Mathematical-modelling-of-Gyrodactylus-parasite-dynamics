#Function for extract experimental descriptors
Experiment_descriptors<-function(empirical_data){
#Fish-parasite combinations/groups
parasite_fish=c("Gt3-OS","Gt3-LA","Gt3-UA","Gt-OS","Gt-LA","Gt-UA","Gb-OS","Gb-LA","Gb-UA") 
levels(empirical_data$Sex_fish)<-c("F","M")
empirical_data$LowerRegion<-empirical_data$LB+empirical_data$Pelvic+empirical_data$Anal+empirical_data$Dorsal
empirical_data$UpperRegion<-empirical_data$UB +Combined_data$Pectoral

# Data across the four recategorized body regions
Data_fourRegions<-empirical_data[,c(1,15,16,8,9,10,12,13,11,14)]
#head(Data_fourRegions,n=4) 

#Data across parasite strains
Gt3_data<-Data_fourRegions[Data_fourRegions$Parasite_strain=="Gt3",]
Gt_data<-Data_fourRegions[Data_fourRegions$Parasite_strain=="Gt",]
Gb_data<-Data_fourRegions[Data_fourRegions$Parasite_strain=="Gb",]


Parasite_fish_data=NULL;fishID=NULL;numF=NULL;pop_obs=NULL;alive_obs=NULL;fishSize=NULL;fishSex=NULL
Size=NULL;Sex=NULL; Parasite_strain=NULL;Strain=NULL;Fish_type=NULL;Fish_stock=NULL
#Data of each parasite strain across fish stocks

Parasite_fish_data[[parasite_fish[1]]]<-split(Gt3_data,Gt3_data$Fish_strain)$"OS"
Parasite_fish_data[[parasite_fish[2]]]<-split(Gt3_data,Gt3_data$Fish_strain)$"LA"
Parasite_fish_data[[parasite_fish[3]]]<-split(Gt3_data,Gt3_data$Fish_strain)$"UA"
Parasite_fish_data[[parasite_fish[4]]]<-split(Gt_data,Gt_data$Fish_strain)$"OS"
Parasite_fish_data[[parasite_fish[5]]]<-split(Gt_data,Gt_data$Fish_strain)$"LA"
Parasite_fish_data[[parasite_fish[6]]]<-split(Gt_data,Gt_data$Fish_strain)$"UA"
Parasite_fish_data[[parasite_fish[7]]]<-split(Gb_data,Gb_data$Fish_strain)$"OS"
Parasite_fish_data[[parasite_fish[8]]]<-split(Gb_data,Gb_data$Fish_strain)$"LA"
Parasite_fish_data[[parasite_fish[9]]]<-split(Gb_data,Gb_data$Fish_strain)$"UA"

for (pf in 1:length(parasite_fish)){
    #Assigning unique ID for  data
    fishID[[pf]]<- unique(Parasite_fish_data[[parasite_fish[pf]]]$Fish_ID) 
    #Total number of fish used for data
    numF[[pf]] <- length(fishID[[pf]]) 
    #Observed data or matrix across 4 regions
    pop_obs[[pf]] <- array(dim = c(numF[[pf]], 4, 9))
    #Array for time steps fish was alive for each combination
    alive_obs[[pf]] <- array(dim = c(numF[[pf]], 9))     
    
   #NB: Fish size & sex  over time
   #Array of fish size across the 9 time steps  for each combination 
    Size[[pf]]<- array(dim = c(numF[[pf]], 9)) 
    Sex[[pf]]<- array(dim = c(numF[[pf]], 9))
    Parasite_strain[[pf]]<-  array(dim = c(numF[[pf]], 9))
    Fish_type[[pf]]<-  array(dim = c(numF[[pf]], 9))

    for (i in 1:numF[[pf]]) {
          pop_obs[[pf]][i,,] <- t(Parasite_fish_data[[parasite_fish[pf]]]
                         [Parasite_fish_data[[parasite_fish[pf]]]$Fish_ID==fishID[[pf]][i], 1:4])
          alive_obs[[pf]][i, ] <- ifelse(is.na(pop_obs[[pf]][i,1,]), 2, 1)
          Size[[pf]][i, ]<- Parasite_fish_data[[parasite_fish[pf]]][Parasite_fish_data[[parasite_fish[pf]]]$Fish_ID==fishID[[pf]][i], 9]
 
         #1=Female fish & 2=Male fish 
          Sex[[pf]][i, ]<-paste(Parasite_fish_data[[parasite_fish[pf]]][Parasite_fish_data[[parasite_fish[pf]]]$Fish_ID==fishID[[pf]][i], 10])    

          Parasite_strain[[pf]][i, ]<-paste(Parasite_fish_data[[parasite_fish[pf]]][Parasite_fish_data[[parasite_fish[pf]]]$Fish_ID==fishID[[pf]][i], 7])   

          Fish_type[[pf]][i ,]<-paste(Parasite_fish_data[[parasite_fish[pf]]][Parasite_fish_data[[parasite_fish[pf]]]$Fish_ID==fishID[[pf]][i], 8])   

                         }
# Experiment descriptors
fishSize[[pf]]<-apply(Size[[pf]],1,unique)
fishSex[[pf]]<- apply(Sex[[pf]],1,unique)
Strain[[pf]]<- apply(Parasite_strain[[pf]],1,unique)
Fish_stock[[pf]]<- apply(Fish_type[[pf]],1,unique)
                 }
    #return data on experiment descriptors (fish size, sex, fish type & strain) for each parasite-fish group
    return(list(fishSize=fishSize,fishSex=fishSex,Strain=Strain,Fish_stock=Fish_stock,
             numF= numF,fishID=fishID,pop_obs=pop_obs,alive_obs=alive_obs,fishID=fishID))
        }  

# Fraction of area of the 4 re-categorized body regions for male and female fish
#area of each body part (depends on size and gender)
Body_area<-function(Area_data){
#Area across the four body location
Area_data$Area_Pectoral<-Area_data$Area_Pectoral*2 #Note there are two pectoral fins on the guppy fish
Area_data$Area_Pelvic<-Area_data$Area_Pelvic*2 #Note there are two pelvic fins on the guppy fish also
Area_data$Area_LowerRegion<-Area_data$Area_LB+Area_data$Area_Anal+Area_data$Area_Pelvic+Area_data$Area_Dorsal
Area_data$Area_UpperRegion<-Area_data$Area_Pectoral+Area_data$Area_UB

Area_fishSex<- data.frame(matrix(NA,nrow=4,ncol=2))
names(Area_fishSex)<-c("Female_areas","Male_areas")
rownames(Area_fishSex)<-c("Tail","LowerRegion","UpperRegion","Head")

Area_fishSex[1,]<-as.vector(tapply(Area_data$Area_Tail,Area_data$Fish_sex,mean))
Area_fishSex[2,]<-as.vector(tapply(Area_data$Area_LowerRegion,Area_data$Fish_sex,mean))
Area_fishSex[3,]<-as.vector(tapply(Area_data$Area_UpperRegion,Area_data$Fish_sex,mean))
Area_fishSex[4,]<-as.vector(tapply(Area_data$Area_Head,Area_data$Fish_sex,mean))

#print(paste("average body area for each male and female fish"))
#Area_fishSex

Prop_area<-function(x)return(x/sum(x)) #normalizing the area to sum up to 1 for each fish sex

Area_normalized<- apply(Area_fishSex,2,Prop_area)

#print(paste("normalized across fish sex"))
#Area_normalized
#Areas_Female<-Area_normalized[,1]
#Areas_male<-Area_normalized[, 2]
    
    #Output:
    #Area_normalized[, 1]= normalized area for female fish
    #Area_normalized[, 2]= normalized area for male fish 
    return(Area_normalized)
}
