#-------------------------------------------------------------------------------#
#                       Model fitting                                           #
#-------------------------------------------------------------------------------#

#----   Loading of data packages  

pacman::p_load(sf,               
               lubridate,        
               tmap,    
               SpatialEpi,    
               ggplot2,   
               rgdal,     
               gganimate,     
               INLA,     
               spdep,     
               tidymodels,    
               readr,     
               utils,
               cvms,
               dplyr,
               klaR
               
)

#---- Load data files 

map = readOGR("data/mwi_admbnda_adm2_nso_20181016.shp") #Downloaded from the internet after a simple Google Search

plot(map)
map_sf <- st_as_sf(map)
map_sf = map_sf[,3]

#--- Assuming that there is cholera outbreak and one is interested to model the trends

#Simulate human population 

seqpop = seq(from =60000, to = 1000000, by = 500)
map_sf$pop = sample(seqpop, size =32, replace = TRUE)


#--- Simulate disease incidences

#data$cases*

seqcases2010 = seq(from = 3, to = 600, by = 3)
map_sf$cases2010 = sample(seqcases2010, size =32, replace = TRUE)

seqcases2011 = seq(from = 3, to = 600, by = 3)
map_sf$cases2011 = sample(seqcases2011, size =32, replace = TRUE)

seqcases2012 = seq(from = 3, to = 600, by = 3)
map_sf$cases2012 = sample(seqcases2012, size =32, replace = TRUE)

# Restructure the data

#map_sf <- gather(map_sf, date,SIR, 1:84)
map_sf2 <- gather(map_sf, date,cases, 4:6)

#-- replace SIR

#map_sf$date = gsub(pattern = "SIR.", "", map_sf$date)
map_sf2$date = gsub(pattern = "cases", "", map_sf2$date)

map_sf2$date = as.numeric(map_sf2$date)

#map_sf2$date = lubridate::y(map_sf2$date)


#---- Neighborhood matrix to define spatial random effect

#---- Neighborhood matrix to define spatial random effect

nb <- poly2nb(map)
malawi_graph = nb2mat(nb, style = 'B', zero.policy = TRUE)

#---- Modelling 

#--- Start with basic GLM

logit_base = glm(cases ~pop,
                 data =map_sf2, family = 'poisson')
summary(logit_base)

#--- INLA model here

map_sf2$idarea <- as.numeric(as.factor(map_sf2$ADM2_EN))
map_sf2$idate <- 1 + map_sf2$date - min(map_sf2$date)

#------- model fitting

res <- inla(cases ~ offset(log(pop))+
              f(idate, model = 'iid') +
              f(idarea, model ='besag', graph = malawi_graph, scale.model =  TRUE),
            family = "poisson", data = map_sf2, verbose = TRUE,
            control.compute=list(config = TRUE)
            )

#-------- Adding a categorical variable for modelling using INLA
#Simulating data with only one class
map_sf2$LC = 1
map_sf2$LC = as.factor(map_sf2$LC)

#------- model fitting

res2 <- inla(cases ~ offset(log(pop))+ LC +
              f(idate, model = 'iid') +
              f(idarea, model ='besag', graph = malawi_graph, scale.model =  TRUE),
            family = "poisson", data = map_sf2, verbose = TRUE,
            control.compute=list(config = TRUE)
)



