
#import the libraies I need to use
library(tidyverse) #data reading, wrangling, and plotting
library(rjson) # dealing with the JSON data format the downloads come in
library(lubridate) #dealing with dates


#the location of all the source data for the different locations around Australia
#https://data.gov.au/data/organization/australian-radiation-protection-and-nuclear-safety-agency-arpansa
source <- c("https://data.gov.au/api/3/action/package_show?id=cedcfdbb-d360-4a99-8369-d5941c2ad66b", #Townsville
            "https://data.gov.au/api/3/action/package_show?id=c31a759c-a4d4-455f-87a7-98576be14f11", #sydney
            "https://data.gov.au/api/3/action/package_show?id=1b55352e-c0d8-48c8-9828-ef12885c9797", #Perth
            "https://data.gov.au/api/3/action/package_show?id=18c7e503-0a11-43c7-a401-a983d5e17a3b", #newcastle
            "https://data.gov.au/api/3/action/package_show?id=fb836013-f300-4f92-aa1e-fb5014aea40e", #melbourne
            "https://data.gov.au/api/3/action/package_show?id=a182fb57-4355-404d-b462-dfd7c92feba9", #kingston
            "https://data.gov.au/api/3/action/package_show?id=db5bf171-3367-4c07-91f8-0d3c2ac62e98", #gold coast
            "https://data.gov.au/api/3/action/package_show?id=f4fd0683-7fdc-46ee-a83b-b1ed9099fb09", #Emerald
            "https://data.gov.au/api/3/action/package_show?id=d549770d-6349-405a-b6f0-c168e4edabb8", #Darwin
            "https://data.gov.au/api/3/action/package_show?id=154d4d3b-2e8d-4dc2-b8ac-8f9805f99826", #Canberra
            "https://data.gov.au/api/3/action/package_show?id=2a1a2e49-de97-450e-9d0a-482adec68b22", #Brisbane
            "https://data.gov.au/api/3/action/package_show?id=532e300b-89a7-4b37-9ec1-f0476ccdce94", #Alice Springs
            "https://data.gov.au/api/3/action/package_show?id=026d4974-9efb-403d-9b39-27aee31a6439"  #Adelaide
            )

#Let's get them all
js <- source %>%  map(~fromJSON(file = .))


##The original loop I managed to make into a partial map
#for(i in 1:length(js)) {
#  for(j in 1: length(js[[i]]$result$resources)){
#   urls <- append(urls, js[[i]]$result$resources[[j]]$url ) 
#  }
#}



#extract urls from the JSON objects so I can get the actual data
#half way there with map
urls <- c()
for(i in 1:length(js)) {
  urls <- append(urls, map_chr(js[[i]]$result$resources, "url") ) 
}


#read in data
#about 2 Gb in this step
df_master <- urls %>% map(read_csv)

#save the object so I can do things again without having to download all of it again.
saveRDS(df_master, file = "data.rds")
#df_master <- readRDS("data.rds")



#clean column names and fix UV index as a numeric
df_master <- df_master %>% map(~rename(., 
                                       timestamp = 1,
                                       lat = 2,
                                       lon = 3,
                                       UV_index = 4))

df_master <- df_master %>% map(~mutate(., UV_index = as.numeric(UV_index)))



#get the location names from the urls
url_tokens <- unlist(str_split(urls, "/"))
locations <- url_tokens[seq(10, length(url_tokens), 10)]
locations <- str_match(locations, "uv-(.*)-\\d\\d\\d\\d.*")[,2]


#there is probably a map to do this
#put the locations against the data
for(i in 1:length(df_master)) {
  df_master[[i]] <- df_master[[i]] %>%  mutate(location = locations[i])
}
#


#make the list into a single tibble
df <- bind_rows(df_master)

#fee up some memory
rm(df_master)

#add some helper columns
df <- df %>% 
  mutate(hour = hour(timestamp) + (minute(timestamp)/60) + (second(timestamp)/60/60),
         date = as.Date(timestamp))

#save the data again for later
saveRDS(df, "df.RDS")




#now to nest data by location and then date. 
#  each date holds all the data for a single day and is 
#  then fitted by a single instance of an nls fit.

df_nested <- 
df %>% 
  group_by(location, date) %>% 
  filter(n()>50) %>% #make sure there is more than 50 observations in a day. 
  nest() %>% 
  ungroup()
  
#to save memory
rm(df)

#to pick up again for later
saveRDS(df_nested, "df_nested.RDS")





#peak shapes
#I tried Gaussian and pseudo-Voigt, but Gaussian was the best.
# Here, "int" is the peak area, "posn" the peak maximum, and 
# "fwhm" is the full-width at half-maximum
gaussian <- function (int, posn, fwhm, x)	{
  int*(2*sqrt(log(2)/pi)/fwhm)*exp(-4*log(2)*((x-posn)/fwhm)^2)
}

lorentzian <- function ( int, posn, fwhm, x)	{
  int*(1/(2*pi))*(fwhm / ((fwhm/2)^2 + (x-posn)^2))
}

pseudo_voigt <- function(int, posn, fwhm, mix, x)	{
  mix*gaussian(int, posn, fwhm, x) + (1 - mix)*lorentzian(int, posn, fwhm,x)
}



#I'm going to loop over everything. There is probably a map way of doing it, but right
# now I want to get the modelling done.

##when you need to debug...
#df_nested <- readRDS("df_nested.RDS")

#store the coefficients of the models in this vector.
#  preallocate it to save copying everything everytime you touch it.
coefs <- vector(mode = "list", length = length(df_nested[[3]]))

#for all of the data
# this takes a long time. I don't know how to map it, but it would probably speed it up a lot.
for(i in 1:length(df_nested[[3]])){
  cat('Processing line', i, 'of', length(df_nested[[3]]),'\n')
  
  #do a fit
  fit <- nls(UV_index ~ gaussian(int, posn, fwhm, hour), data=df_nested[[3]][i][[1]],
             start=list(int=3, posn=12, fwhm=1), 
             lower=list(int=0, posn= 6, fwhm=0.01),
             upper=list(int=100,posn=18,fwhm=10),
             control = nls.control(maxiter = 1000,
                                   warnOnly = TRUE),
             algorithm="port")  
  
  #get the coefficients
  coefs[[i]] <- coef(fit)
  
  #   add the fit to the data
  df_nested[[3]][i][[1]] <- df_nested[[3]][i][[1]] %>% 
    mutate(fit = predict(fit)) 
}

#check them all out!
warnings()

#I want to get just the date, location, and peak parameters to visualise.
df_summary <- df_nested %>% 
  select(location, date) %>% 
  mutate(fit = coefs) %>% 
  unnest_wider(fit) %>% 
  mutate(height = int*(2*sqrt(log(2)/pi)/fwhm)) #calculate peak height from its integrated area and width




#make that plot
df_summary %>% 
  mutate(year = year(date)) %>% 
  filter(!location %in% c("emerald","gold-coast", "kingston","newcastle"),
         year >= 2015) %>% 
  mutate(location = str_replace(location, "-"," "),
         location = str_to_title(location)) %>% 
  ggplot(aes(x = date, y = height, colour = posn)) +
  geom_line() +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_continuous(type="viridis", 
                         limits = c(11,14), 
                         breaks=c(11,12,13,14), 
                         labels=c("11:00","12:00", "13:00","14:00"))+
  coord_cartesian(ylim = c(0,16))+
  facet_wrap(~ location)+
  labs(x = "Date of observation",
       y = "Peak height of fitted Gaussian (UV index)",
       colour = "Peak time",
       title = "Variation in daily maximum UV index in some Australian cities",
       subtitle = "The peak height, and time of that height, were taken from a Gaussian fit to each day's data",
       caption = "Data from the Australian Radiation Protection and Nuclear Safety Agency - https://data.gov.au/\nModelling and visualisation by Matthew Rowles | @xray_matt") +
  theme_dark()+
  theme(panel.background = element_rect(fill = "grey55"),
        plot.background = element_rect(fill = "grey25"),
        legend.background = element_rect(fill = "grey25"),
        axis.text = element_text(colour="grey75"),
        axis.title = element_text(colour="grey75"),
        plot.title = element_text(colour="grey95"),
        plot.subtitle = element_text(colour="grey95"),
        plot.caption = element_text(colour="grey95"),
        legend.text = element_text(colour="grey75"),
        legend.title = element_text(colour="grey75"),
        #legend.position = "bottom"
        )




