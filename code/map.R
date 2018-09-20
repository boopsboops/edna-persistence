# load ggmap and ggplot
library(tidyverse)
library(ggmap)
library(OpenStreetMap)
library(ggsn)
library(ggplot2)


# read gps data
gps <- data.frame(site_code=c("E1", "harbour"), latitude=c(50.033,50.370), longitude=c(-4.367,-4.133)) 


# lat is flat
df <- projectMercator(lat=c(50.033,50.370), long=c(-4.367,-4.133))
df <- data.frame(site_code=c("E1", "harbour"), longitude=df[,1], latitude=df[,2]) 

# get a UK map
uk <- get_map(location=c(-3.7,51.5), zoom=6, maptype="toner-lite", source="stamen") 

# devon map
dev <- get_map(location=c(-4.3,50.25), zoom=10, maptype="toner-lite", source="stamen") 

# get the bounding box 
rec <- as_tibble(attr(dev,"bb"))

# plot
uk.plot <- uk %>% #
    ggmap + #
    geom_rect(xmin=rec$ll.lon, xmax=rec$ur.lon, ymin=rec$ll.lat, ymax=rec$ur.lat, col="red", fill=NA) +
    xlab("Longitude") + #
    ylab("Latitude") +
    theme_bw()

dev.plot <- dev %>% #
    ggmap + #
    geom_point(data=gps, aes(x = longitude, y=latitude), color="red", shape=3, alpha=1, size=5, stroke=1.25) + #
    ggsn::scalebar(dist=10, dd2km=TRUE, model="WGS84", y.min=50, y.max=50.3, x.min=-4.7, x.max=-3.9, st.size=4, location="bottomleft") +
    xlab("Longitude") + #
    ylab("Latitude") +
    theme_bw()

# save 
#ggsave(filename="../temp/uk.pdf", width=7, height=7, units="in", plot=uk.plot, dpi=600, device="pdf")
#ggsave(filename="../temp/devon.pdf", width=7, height=7, units="in", plot=dev.plot, dpi=600, device="pdf")
