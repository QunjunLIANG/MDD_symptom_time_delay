#####################################
#
# This function used to generate surface 
# plot based on Schafer 400

# define a function for plotting
DrawSurfaceOnSchaefer400 <- function(value, 
                                     p_title, legend_title,
                                     legent_pos = 'bottom',
                                     boundary_width=0.5, 
                                     highColor = "darkred", midColor = 'white',
                                     lowColor = "steelblue") {
  
  library(ggseg)
  library(ggsegSchaefer)
  library(ggeasy)
  
  ## Load in atlas data provided by ggseg package
  atlas      = as_tibble(schaefer7_400)
  
  ## Select atlas region names and hemisphere so that we can add the values
  ## we want to plot:
  region     = atlas$region
  hemi       = atlas$hemi
  data       = distinct(na.omit(data.frame(region,hemi))) #remove NA and duplicate regions
  
  ## add the value to each region
  data['value'] <- value
  
  atlas_data <- atlas %>% left_join(data) # add projection map to the atlas
  
  p_this <- ggplot() + geom_brain(
    atlas       = atlas_data,
    mapping     = aes(fill=value),
    position    = position_brain("horizontal"),
    color       ='black',
    size        = boundary_width,
    show.legend = T) +
    ggtitle(p_title) +
    scale_fill_gradient2(low = lowColor, high = highColor, mid = midColor) +
    theme_void() + easy_move_legend(to = legent_pos) +
    easy_center_title() + easy_add_legend_title(legend_title) +
    easy_text_size(which = c("plot.title", "legend.title", "strip.text"),
                   size = 15)
  return(p_this)
}

DrawSurfaceOnSchaefer400_binary <- function(value, 
                                            p_title, legend_title,
                                            legent_pos = 'bottom',
                                            boundary_width=0.5, 
                                            highColor = "darkred", midColor = 'white',
                                            lowColor = "steelblue") {
  
  library(ggseg)
  library(ggsegSchaefer)
  library(ggeasy)
  
  ## Load in atlas data provided by ggseg package
  atlas      = as_tibble(schaefer7_400)
  
  ## Select atlas region names and hemisphere so that we can add the values
  ## we want to plot:
  region     = atlas$region
  hemi       = atlas$hemi
  data       = distinct(na.omit(data.frame(region,hemi))) #remove NA and duplicate regions
  
  ## add the value to each region
  data['value'] <- value
  
  atlas_data <- atlas %>% left_join(data) # add projection map to the atlas
  
  p_this <- ggplot() + geom_brain(
    atlas       = atlas_data,
    mapping     = aes(fill=value),
    position    = position_brain("horizontal"),
    color       ='black',
    size        = boundary_width,
    show.legend = T) +
    ggtitle(p_title) +
    scale_fill_gradient2(low = lowColor, high = highColor, mid = midColor) +
    theme_void() + easy_move_legend(to = legent_pos) +
    easy_center_title() + easy_add_legend_title(legend_title) +
    easy_text_size(which = c("plot.title", "legend.title", "strip.text"),
                   size = 15)
  return(p_this)
}

DrawComparableSurfaceOnSchaefer400 <- function(value_1, value_2, 
                                               value_1_name, value_2_name, 
                                               p_title, legend_title, 
                                               legent_pos = 'bottom',
                                               boundary_width=0.5, midColor = 'white',
                                               highColor = "darkred", lowColor = "steelblue") {
  
  library(ggseg)
  library(ggsegSchaefer)
  library(ggeasy)
  
  ## Load in atlas data provided by ggseg package
  atlas      = as_tibble(schaefer7_400)
  
  ## Select atlas region names and hemisphere so that we can add the values
  ## we want to plot:
  region     = atlas$region
  hemi       = atlas$hemi
  data       = distinct(na.omit(data.frame(region,hemi))) #remove NA and duplicate regions
  
  ## add the value to each region
  ## HAMA
  data_hama <- data
  data_hama['value'] <- value_1
  atlas_hama <- atlas %>% left_join(data_hama) # add projection map to the atlas
  atlas_hama['scale'] <-  value_1_name
  
  ## HAMD
  data_hamd <- data
  data_hamd['value'] <- value_2
  atlas_hamd <- atlas %>% left_join(data_hamd) # add projection map to the atlas
  atlas_hamd['scale'] <- value_2_name
  
  ## combine two data tables
  atlas_data <- rbind(atlas_hamd, atlas_hama)
  
  ## plot
  p_this <- ggplot() + geom_brain(
    atlas       = atlas_data,
    mapping     = aes(fill=value),
    position    = position_brain("horizontal"),
    color       ='black',
    size        = boundary_width,
    show.legend = T) +
    ggtitle(p_title) +
    facet_grid(scale~.) +
    #scale_fill_viridis_c(option= fill_color) +
    scale_fill_gradient2(low = lowColor, high = highColor, mid = midColor) +
    theme_void() + easy_move_legend(to = legent_pos) +
    easy_center_title() + easy_add_legend_title(legend_title) +
    easy_text_size(which = c("plot.title", "legend.title", "strip.text"),
                   size = 15)
  return(p_this)
}

DrawMultipleSurfaceOnSchaefer400 <- function(dat_value, 
                                             dat_name,
                                             p_title, legend_title, 
                                             legent_pos = 'bottom',
                                             boundary_width=0.5, 
                                             highColor = "darkred", 
                                             lowColor = "steelblue",
                                             midColor = 'white') {
  
  library(ggseg)
  library(ggsegSchaefer)
  library(ggeasy)
  
  ## Load in atlas data provided by ggseg package
  atlas      = as_tibble(schaefer7_400)
  
  ## Select atlas region names and hemisphere so that we can add the values
  ## we want to plot:
  region     = atlas$region
  hemi       = atlas$hemi
  data       = distinct(na.omit(data.frame(region,hemi))) #remove NA and duplicate regions
  
  # tidy the data ----------------------------------------------------------------------
  
  ## firstly, generate dataframe for value 1
  data_1 <- data
  data_1['value'] <- dat_value[,1] 
  atlas_data <- atlas %>% left_join(data_1) # add projection map to the atlas
  atlas_data['scale'] <-  dat_name[1]
  
  ## add the data in loop
  for (i in 2:ncol(dat_value)) {
    data_tmp <- data
    data_tmp['value'] <- dat_value[,i] 
    atlas_tmp <- atlas %>% left_join(data_tmp) # add projection map to the atlas
    atlas_tmp['scale'] <-  dat_name[i]
    ## merge the atlas data
    atlas_data <- rbind(atlas_data, atlas_tmp)
  }
  
  ## plot
  p_this <- ggplot() + geom_brain(
    atlas       = atlas_data,
    mapping     = aes(fill=value),
    position    = position_brain("horizontal"),
    color       ='black',
    size        = boundary_width,
    show.legend = T) +
    ggtitle(p_title) +
    facet_grid(scale~.) +
    #scale_fill_viridis_c(option= fill_color) +
    scale_fill_gradient2(low = lowColor, high = highColor, mid = midColor) +
    theme_void() + easy_move_legend(to = legent_pos) +
    easy_center_title() + easy_add_legend_title(legend_title) +
    easy_text_size(which = c("plot.title", "legend.title", "strip.text"),
                   size = 15)
  return(p_this)
}
