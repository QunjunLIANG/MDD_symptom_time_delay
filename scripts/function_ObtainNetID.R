#######################################
#
# This function is used to obtain the 
# network identification of Schafer 400

ObtainNetID <- function(id_path, anno_path){
  
  net_id <- readr::read_csv(id_path, col_names = 'network_ind') 
  net_parcel <- readr::read_csv(anno_path, col_names = 'region')
  net_ana <- cbind(net_id, net_parcel)
  net_ana$network_ind <- factor(net_ana$network_ind)
  net_ana['parcel'] <- paste0('parcel', 1:400)
  net_ana <- net_ana %>% mutate(network = ifelse(str_detect(region, pattern = 'Vis'),'visual',
                                                 ifelse(str_detect(region, pattern = 'SomMot'),'somMot',
                                                        ifelse(str_detect(region, pattern = 'DorsAttn'),'dorsalAttn',
                                                               ifelse(str_detect(region, pattern = 'SalVentAttn'),'salience',
                                                                      ifelse(str_detect(region, pattern = 'Limbic'),'limbic',
                                                                             ifelse(str_detect(region, pattern = 'Cont'),'control',
                                                                                    'DMN')))))))
  
  net_ana$network <- factor(net_ana$network, levels = c("somMot","visual",'dorsalAttn',"salience","control","limbic","DMN"))
  
  return(net_ana)
}