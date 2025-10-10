options(warn=-1)
options(dplyr.summarise.inform = FALSE)

`%notin%` <- Negate(`%in%`)

# Function to allow rbinding dataframes with foreach even when some dataframes 
# may not have any rows
foreach_rbind <- function(d1, d2) {
  if (is.null(d1) & is.null(d2)) {
    return(NULL)
  } else if (!is.null(d1) & is.null(d2)) {
    return(d1)
  } else if (is.null(d1) & !is.null(d2)) {
    return(d2)
  } else  {
    return(rbind(d1, d2))
  }
}

#--------------------------------------------------------------------------------
match_wocat <- function(df, pid) {
  
  registerDoParallel(4)
  
  options("optmatch_max_problem_size"=Inf)
    
  # Note custom combine to handle iterations that don't return any value
  #test nested foreach loops
  ret <- foreach (this_lc=unique(df$land_cover),
                  .packages=c('optmatch', 'dplyr'),
                  .combine=foreach_rbind, .inorder=FALSE) %dopar% {
                    this_d <- df
                    d_wocat <- filter(this_d, status)
                    # Filter out climates and land covers that don't appear in the wocat
                    # sample, and drop these levels from the factors
                    this_d <- filter(this_d,
                                     land_cover %in% unique(d_wocat$land_cover))
##                                     wwfbiom %in% unique(d_wocat$wwfbiom),
##                                     wwfecoreg %in% unique(d_wocat$wwfecoreg)
                    
                    this_d$land_cover <- droplevels(this_d$land_cover)
##                    this_d$wwfbiom <- droplevels(this_d$wwfbiom)
##                    this_d$wwfecoreg <- droplevels(this_d$wwfecoreg)
                    # table(this_d$status)
                    dat <- dplyr::select(this_d, lat, lon, UID, status, land_cover,
                                         #wwfbiom, wwfecoreg, 
                                         elevation, slope, mean_temp, max_temp, min_temp, prec,
                                         d2road, d2city, popden, tt2city, popcnt) 
                    ps_A <- glm(status ~ mean_temp + max_temp + min_temp + prec + elevation + slope + d2road + d2city + popden + popcnt + tt2city, data = dat)
                    dat$propensity_scoreA <- fitted(ps_A)
                    f <- status~ mean_temp + max_temp + min_temp + prec + elevation + slope + d2road + d2city + popden + popcnt + tt2city
                    this_d <- bind_cols(this_d,propensity_scoreA=dat$propensity_scoreA)
                    # Can't stratify by land cover or climate if they only have one level
                    if (nlevels(this_d$land_cover) >= 2) {
                      f <- update(f, ~ . + strata(land_cover))
                    } else {
                      f <- update(f, ~ . - land_cover)
                    }
##                    if (nlevels(this_d$wwfbiom) >= 2) {
##                      f <- update(f, ~ . + strata(wwfbiom))
##                    } else {
##                      f <- update(f, ~ . - wwfbiom)
##                    }
##                    if (nlevels(this_d$wwfecoreg) >= 2) {
##                      f <- update(f, ~ . + strata(wwfecoreg))
##                    } else {
##                      f <- update(f, ~ . - wwfecoreg)
##                    }
                    if (nrow(d_wocat) > 2) {
                      model <- glm(f, data=this_d)
                      dists <- match_on(model, data=this_d)  #***NEED PROPENSITY SCORE BEFORE HERE
                    } else {
                      # Use Mahalanobis distance if there aren't enough points to run a glm
                      dists <- match_on(f, data=this_d)
                    }
                    #potentially drop caliper line; will cut down dists matrix but not the speed issue
                    # dists <- caliper(dists, 2)
                    # If the controls are too far from the treatments (due to the caliper) 
                    # then the matching may fail. Can test for this by seeing if subdim 
                    # runs successfully
                    subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                                             error=function(e)return(FALSE))
                    if (subdim_works) {
                      m <- fullmatch(dists, min.controls=1, max.controls=1, data=this_d)
                      prematch_d <- this_d
                      this_d$matched <- m
                      this_d <- this_d[matched(m), ]
                    } else {
                      this_d <- data.frame()
                    }
                    # Need to handle the possibility that there were no matches for this 
                    # treatment, meaning this_d will be an empty data.frame
                    if (nrow(this_d) == 0) {   #if matching w/ ecoreg return no results, match again w/o ecoreg and check matching results
                      this_d<-df
                      d_wocat <- filter(this_d, status)
                      this_d <- filter(this_d,
                                       land_cover %in% unique(d_wocat$land_cover))
##                                       wwfbiom %in% unique(d_wocat$wwfbiom))
                      
                      this_d$land_cover <- droplevels(this_d$land_cover)
##                      this_d$wwfbiom <- droplevels(this_d$wwfbiom)
                      f <- status ~ mean_temp + max_temp + min_temp + prec + elevation + slope + d2road + d2city + popden + popcnt + tt2city
                      if (nlevels(this_d$land_cover) >= 2) {
                        f <- update(f, ~ . + strata(land_cover))
                      } else {
                        f <- update(f, ~ . - land_cover)
                      }
##                      if (nlevels(this_d$wwfbiom) >= 2) {
##                        f <- update(f, ~ . + strata(wwfbiom))
##                      } else {
##                        f <- update(f, ~ . - wwfbiom)
##                      }
                      if (nrow(d_wocat) > 2) {
                        model <- glm(f, data=this_d)
                        dists <- match_on(model, data=this_d)
                      } else {
                        dists <- match_on(f, data=this_d)
                      }
                      subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                                               error=function(e)return(FALSE))
                      if (subdim_works) {
                        m <- fullmatch(dists, min.controls=1, max.controls=1, data=this_d)
                        prematch_d <- this_d
                        this_d$matched <- m
                        this_d <- this_d[matched(m), ]
                      } else {
                        this_d <- data.frame()
                      }
                      if(nrow(this_d)==0){
                        log_mes <- paste(pid,"-Matching without wwfecoreg:Failed\n",sep="")
                        dir.create(paste(paste(f.path,"WDPA_matching_log/",iso3,sep="")))
                        cat(log_mes,file=paste(f.path,"WDPA_matching_log/",iso3,"/",iso3,"_pa_",pid,"_matching_used_covar_log_wk", gediwk,".txt",sep=""),append=TRUE)
                        return(NULL)
                      } else{
                        log_mes <- paste(pid,"-Matching without wwfecoreg:Succeed\n",sep="")
                        dir.create(paste(paste(f.path,"WDPA_matching_log/",iso3,sep="")))
                        cat(log_mes,file=paste(f.path,"WDPA_matching_log/",iso3,"/",iso3,"_pa_",pid,"_matching_used_covar_log_wk", gediwk,".txt",sep=""),append=TRUE)
                        match_results <- list("match_obj" = m, "df" = this_d, "func"=f, "prematch_d"=prematch_d)
                        return(match_results)
                      }
                    } else {
                      log_mes <- paste(pid,"-Matching with wwfecoreg:Succeed\n",sep="")
                      match_results <- list("match_obj" = m, "df" = this_d, "func"=f, "prematch_d"=prematch_d)
                      dir.create(paste(paste(f.path,"WDPA_matching_log/",iso3,sep="")))
                      cat(log_mes,file=paste(f.path,"WDPA_matching_log/",iso3,"/",iso3,"_pa_",pid,"_matching_used_covar_log_wk", gediwk,".txt",sep=""),append=TRUE)
                      return(match_results)
                    }
                  }
  
  stopImplicitCluster()
  return(ret)
}

#--------------------------------------------------------------------------------
propensity_filter <- function(pa_df, d_control_local){
  pa_df <-pa_df[complete.cases(pa_df), ]  #filter away non-complete cases w/ NA in control set
  d <- dplyr::bind_rows(d_control_local, pa_df)
  ## bring in matching algorithm from STEP5 here to loop through each PA in d_PAs
  #filter controls based on propensity scores 
  d_all <- dplyr::select(d, lat, lon, UID, status, land_cover, 
                         #wwfbiom, wwfecoreg, 
                         elevation, slope, mean_temp, max_temp, min_temp, prec, 
                         d2road, d2city,  popden, tt2city, popcnt) 
  
  d_all$status <- ifelse(d_all$status==TRUE,1,0)
  
  #calculate the propensity scores & filter out controls not overlapping w/ treatment propensity scores
  ps <- glm(status ~ mean_temp + max_temp + min_temp + prec + elevation + slope + d2road + d2city + popden + popcnt + tt2city,data = d_all)
  # boxplot(ps)  #check the distribution of propensity scores for treatment and controls
  #filter out the controls with propensity scores outside of the overlapping region
  d_all$propensity_score <- fitted(ps)
  d_sep <- d_all %>% dplyr::group_by(status)
  d_sep_range <- d_all %>% dplyr::group_by(status)%>% 
    dplyr::summarise(propmin= min(propensity_score), promax=max(propensity_score)) 
  # cat(iso3, "Filtering the control sites by overlaping with the treatment PS\n")
  d_filtered <- d_sep %>% 
    filter(status==1 | between(propensity_score,d_sep_range$propmin[2],d_sep_range$promax[2])) %>% 
    ungroup() 
  
  d_filtered$status <- ifelse(d_filtered$status==1,TRUE,FALSE)
  
  return(d_filtered)
}

#--------------------------------------------------------------------------------
matched2ras <- function(matched_df){
  cat(iso3,"converting the matched csv to a raster stack for extraction\n")
  
  matched_pts <- SpatialPointsDataFrame(coords=matched_df[,c("lon","lat")],
                                        proj4string=CRS("+init=epsg:4326"), data=matched_df) %>% 
    spTransform(., CRS("+init=epsg:6933"))

  matched_pts <- vect(matched_pts)
  
  # Ensure fields are in appropriate formats
  matched_pts$UID <- as.integer(matched_pts$UID)
  matched_pts$pa_id <- as.integer(matched_pts$pa_id)
  #matched_pts$pa_id <- as.integer(rep(id_pa,nrow(matched_pts)))
  matched_pts$status <- as.logical(matched_pts$status)
##  matched_pts$wwfbiom <- as.numeric(matched_pts$wwfbiom)
##  matched_pts$wwfecoreg <- as.numeric(matched_pts$wwfecoreg)
 
  # Define the raster extent for cropping
  buffer_ext <- ext(buffer(matched_pts, 10000))
  r <- crop(MCD12Q1, buffer_ext)
  continent <- crop(world_region, buffer_ext)
  
  # Assign names to the raster layers
  names(r) <- "pft"
  names(continent) <- "region"
  
  # List of fields to rasterize
  fields <- c("status", "pa_id", "UID", "land_cover") #"wwfbiom", "wwfecoreg", 
  rasters <- list()
  
  # Rasterize each field
  for (field in fields) {
    r_field <- rasterize(matched_pts, r, field = field)
    rasters[[field]] <- r_field
  }
  
  # Combine rasters into a SpatRaster stack
  rasters <- c(rasters, list(pft = r, region = continent))
  matched_ras <- rast(rasters)
  
  return(matched_ras)
}

#--------------------------------------------------------------------------------
convertFactor <- function(matched0, exgedi){
  exgedi$pft <- as.character(exgedi$pft)
  
  exgedi$pft <- factor(exgedi$pft, levels=sequence(6),
                       labels = c("ENT",
                                  "EBT",
                                  "ENT",
                                  "DBT",
                                  "GS",
                                  "GS"))
  exgedi$region <- as.character(exgedi$region)
  exgedi$region <- factor(exgedi$region, levels=c(1:7),
                          labels = c("Eu",
                                     "As",
                                     "Au",
                                     "Af",
                                     "As",
                                     "SA",
                                     "US"))
  
  exgedi$stratum <- paste(exgedi$pft, exgedi$region,sep="_")
    
  return(exgedi)
}

#--------------------------------------------------------------------------------
subdfExport <- function(filtered_df){
  #export invidual pa results
  spt2 <- split(filtered_df, filtered_df$pa_id)
  
  dfl <- lapply(names(spt2), function(x){
    
    if(dim(spt2[[x]])[1]>0){
      control_sub <- spt2[[x]][spt2[[x]]$status==0,]
      treat_sub <-  spt2[[x]][spt2[[x]]$status==1,]
      ncontrol <- nrow(control_sub)
      ntreat <- nrow(treat_sub)
      if (ntreat-ncontrol > 0){
        newtreatid <- sample(ntreat, ncontrol)
        newtreat <- treat_sub[newtreatid,]
        spt2_new <- rbind(newtreat, control_sub)
      } else if (ntreat - ncontrol< 0){
        newcontrolid <- sample(ncontrol, ntreat)
        newcontrol <- control_sub[newcontrolid,]
        spt2_new <- rbind(newcontrol, treat_sub)
      } else if (ntreat-ncontrol==0){
        spt2_new <- spt2[[x]]
      } else if (ntreat==0 || ncontrol==0){
        spt2_new=NA
      }
##      biom <- spt2_new$wwfbiom %>% unique() %>% as.character() %>% gsub('\\b(\\pL)\\pL{4,}|.','\\U\\1',.,perl = TRUE)
##      if(length(biom)>1){
##        biom <- paste(c(biom), collapse="&")
##      }
      # print(biom)
##      write.csv(spt2_new, file=paste(f.path,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",iso3,"_PA_",unique(spt2_new$pa_id),"_",biom,".csv", sep=""))
        write.csv(spt2_new, file=paste(f.path,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",iso3,"_PA_",unique(spt2_new$pa_id),".csv", sep=""))

      return(spt2_new)
    }
  })
  
  total_df <- do.call("rbind", dfl) 
  cat("Exported individual PAs results for ", iso3, "\n")
  
  return(total_df)
}

getmode <- function(v,na.rm) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#--------------------------------------------------------------------------------
extract_gedi <- function(matched, mras, iso3){
  #f.path <- '~/GEDI_PA/Matching_Layers/SEN/SEN_Tiles/'
  #f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"
  f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
  f.path_l2 <- paste(f.path,"WDPA_gedi_L2A_tiles/",sep="")
  f.path_l4 <- paste(f.path,"WDPA_gedi_L4A_tiles/",sep="")

######Vero## select a subset of tiles that overlap with the points from matched_df    
    s3_get_files(c(paste(f.path,"GRID1deg_poly_52N52S/GRID1deg_poly_52N52S.shp",sep=""),
              paste(f.path,"GRID1deg_poly_52N52S/GRID1deg_poly_52N52S.shx",sep=""),
              paste(f.path,"GRID1deg_poly_52N52S/GRID1deg_poly_52N52S.prj",sep=""),
              paste(f.path,"GRID1deg_poly_52N52S/GRID1deg_poly_52N52S.dbf",sep="")),confirm = FALSE)
    
    #G1deg_tiles <- terra::vect(s3_get(paste(f.path,"GRID1deg_poly_52N52S/GRID1deg_poly_52N52S.shp",sep="")))
    G1deg_tiles <- st_read(s3_get(paste(f.path,"GRID1deg_poly_52N52S/GRID1deg_poly_52N52S.shp",sep="")))
    G1deg_tiles <- vect(G1deg_tiles)
    matched_points <- vect(matched, geom=c("lon","lat"), crs="epsg:4326", keepgeom=FALSE)       

    intersecting <- G1deg_tiles[matched_points]
        if (nrow(intersecting) == 0) {stop('no intersecting G1deg_tiles found')
            } else {
            tileindex <- intersecting$id
            }
    #print(length(tileindex))

    tileindex_df <- read.csv(s3_get(paste(f.path,"vero_1deg_tileindex/tileindex_",iso3,".csv", sep="")))
    iso3_tiles <- intersect(tileindex_df$tileindex, tileindex)

    print(paste("number of overlapping GEDI tiles =",length(iso3_tiles),sep=" "))
######Vero## select a subset of tiles that overlap with the points from matched_df    
    
  all_gedil2_f <- c()
  all_gedil4_f <- c()
  for(i in 1:length(iso3_tiles)){
    iso3_tile_in <- paste("tile_num_",iso3_tiles[i],sep="")
    gedi_file_l2 <- paste(iso3_tile_in,"_L2A.gpkg",sep="")
    gedi_file_l4 <- paste(iso3_tile_in,"_L4A.gpkg",sep="")
    all_gedil2_f <- c(all_gedil2_f, gedi_file_l2)
    all_gedil4_f <- c(all_gedil4_f, gedi_file_l4)
    }
  
  # Initialize an empty list to store results
  results_list <- list()
#  iso_matched_gedi_df <- NULL # Initialize before loop

  # Iterate over the sequence of indices for your files
  for (this_tile in seq_along(all_gedil2_f)) {
                cat("Reading in no. ", this_tile, "GPKG file of ", length(all_gedil2_f), "GEDI tiles for iso3", iso3, "\n")
                
                #f.path <- "/projects/my-public-bucket/GEDI_global_PA_v2/"
                f.path <- "s3://maap-ops-workspace/shared/leitoldv/GEDI_global_PA_v2/"
      
                # Read GEDI L4A data
                gedil4_f_path <- paste(f.path, "WDPA_gedi_L4A_tiles/", all_gedil4_f[this_tile], sep = "")
                gedil4_f <- vect(s3_get(gedil4_f_path))
                
                # Read GEDI L2A data
                gedil2_f_path <- paste(f.path, "WDPA_gedi_L2A_tiles/", all_gedil2_f[this_tile], sep = "")
                gedil2_f <- vect(s3_get(gedil2_f_path))

                gedi_l2_sub <- gedil2_f[,c("shot_number", "lat_lowestmode", "lon_lowestmode", "rh98", "filename")] # , "rh25", "rh50", "rh75", 
                
                rm(gedil2_f)
      
                # Check if GEDI L4A data is empty
                if (nrow(gedil4_f) < 1) {
                    cat("Error: No data for GEDI L4A\n")
                    #gedi_l24 <- gedil2_f
                    gedi_l24 <- gedi_l2_sub
                    gedi_l24$agbd <- NA
                    gedi_l24$agbd_se <- NA
                    gedi_l24$agbd_t <- NA
                    gedi_l24$agbd_t_se <- NA
                } else {
                    # Select relevant columns from GEDI L4A
                    gedi_l4_sub <- gedil4_f[, c("shot_number", "agbd", "agbd_se")]
                    
                    rm(gedil4_f)

                    # Join with GEDI L2A data
                    gedi_l24 <- merge(gedi_l2_sub, gedi_l4_sub, by = "shot_number")
                    
                    gedi_l24$year <- as.integer(sub(".*GEDI04_A_(\\d{4}).*", "\\1", gedi_l24$filename))
                    gedi_l24$filename <- NULL  # Drop filename column
                
                }
            
                cat("GEDI tile # ", this_tile, "GPKG file has", nrow(gedi_l24), "rows and ", ncol(gedi_l24), "columns", "\n")
                
                rm(gedi_l2_sub)
                rm(gedi_l4_sub)
            
                # Initialize empty spatial object for the current iteration
                gedi_l24_sp <- NULL
            
                # Convert to spatial points data frame if there is data
                if (nrow(gedi_l24) > 0) {
                gedi_l24_sp <- gedi_l24

                rm(gedi_l24)
                
                gedi_l24_sp <- project(gedi_l24_sp, "epsg:6933")
                
                #matched_gedi <- terra::extract(mras,vect(gedi_l24_sp), df=TRUE)
                #matched_gedi_metrics <- cbind(matched_gedi, gedi_l24_sp@data)
                matched_gedi <- terra::extract(mras, gedi_l24_sp, df=TRUE)
                matched_gedi_metrics <- cbind(matched_gedi, gedi_l24_sp)
                print(nrow(matched_gedi_metrics))
                matched_gedi_metrics_filtered <- matched_gedi_metrics %>% dplyr::filter(!is.na(status)) %>%
                                                    convertFactor(matched0 = matched,exgedi = .)
                print(nrow(matched_gedi_metrics_filtered))

            #iso_matched_gedi_df <- rbind(matched_gedi_metrics_filtered, iso_matched_gedi_df)
            cat("dataframe #", this_tile, "has", nrow(matched_gedi_metrics_filtered), "rows and ", ncol(matched_gedi_metrics_filtered), "columns", "\n")

            rm(gedi_l24_sp)
            rm(matched_gedi)
            rm(matched_gedi_metrics)
            #rm(matched_gedi_metrics_filtered)
         }
        
        # Store results in a list
        #results_list[[this_tile]] <- iso_matched_gedi_df
        results_list[[this_tile]] <- matched_gedi_metrics_filtered
        #print(nrow(results_list[[this_tile]]))
    }
    
    # Combine all results
    if (!is.null(results_list)) {
#    if (!is.null(iso_matched_gedi_df)) {
        #iso_matched_gedi_df <- do.call(rbind, results_list)
        iso_matched_gedi_df <- dplyr::bind_rows(results_list)
        cat("output dataframe has", nrow(iso_matched_gedi_df), "rows and ", ncol(iso_matched_gedi_df), "columns", "\n")
    }
    
    cat("Done GEDI processing for PA ", id_pa, "\n")
    return(iso_matched_gedi_df)
}

#--------------------------------------------------------------------------------
SplitRas <- function(raster,ppside){
  h        <- ceiling(ncol(raster)/ppside)
  v        <- ceiling(nrow(raster)/ppside)
  agg      <- aggregate(raster,fact=c(h,v))
  agg[]    <- 1:ncell(agg)
  agg_poly <- rasterToPolygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for(i in 1:ncell(agg)){
    e1          <- extent(agg_poly[agg_poly$polis==i,])
    r_list[[i]] <- crop(raster,e1)
  }
  return(r_list)
}

#--------------------------------------------------------------------------------
rasExtract2020 <- function(l4_sp){
  cat(iso3,"converting the matched csv to a raster stack for extraction\n")
  tif2020 <- c("pop_cnt_2020", "pop_den_2020", "lc2019", "tt2cities_2015", 
               "wc_prec_2010-2018", "wc_tavg_2010-2018", "wc_tmax_2010-2018",
               "wc_tmin_2010-2018","dem","slope","d2roads","dcities") #"wwf_biomes","wwf_ecoreg",
  for (t in 1:length(tif2020)){
    # print(tif2020[t])
    covar2020 <- raster(paste(f.path, "WDPA_input_vars_iso3_v2/",iso3,"/",tif2020[t],".tif", sep=""))
    ras_ex <- raster::extract(covar2020, l4_sp@coords, method="simple", factors=F)
    nm <- names(covar2020)
    l4_sp <- cbind(l4_sp, ras_ex)
    names(l4_sp)[t+6] <- tif2020[t]
  }
  return(l4_sp)
}

