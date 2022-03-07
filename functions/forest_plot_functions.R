# get correlates of protection
gilbert_cp <- readxl::read_excel("./data/gilbert21_mod_correlate_protection.xlsx")
feng_cp <- readxl::read_excel("./data/feng21_chAd_correlate_protection.xlsx")

common_cols <- intersect(colnames(gilbert_cp), colnames(feng_cp))

cps <- rbind(gilbert_cp[,common_cols], feng_cp[,common_cols] %>% filter(measurement != "NF50 (Live-virus)"))
cps$Log2Titer <- log2(cps$Titer/10)
cps <- cps %>% filter(Titer > 5)

x_labels <- 2^seq(from = 0, to = 15, by = 2)*10


# point size. original: 3, hline = 1.5
figure_factor_increase <-  1.6 # for 3x vax only 0.8
point_size <-1.2/figure_factor_increase
axis_text_size_p <- (4.5/figure_factor_increase)
hline_size <- 2/figure_factor_increase
arrow_line_size <- 0.5/figure_factor_increase
arrow_head_size <- 0.05/figure_factor_increase

mapColors <- read.csv(file = './data/map/ag_colors.csv', row.names = 'Antigen', header = TRUE)

colors <- c(mapColors[c("B.1.617.2","B.1.1.7","B.1.351", "WA1", "WA1", "WA1", "WA1", "WA1"),"Color"], "grey", "#BB0D3A")
names(colors) <- c("Delta", "Alpha", "Beta", "Wu-1", "D614G","WT", "B.1","WA1", "Wu-1?", "Omicron")
fill_colors <- c("#BB0D3A", "#F2ACB9")
names(fill_colors) <- c("n", "y")

shapes <- c(21:25, 0:3, 5, 6)
shapes_fold <- c(21:25, 15:20)
names(shapes) <- c("VeroE6", "VeroE6-TMPRSS2", "Unknown", "Vero","HEK293T-ACE2",
               "H1299-ACE2", "HT1080/ACE2", "S-Fuse", "Caco-2", "Huh 7")

names(shapes_fold) <- names(shapes)

#shapes webplot
shapes_webplot <- c(23, 21, 22, 22, 22)
names(shapes_webplot) <- c("y", "n", "y (titers)", "y (Titers)", "y (Delta titers)")

# ----------------------------------------------------------- set up required functions
conv_to_logfold <- function(drop) {
  numeric_drop <- as.numeric(gsub("[>x~]", "", drop))
  log2(numeric_drop)
}

check_factor_correct <- function(forest_data) {
  
  data <- as.data.frame(forest_data)
  cols_to_check <- c("standardise_encounters", "vacc_type", "vaccine_manufacturer", "Comparator antigen", "standardised_cell", "SAVE_lab", "standardised_assay")
  
  factor_check <- lapply(cols_to_check, function(x){
    vals <- is.na(unique(as.character(data[,x])))
    
    if("TRUE" %in% vals){
      message(paste0("Factorising created NA in ", x, ". Check the factor levels."))
    } else {
    }
  })
  
  return()
}

factorise <- function(forest_data){

  forest_data$SAVE_lab[is.na(forest_data$SAVE_lab)] <- "n"
  forest_data$SAVE_lab <- as.factor(forest_data$SAVE_lab)
  
  
  # set levels of factors to be ordered by
  forest_data$standardised_cell <- factor(forest_data$standardised_cell, levels = rev(c("Vero", "VeroE6","Vero-TMPRSS2", "VeroE6-TMPRSS2","VeroE6-ACE2/TMPRSS2", "293T-ACE2/VERO", "293T-ACE2/TMPRSS2", "293T-ACE2", "293T",
                                                                                        "H1299-ACE2", "HT1080/ACE2", "S-Fuse", "Caco-2", "Huh 7","MDA-MB-231-hACE2", "NA","Unknown")))
  
  forest_data$standardised_assay <- factor(forest_data$standardised_assay, levels = rev(c("Live-virus", "Pseudovirus", "Surrogate Virus")))
 
  forest_data$standardise_encounters <- as.character(forest_data$standardise_encounters)
  forest_data$standardise_encounters <- factor(forest_data$standardise_encounters, levels = c(rev(c("WT conv", "Alpha conv", "Beta conv", "Gamma conv", "Delta conv")),c("conv","vacc+inf", "inf+vacc", "3.0", "2.0")))
  
  forest_data$vacc_type <- factor(forest_data$vacc_type, levels = rev(c("mRNA", "az+mRNA", "az", "j&j","heterologous", "other", "inf+vacc", "vacc+inf","conv", "WT conv", "Alpha conv", "Beta conv", "Gamma conv", "Delta conv")))
  forest_data$vaccine_manufacturer <- factor(forest_data$vaccine_manufacturer, levels = rev(c("pfizer", "moderna", "mRNA", "az", "j&j", "heterologous", "other", "inf+vacc", "vacc+inf","conv",  "WT conv", "Alpha conv", "Beta conv", "Gamma conv", "Delta conv")))
  forest_data$`Comparator antigen` <- factor(forest_data$`Comparator antigen`, levels = c("WT","Wu-1", "WA1",
                                                                                          "D614G","B.1", "Wu-1?",
                                                                                          "Alpha", "Beta", "Gamma", "Delta"))
  
  check_factor_correct(forest_data)
  return(forest_data)
}


melt_variant_comparator <- function(data, variant) {
  
  #1st get the numerical titre drops, if no variant titers are given
  variant_drops <- data[(!(is.na(data[,paste0("Titre drop ", variant)])) & is.na(data[,paste0("Titers", variant)])),]
  variant_drops <- as.data.frame(variant_drops)
  variant_drops$`Comparator antigen` <- variant
  variant_drops$`TitersHAg` <- NA
  variant_drops$`numerical Titre drop` <- as.numeric(variant_drops[,paste0("Titre drop ", variant)])
  variant_drops$`Titre drop` <- as.numeric(variant_drops[,paste0("Titre drop ", variant)])
  variant_drops$log_fold_change <- log2(variant_drops$`numerical Titre drop`)
  
  variant_df <- data %>% filter(!is.na(data[,paste0("Titers", variant)]))
  variant_df <- as.data.frame(variant_df)
  variant_df$`Comparator antigen` <- variant
  variant_df$TitersHAg <- as.numeric(variant_df[,paste0("Titers", variant)])
  variant_df$`numerical Titre drop` <- variant_df$TitersHAg/variant_df$TitersOmicron
  variant_df$log_fold_change <- log2(variant_df$`numerical Titre drop`)
  variant_df$`Titre drop` <- unlist(lapply(c(1:nrow(variant_df)), function(x){
    if(is.na(variant_df$Uncertainty[x])) {
      round(variant_df$`numerical Titre drop`[x],2)
    } else{
      paste0(variant_df$Uncertainty[x], round(variant_df$`numerical Titre drop`[x],2))
    }
  }
  ))
  
  variant_df$Log2HAg <- log2(variant_df$TitersHAg/10)
  
  
  data <- rbind(data,variant_df, variant_drops)
  data <- unique(data)

  return(data)
}

# rowlabel function
character_rowlabel <- function(forest_data) {
  # get proper rowlabel with monospacing
  forest_data$time[is.na(forest_data$time)] <- ""
  study_length<- max(nchar(forest_data$Study)) + 1
  sera_length<- max(nchar(forest_data$Standardised_sera_names))
  timing_length <- max(nchar(forest_data$time)) + 1
  
  
  total_length <- study_length + sera_length + timing_length
  #spaces_to_timing <- rowlabel_length - timing_length
  
  forest_data$rowlabel_old <- forest_data$rowlabel
  forest_data$rowlabel <-  unlist(lapply(1:nrow(forest_data), function(x) {
    study <- forest_data$`Study`[x]
    serum <- forest_data$Standardised_sera_names[x]
    time <- forest_data$time[x]
    
    
    spaces_to_time <- study_length - nchar(study)
    spaces_to_end <- sera_length -nchar(serum)
    if(time == "") {
      spaces_to_study <- timing_length
      
      paste0(study, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
     
    } else {
      spaces_to_study <- timing_length - nchar(time)
      
      paste0(study, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
      
    }
    
    
  }))
  
  return(forest_data)
}

character_rowlabel_standard_time <- function(forest_data) {
  # get proper rowlabel with monospacing
  forest_data$begin_d <- round(forest_data$begin_d/30,2)
  forest_data$end_d <- round(forest_data$end_d/30,2)
  
  forest_data$time <- unlist(lapply(1:nrow(forest_data), function(x) {
    if(is.na(forest_data$time[x])) {
      return("")
    }else if(!(is.na(forest_data$end_d[x]))) {
      temp <- paste0(forest_data$begin_d[x], "-", forest_data$end_d[x], "m ", str_split_fixed(forest_data$full_time[x], " ", 2)[1,2])
    } else {
      temp <- paste0(forest_data$begin_d[x], "m ", str_split_fixed(forest_data$full_time[x], " ", 2)[1,2])
    }
    if(forest_data$shorter[x]) {
      temp <- paste0("<",temp)
    } else if(forest_data$longer[x]) {
      temp <- paste0(">", temp)
    }
    return(temp)
  }))
  
 
  study_length<- max(nchar(forest_data$Study)) + 1
  sera_length<- max(nchar(forest_data$Standardised_sera_names))
  timing_length <- max(nchar(forest_data$time)) + 1
  
  
  total_length <- study_length + sera_length + timing_length
  #spaces_to_timing <- rowlabel_length - timing_length
  
  forest_data$rowlabel_old <- forest_data$rowlabel
  forest_data$rowlabel <-  unlist(lapply(1:nrow(forest_data), function(x) {
    study <- forest_data$`Study`[x]
    serum <- forest_data$Standardised_sera_names[x]
    time <- forest_data$time[x]
    
    
    spaces_to_time <- study_length - nchar(study)
    spaces_to_end <- sera_length -nchar(serum)
    if(time == "") {
      spaces_to_study <- timing_length
      
      paste0(study, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
      #  paste0(label, paste0(rep(" ", rowlabel_length - nchar(label)), collapse = ""))
    } else {
      spaces_to_study <- timing_length - nchar(time)
      
      paste0(study, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
      # paste0(label, paste0(rep(" ", spaces_to_timing - nchar(label)), collapse = ""), time)
    }
    
    
  }))
  
  return(forest_data)
}


# assay labelling function
add_assay_cell_type_label <- function(data) {
  data$label_assay_cell_type <- paste0(as.character(data$standardised_assay), " (", as.character(data$standardised_cell), ")")
  
  return(data)
}

add_assay_type_label <- function(data) {
  data$label_assay_type <- as.character(data$standardised_assay)
  
  return(data)
}

standardise_time_to_mean_days <- function(data) {
  time_df <- data.frame("full_time" = data$time)
  time_df <- time_df %>%
    mutate("small_time" = gsub("<|>", "", str_split_fixed(full_time, " ", 2)[,1]),
           "shorter" = grepl("<", full_time),
           "longer" = grepl(">", full_time),
           "unit" = str_sub(small_time, start = -1),
           "timeframe" = str_sub(small_time, end = -2))
  
  time_df[time_df$small_time == "pre-boost", c("timeframe", "unit")] <- ""
  time_df$begin <- unlist(lapply(time_df$timeframe, function(x) strsplit(x, "-")[[1]][1]))
  time_df$end <- unlist(lapply(time_df$timeframe, function(x) strsplit(x, "-")[[1]][2]))
  
  unit_factor <- c("d" =1, "m" = 30, "w" = 7)
  
  time_df <- time_df %>% 
    mutate(time_factor = unit_factor[unit],
           begin_d = as.numeric(begin)*time_factor,
           end_d = as.numeric(end)*time_factor)
  
  time_df$mean_days <- rowMeans(time_df[,c("begin_d", "end_d")], na.rm = TRUE)
  
  data <- cbind(data, time_df)
    
}

# rowlabel function
character_rowlabel_assay_cell <- function(forest_data) {
  
  forest_data <- add_assay_cell_type_label(forest_data)
  # get proper rowlabel with monospacing
  forest_data$time[is.na(forest_data$time)] <- ""
  study_length<- max(nchar(forest_data$Study)) + 2
  assay_length <- max(nchar(forest_data$label_assay_cell_type)) + 2
  sera_length<- max(nchar(forest_data$Standardised_sera_names))
  timing_length <- max(nchar(forest_data$time)) + 2
  
  
  total_length <- study_length + sera_length + timing_length + assay_length
  #spaces_to_timing <- rowlabel_length - timing_length
  
  forest_data$rowlabel_old <- forest_data$rowlabel
  forest_data$rowlabel <-  unlist(lapply(1:nrow(forest_data), function(x) {
    study <- forest_data$`Study`[x]
    assay <- forest_data$label_assay_cell_type[x]
    serum <- forest_data$Standardised_sera_names[x]
    time <- forest_data$time[x]
    
    spaces_to_assay <- study_length - nchar(study)
    spaces_to_time <- assay_length - nchar(assay)
    spaces_to_end <- sera_length -nchar(serum)
   
    if(time == "") {
      spaces_to_study <- timing_length
      paste0(study, paste0(rep(" ", spaces_to_assay), collapse = ""), assay, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
    } else {
      spaces_to_study <- timing_length - nchar(time)
      paste0(study, paste0(rep(" ", spaces_to_assay), collapse = ""), assay, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
    }
    
    
  }))
  
  return(forest_data)
}

# rowlabel function
character_rowlabel_assay <- function(forest_data) {
  
  forest_data <- add_assay_type_label(forest_data)
  # get proper rowlabel with monospacing
  forest_data$time[is.na(forest_data$time)] <- ""
  study_length<- max(nchar(forest_data$Study)) + 2
  assay_length <- max(nchar(forest_data$label_assay_type)) + 2
  sera_length<- max(nchar(forest_data$Standardised_sera_names))
  timing_length <- max(nchar(forest_data$time)) + 2
  
  
  total_length <- study_length + sera_length + timing_length + assay_length
  #spaces_to_timing <- rowlabel_length - timing_length
  
  forest_data$rowlabel_old <- forest_data$rowlabel
  forest_data$rowlabel <-  unlist(lapply(1:nrow(forest_data), function(x) {
    study <- forest_data$`Study`[x]
    assay <- forest_data$label_assay_type[x]
    serum <- forest_data$Standardised_sera_names[x]
    time <- forest_data$time[x]
    
    spaces_to_assay <- study_length - nchar(study)
    spaces_to_time <- assay_length - nchar(assay)
    spaces_to_end <- sera_length -nchar(serum)
    
    if(time == "") {
      spaces_to_study <- timing_length
      paste0(study, paste0(rep(" ", spaces_to_assay), collapse = ""), assay, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
    } else {
      spaces_to_study <- timing_length - nchar(time)
      paste0(study, paste0(rep(" ", spaces_to_assay), collapse = ""), assay, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
    }
    
    
  }))
  
  return(forest_data)
}


character_rowlabel_study_assay <- function(forest_data) {
  abbreviations_assay <- c("PV", "LV", "SV")
  names(abbreviations_assay) <- c("Pseudovirus", "Live-virus", "Surrogate Virus")
  
  
  forest_data$Study <- paste0(forest_data$Study, " (", abbreviations_assay[forest_data$standardised_assay], ")")
  # get proper rowlabel with monospacing
  forest_data$time[is.na(forest_data$time)] <- ""
  study_length<- max(nchar(forest_data$Study)) + 1
  sera_length<- max(nchar(forest_data$Standardised_sera_names))
  timing_length <- max(nchar(forest_data$time)) + 2
  
  
  total_length <- study_length + sera_length + timing_length
  #spaces_to_timing <- rowlabel_length - timing_length
  
  forest_data$rowlabel_old <- forest_data$rowlabel
  forest_data$rowlabel <-  unlist(lapply(1:nrow(forest_data), function(x) {
    study <- forest_data$`Study`[x]
    serum <- forest_data$Standardised_sera_names[x]
    time <- forest_data$time[x]
    
    
    spaces_to_time <- study_length - nchar(study)
    spaces_to_end <- sera_length -nchar(serum)
    if(time == "") {
      spaces_to_study <- timing_length
      
      paste0(study, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
      #  paste0(label, paste0(rep(" ", rowlabel_length - nchar(label)), collapse = ""))
    } else {
      spaces_to_study <- timing_length - nchar(time)
      
      paste0(study, paste0(rep(" ", spaces_to_time), collapse = ""), time, paste0(rep(" ", spaces_to_study), collapse = ""), serum, paste0(rep(" ", spaces_to_end), collapse = ""))
      # paste0(label, paste0(rep(" ", spaces_to_timing - nchar(label)), collapse = ""), time)
    }
    
    
  }))
  
  return(forest_data)
}


reorder_data <- function(data, rev = TRUE){
  if(rev) {
    data$Row_long <- factor(data$Row_long, levels = rev(data$Row_long))
  } else {
    data$Row_long <- factor(data$Row_long, levels = data$Row_long)
  }
  
  data <- data[order(data$Row_long),]
}

get_hdiv_subgroup <- function(data, hline_by, super_pos) {
  hline_rows <- get_hdiv_positions(data[c(1:super_pos[2]), ], hline_by)
  
  if(length(super_pos) > 2) {
    for(p in 2:(length(super_pos)-1)) {
      hline_rows <- c(hline_rows, super_pos[p] + get_hdiv_positions(data[c((super_pos[p]+1): super_pos[p+1]), ], hline_by))
    }
    hline_rows <- c(hline_rows, super_pos[p+1] + get_hdiv_positions(data[c((super_pos[p+1]+1): nrow(data)), ], hline_by))
  }
  
  combined <- sort(unique(c(hline_rows, super_pos)))
 
  return(combined) 
}

get_hdiv_positions <- function(data, hline_by) {
  data <- as.data.frame(data)
 
  if(length(hline_by)>1) {
    hline_positions1 <- match(unique(data[,hline_by[1]]), data[,hline_by[1]])
    
    hline_positions <- unlist(lapply(hline_by[2:length(hline_by)], function(x) get_hdiv_subgroup(data, x, hline_positions1)))
    
  } else {
    hline_positions <- match(unique(data[,hline_by]), data[,hline_by])
  }
   # hline_positions <- unlist(lapply(hline_by, function(x) {
  #    match(unique(data[,x]), data[,x])
  #  }))
  
 # hline_positions <- match(unique(data[,hline_by]), data[,hline_by])
  
  return(hline_positions)
}

add_hdiv_to_data <- function(data, hline_by) {
  
  # # test adding empty rows
  hlines <- get_hdiv_positions(data, hline_by[1])
  
  rows <- c()
  if(length(hlines) >1) {
    for(pos in 2:length(hlines)) {
      rows <- c(rows, hlines[pos-1]:hlines[pos])
    }
  } 
  rows <- c(rows, last(hlines):nrow(data))
  
  ## adding empty row at hline
  test <- data[rows,]
  
  new_hdivs <- get_hdiv_positions(test, hline_by)
  if(length(new_hdivs) > 1) {
    test[new_hdivs[2:length(new_hdivs)],] <- NA
    hlines <- new_hdivs[2:length(new_hdivs)]
  } else {
    hlines <- NA
  }
  
  test$Row_long <- c(1:nrow(test))
  test$Row_long <- factor(test$Row_long, levels = test$Row_long)
  
  return(list("data" = test, "hlines" = hlines))
}

hdiv_to_data_multiple_groupings <- function(data, hline_by) {
  res <- add_hdiv_to_data(data, hline_by[1])
  data <- res$data
  hline_positions <- res$hlines
  
  
  temp_res <- add_hdiv_to_data(data[c(1:hline_positions[1]),], hline_by[2])
  data_sub <- temp_res$data
  lines <- temp_res$hlines
  
  # now do it for the subgroupings and cbind it
    if(length(hline_positions)>1) {
      for(p in 1:(length(hline_positions)-1)) {
        temp_res <- add_hdiv_to_data(data[c((hline_positions[p]+1): hline_positions[p+1]),], hline_by[2])
        data_sub <- rbind(data_sub, temp_res$data)
        lines <- c(lines, temp_res$hlines)
      }
      temp_res <- add_hdiv_to_data(data[c((hline_positions[p+1]+1): nrow(data)),], hline_by[2])
      data_sub <- rbind(data_sub, temp_res$data)
      lines <- c(lines, temp_res$hlines)
    } else {
      temp_res <- add_hdiv_to_data(data[c((hline_positions[1]+1): nrow(data)),], hline_by[2])
      data_sub <- rbind(data_sub, temp_res$data)
      lines <- c(lines, temp_res$hlines)
    } 
  
  
  lines <- as.integer(which(rowSums(is.na(data_sub))==ncol(data_sub)-1))

  data_sub$Row_long <- c(1:nrow(data_sub))
  data_sub$Row_long <- factor(data_sub$Row_long, levels = data_sub$Row_long)
  return(list("data" = data_sub, "hlines" = lines))
  
}

pretty_plot_names <- function(data){
  
  # set standard encounter names
  pretty_names_encounters <- c("2.0" =  "2x Vax", "3.0" = "3x Vax","inf+vacc"= "Inf + 2x Vax", "vacc+inf" = "2x Vax + Inf", 
                               "WT conv" =  "WT conv" ,"Alpha conv" =  "Alpha conv", 
                               "Beta conv" = "Beta conv","Gamma conv"= "Gamma conv", "Delta conv"= "Delta conv", "conv" = "Conv")
  
  pretty_names_vax_manuf <- c(rev(c("inf+vacc"= "Inf + 2x Vax", "vacc+inf" = "2x Vax + Inf","WT conv" =  "WT conv" ,"Alpha conv" =  "Alpha conv", 
                                    "Beta conv" = "Beta conv","Gamma conv"= "Gamma conv", "Delta conv"= "Delta conv", "conv" = "Conv" )), 
                              c("heterologous" = "Het", "other" = "Other", "j&j" = "J&J", "az" = "AZ", "mRNA" = "mRNA", "moderna" = "Moderna", "pfizer" = "Pfizer"))
  
  pretty_names_vax_type <- c(rev(c("inf+vacc"= "Inf + 2x Vax", "vacc+inf" = "2x Vax + Inf","WT conv" =  "WT conv" ,"Alpha conv" =  "Alpha conv", 
                                    "Beta conv" = "Beta conv","Gamma conv"= "Gamma conv", "Delta conv"= "Delta conv", "conv" = "Conv" )), 
                              c("heterologous" = "Het", "other" = "Other", "j&j" = "J&J", "az" = "AZ","az+mRNA"="AZ + mRNA", "mRNA" = "mRNA"))
  
  pretty_names_assay <- c("Live-virus" = "LV", "Pseudovirus" = "PV", "Surrogate Virus" = "SV")
  
 
  data$standardise_encounters <- pretty_names_encounters[as.character(data$standardise_encounters)]
  data$standardise_encounters <- factor(data$standardise_encounters, level = rev(pretty_names_encounters))
 
  data$vacc_type <- pretty_names_vax_type[as.character(data$vacc_type)]
  data$vacc_type <- factor(data$vacc_type, level = pretty_names_vax_type)
  
  
  data$vaccine_manufacturer <- pretty_names_vax_manuf[as.character(data$vaccine_manufacturer)]
  data$vaccine_manufacturer <- factor(data$vaccine_manufacturer, level = pretty_names_vax_manuf)
  
  data$standardised_assay <- pretty_names_assay[as.character(data$standardised_assay)]
  data$standardised_assay <- factor(data$standardised_assay, level = rev(pretty_names_assay))

  check_factor_correct(data)
  
  return(data)
  
}

add_horizontal_labelling <- function(base_plot, data, hline_by, hline_positions) {
  
  y_adjust <- 0.5
  
  breaks <- c(hline_positions, nrow(data) +1)
  
  bs <- c()
  for(x in 1: length(breaks)) {
    if(x > 1) {
      if(breaks[x] == breaks[x-1]+1) {
        bs <- c(bs, breaks[x]-y_adjust)
      }
    } 
  }
  
  breaks <- bs
  
  base_plot + coord_cartesian(clip = "off") -> base_plot
  
  for(b in breaks) {
    base_plot <- base_plot + 
      annotation_custom(grob = linesGrob(gp = gpar(col = "grey60")), 
                        xmin =-50, 
                        xmax = 30,
                        ymin = b,
                        ymax = b)
  }
  
  
  return(base_plot)
}


vertical_labelling_plot <- function(base_plot, data, hline_by, hline_positions, axis_text_size) {
  
  l_hline_by <- length(hline_by)
  y_adjust <- 1
  
  text_size <- axis_text_size*2
  breaks <- c(0, hline_positions, nrow(data) +1)
  
  if(l_hline_by > 1) {
    outer_by <- hline_by[1]
    
    hline_by <- hline_by[2]
    if(hline_by != "standardised_cell") {
      text_size <- text_size*0.7
    } else {
     
      text_size <- text_size*0.5
    }
    
    y_adjust <- 1
    
    bs <- c()
    
    outer <- c()
    for(x in 1: length(breaks)) {
      if(x > 1) {
        if(breaks[x] != breaks[x-1]+1) {
          bs <- c(bs, breaks[x])
        } else {
          outer <- c(outer, breaks[x])
        }
      } else {
        bs <- c(bs, breaks[x])
      }
    }
    
    outer <- c(outer, nrow(data)+1)
    breaks <- bs
  } else {
    outer_by <- ""
  }
  
  
  labels <- data[breaks[2:length(breaks)]-1,hline_by] %>% pull(hline_by)
  if(outer_by %in% c("vacc_type", "vaccine_manufacturer")) {
    labels[grepl("conv|nf", labels)] <- ""
  }
  
  base_plot + scale_x_continuous(breaks = NULL, limits = c(0,2.5)) + 
    theme_void() -> test_plot
  
  test_plot + 
    coord_cartesian(clip = "off") +
    annotate(geom = "text", label = labels, 
             x = rep(2.4, length(labels)), 
             y = breaks[2:length(breaks)]-y_adjust, 
             size = text_size/figure_factor_increase,
             hjust = 1) +
    theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) -> test_plot
  
  if(l_hline_by >1) {
   
    outer_label <- levels(data %>% pull(outer_by))
   
    outer_label <- outer_label[outer_label %in% as.character(unique(data %>% pull(outer_by)))]
      
  test_plot +
     annotate(geom = "text", label = outer_label,
             x = rep(0, length(outer_label)),
            y = outer-2,
           size = text_size,
          hjust = 0) -> test_plot
  }
  
  for(pos in 1:(length(breaks)-1)) {
    test_plot <- test_plot + geom_segment(x = 2.5, xend = 2.5, y = breaks[pos]+1, yend = breaks[pos+1]-1, color = "grey80")
    
  }
  
  return(test_plot)
}

add_labelling_plot <- function(full_plot, label_plot) {
  
  plot <- full_plot + 
    theme(plot.margin = unit(c(0.5, 0, 0.5,0), "lines"))
  
  plot <- ggplotGrob(plot)
  g <- plot
  
  #grob out of label plot
  label_grob <- ggplotGrob(label_plot)
  g1 <- gtable_filter(label_grob, "panel")
  index <- subset(g$layout, name == "panel")
  g <- gtable_add_cols(g, unit(1.5, "strwidth", "line # 1") + unit(1.5, "cm"), pos = 0)
  
  # do one plot with xlim small and write only on the panel
  # then add this panel as column to the left
  gt <- gtable_add_grob(g, g1, t = index$t, l=1,  b=index$b, r=1)
  
  return(gt)
}

calc_mean_per_grouping <- function(data, hline_by) {
  
  temp_data <- data
  temp_data$`Comparator antigen`[!(temp_data$`Comparator antigen` %in% c("Alpha", "Beta", "Delta", "Gamma"))] <- "WT"
  
  # means <- data %>%
  #   group_by(across(all_of(c("Comparator antigen", hline_by)))) %>%
  #   summarise(mean_fold_drop = Rmisc::CI(log_fold_change)[["mean"]], 
  #             mean_fold_drop_lower = Rmisc::CI(log_fold_change)[["lower"]], 
  #             mean_fold_drop_upper = Rmisc::CI(log_fold_change)[["upper"]],
  #             gmt = mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean,
  #             gmt_lower = mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_lower,
  #             gmt_upper = mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_upper,
  #             gmt_o = mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean,
  #             gmt_lower_o = mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_lower,
  #             gmt_upper_o = mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_upper
  #          ) 
  # 
  means <- temp_data %>%
    group_by(across(all_of(c("Comparator antigen", hline_by)))) %>%
    summarise(log_fold_change= mean(log_fold_change, na.rm = TRUE),
              Log2HAg = mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean,
              Log2Omi = mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean,
              rowlabel = "Mean",
              Webplotdigitizer = "n"
    ) %>% ungroup()
  
  data <- bind_rows(data, means)
  return(data)
  
}

do_forest_titer_drop <- function(forest_data, hline_by = "", row_label, show_mean, axis_text_size) {
  
  comp_antigen <- as.character(unique(forest_data$`Comparator antigen`))
  if("B.1" %in% comp_antigen | "D614G" %in% comp_antigen) {
    comp_antigen <- c(comp_antigen, "WT")
  }
  if(show_mean) {
    forest_data <- calc_mean_per_grouping(forest_data, hline_by)
    
    forest_data <- forest_data %>% filter(`Comparator antigen` %in% comp_antigen) 
    
    forest_data <- forest_data %>%
      arrange(desc(log_fold_change)) %>%
      arrange(across(all_of(c(hline_by))))
    
    forest_data$Row_long <- 1:nrow(forest_data)
    
    forest_data <- reorder_data(forest_data, rev = FALSE)
  }
  
  if(length(hline_by) >0 & hline_by != "") {
    
    if(length(hline_by) == 1) {
      res <- add_hdiv_to_data(forest_data, hline_by)
      forest_data <- res$data
      hline_positions <- res$hlines
      
      
    } else {
      hline_size <- hline_size/2
      res <- hdiv_to_data_multiple_groupings(forest_data, hline_by)
      forest_data <- res$data
      hline_positions <- res$hlines
    
    }
   
    
  } else {
    hline_positions <- NA
  }
  
  xmin <- 8
  xmax <- -2
  if("Delta" %in% unique(as.character(forest_data$`Comparator antigen`))) {
    xmin <- 8.5
    xmax <- -4.5
  }
  if("Beta" %in% unique(as.character(forest_data$`Comparator antigen`))) {
    xmax <- -6
  }
 # print(hline_positions)
  forest_data$shape_col <- as.factor(forest_data$Webplotdigitizer)
  save_data <- forest_data %>% filter(SAVE_lab == "y")
  
  # colors
  fill_colors <- c("lightblue", colors[unique(as.character(forest_data$`Comparator antigen`))])
  names(fill_colors) <- c("y", rep("n", length(fill_colors)-1))
  
  if(show_mean) {
    mean_data <- forest_data %>% filter(rowlabel == "Mean")
  }
  if(!row_label) {
    forest_data$rowlabel[!(is.na(forest_data$rowlabel))] <- ""
  }
  
  forest_data %>% ggplot(
    aes(
      x = log_fold_change,
      y = Row_long,
      color = `Comparator antigen`,
      fill = SAVE_lab,
      shape = shape_col
    )
  ) + scale_y_discrete(
    labels = forest_data$rowlabel[!is.na(forest_data$rowlabel)],
    breaks = forest_data$Row_long[!is.na(forest_data$rowlabel)]
  ) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(
        hjust = 0,
        family = "mono",
        size = axis_text_size
      ),
      axis.text.x = element_text(
        family = "mono",
        angle = 90,
        size = 8
      )
    ) +labs(
      x = "Fold drop",
      y = ""
    ) +
    scale_x_reverse(
      labels = function(x) 2^x,
      breaks = floor(xmin):floor(xmax)
    ) +
    coord_cartesian(
      xlim = c(xmin, xmax)
    ) + geom_point(size = point_size, color = "white") -> base_plot
    
    
    base_plot + 
    geom_segment(
      xend =- (forest_data$log_fold_change + forest_data$arrow_length),
      y = c(1:nrow(forest_data)),
      yend = c(1:nrow(forest_data)),
      arrow = arrow(length = unit(arrow_head_size, "inches"), type = "closed"),
      color = "grey58",
      size = arrow_line_size
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "solid"
    ) +
    geom_hline(yintercept = hline_positions, size =hline_size, color = "white") + 
    geom_point(size = point_size) + 
    scale_color_manual(values = colors)+
    geom_point(data = save_data, aes(
      x = log_fold_change,
      y = Row_long,
      color = `Comparator antigen`,
      fill = SAVE_lab,
    ), size = point_size) + #, shape = 21) +
    scale_fill_manual(values = fill_colors) +
    scale_shape_manual(values = shapes_webplot) -> gp
    
    if(!row_label) {
      gp <- gp + scale_y_discrete(
        breaks = NULL
      )
    }
    
    if(show_mean) {
      gp +
        geom_point(data = mean_data, 
                   aes(x = log_fold_change, 
                       y= Row_long),
                   shape = 21, stroke = 2,
                       size =6, color = "red", fill = fill_colors[unique(as.character(mean_data$`Comparator antigen`))]) -> gp
    }
    
    # ad horizontal labelling
    if(length(hline_by) > 1) {
      gp <- add_horizontal_labelling(gp, forest_data, hline_by[1], hline_positions)
    }
  
  label_plot <- vertical_labelling_plot(base_plot = base_plot, data = forest_data, hline_by = hline_by, hline_positions = hline_positions, axis_text_size = axis_text_size_p)
  gp <- add_labelling_plot(gp, label_plot)
  
  return(gp)
}

do_ve_plot <- function(data, hline_by = "", row_label, show_mean, axis_text_size, order_by_omicron = FALSE) {
  
  comp_antigen <- as.character(unique(data$`Comparator antigen`))
  if("B.1" %in% comp_antigen | "D614G" %in% comp_antigen) {
    comp_antigen <- c(comp_antigen, "WT")
  }
  if(show_mean) {
    data <- calc_mean_per_grouping(data, hline_by)
    data <- data %>% filter(`Comparator antigen` %in% comp_antigen) 
    
    if(order_by_omicron){
      data <- data %>%
        arrange(Log2Omi) %>%
        arrange(across(all_of(c(hline_by))))
    } else {
      data <- data %>%
        arrange(Log2HAg) %>%
        arrange(across(all_of(c(hline_by))))
    }
    
    
    data$Row_long <- 1:nrow(data)
    
    data <- reorder_data(data, rev = FALSE)
    
  }
  
  if(length(hline_by) >0 & hline_by != "") {
    
    if(length(hline_by) == 1) {
      res <- add_hdiv_to_data(data, hline_by)
      data <- res$data
      hline_positions <- res$hlines
    } else {
      res <- hdiv_to_data_multiple_groupings(data, hline_by)
      data <- res$data
      hline_positions <- res$hlines
    }
    
    
  } else {
    hline_positions <- NA
  }
  
  #shape
  data$shape_col <- as.factor(data$Webplotdigitizer)
  
  if(show_mean) {
    mean_data <- data %>% filter(rowlabel == "Mean")
  }
  
  if(!row_label) {
    data$rowlabel[!(is.na(data$rowlabel))] <- ""
  }
  
  base_plot <- data %>% ggplot(
    aes(
      x = Log2Omi,
      y = Row_long,
      shape = shape_col
    )
  ) + 
    geom_point(size = point_size, color = "white") +
    scale_x_continuous(
      labels = x_labels,
      breaks = seq(from = 0, to = 15, by = 2),
      limits = c(-3.5, 15.5),
      expand = c(0,0),
      
      sec.axis = dup_axis(breaks = cps$Log2Titer, labels = c("a", "d", "f  g", "b", "c", "e", ""), name = "")
      #    breaks = 0:6
    ) +
    scale_y_discrete(
      labels = data$rowlabel[!is.na(data$rowlabel)],
      breaks = data$Row_long[!is.na(data$rowlabel)]
    ) +
    theme(
      legend.position = "none",# "right",
      axis.text.x.bottom = element_text(angle = 90),
      axis.text.y = element_text(
        hjust = 0,
        family = "mono",
        size = axis_text_size
      ),
      axis.text.x = element_text(
        family = "mono"
      )
    )
    
  base_plot +
    geom_vline(
      xintercept = 0,
      linetype = "solid"
    ) +
    geom_vline(
      xintercept = cps$Log2Titer,
      linetype = "dashed",
      color = "grey56",
      alpha  = 0.8
    ) +
    geom_segment(
      xend = data$Log2HAg,
      y = c(1:nrow(data)),
      yend = c(1:nrow(data)),
      color = "grey44",
      size = arrow_line_size
    ) +
    geom_segment(
      xend = data$Log2Omi - 2*data$arrow_length,
      y = c(1:nrow(data)),
      yend = c(1:nrow(data)),
      arrow = arrow(length = unit(arrow_head_size, "inches"), type = "closed"),
      color = "grey58",
      size = arrow_line_size
    ) +
    geom_hline(yintercept = hline_positions, size =hline_size, color = "white") + 
    geom_point(aes(fill = SAVE_lab), size = point_size, color = colors["Omicron"])+#, shape = 21) + 
    geom_point(data = data, aes(
      x = Log2HAg,
      y = Row_long,
      color = `Comparator antigen`
    ), shape = 16) +
    scale_color_manual(values = colors)+
    scale_fill_manual(values = fill_colors) + 
    scale_shape_manual(values = shapes_webplot) +
    labs(
      x = "Titer",
      y = ""
    )  ->plot
  
  if(show_mean) {
   
    plot +
      geom_point(data = mean_data, 
                 aes(x = Log2HAg,
                     y = Row_long,
                     color = `Comparator antigen`),
                 shape = 21, stroke = 2,
                 size =5) +
      geom_point(data = mean_data, 
                 aes(x = Log2Omi,
                     y = Row_long),
                 shape = 21, stroke = 2, color = colors["Omicron"],
                 size =5) -> plot
  }
  
  if(!row_label) {
    plot <- plot + scale_y_discrete(breaks = NULL)
  }
  
  
  if(length(hline_by) > 1) {
    plot <- add_horizontal_labelling(plot, data, hline_by[1], hline_positions)
  }
  
  label_plot <- vertical_labelling_plot(base_plot = base_plot, data = data, hline_by = hline_by, hline_positions = hline_positions, axis_text_size = axis_text_size_p)
  plot <- add_labelling_plot(plot, label_plot)
  
  
  return(plot)
}


split_up_forest_plot_drops <- function(data, hline_by = "", row_label = TRUE, show_means = FALSE, axis_text_size = axis_text_size_p) {
  show_means_var <- show_means
  if(length(hline_by) > 1) {
    show_means_var <- FALSE
  }
  
  if("standardised_cell" %in% hline_by) {
    show_means <- FALSE
  }
  
  wt <- data %>% filter(`Comparator antigen` %in% c("WT","Wu-1", "WA1","D614G","B.1", "Wu-1?")) %>% do_forest_titer_drop(hline_by, row_label, show_means, axis_text_size)
  beta <- data %>% filter(`Comparator antigen` %in% c("Beta")) %>% do_forest_titer_drop(hline_by, row_label, show_means_var, axis_text_size)
  alpha <- data %>% filter(`Comparator antigen` %in% c("Alpha")) %>% do_forest_titer_drop(hline_by, row_label, show_means_var, axis_text_size)
  delta <- data %>% filter(`Comparator antigen` %in% c("Delta")) %>% do_forest_titer_drop(hline_by, row_label, show_means_var, axis_text_size)
  
  return(list("WT" = wt, "Beta" = beta, "Alpha" = alpha, "Delta" = delta))
  
} 

split_up_forest_plot_titer <- function(data, hline_by = "", row_label = TRUE, show_means = FALSE, axis_text_size = axis_text_size_p, order_by_omi = FALSE) {
  show_means_var <- show_means
  if(length(hline_by) > 1) {
    show_means_var <- FALSE
  }
  if("standardised_cell" %in% hline_by) {
    show_means <- FALSE
  }
  
  wt <- data %>% filter(`Comparator antigen` %in% c("WT","Wu-1", "WA1","D614G","B.1", "Wu-1?")) %>% do_ve_plot(hline_by, row_label, show_means, axis_text_size, order_by_omicron = order_by_omi)
  beta <- data %>% filter(`Comparator antigen` %in% c("Beta")) %>% do_ve_plot(hline_by, row_label, show_means_var, axis_text_size, order_by_omicron = order_by_omi)
  alpha <- data %>% filter(`Comparator antigen` %in% c("Alpha")) %>% do_ve_plot(hline_by, row_label, show_means_var,axis_text_size, order_by_omicron = order_by_omi)
  delta <- data %>% filter(`Comparator antigen` %in% c("Delta")) %>% do_ve_plot(hline_by, row_label, show_means_var,axis_text_size, order_by_omicron = order_by_omi)
  
  return(list("WT" = wt, "Beta" = beta, "Alpha" = alpha, "Delta" = delta))
  
} 

#old width and height:  width = 12, height = 18
save_in_scaled_format <- function(data, which_plot, width = 8, height = 12, to_save= "plots", hline_by = "", single_plots = TRUE, axis_text_size=axis_text_size_p, order_by_omicron = FALSE) {
  nrow_wt <- nrow(data %>% filter(`Comparator antigen` %in% c("WT","Wu-1", "WA1","D614G","B.1", "Wu-1?")))
  nrow_alpha <- nrow(data %>% filter(`Comparator antigen` %in% c("Alpha")))
  nrow_beta <- nrow(data %>% filter(`Comparator antigen` %in% c("Beta")))
  nrow_delta <-  nrow(data %>% filter(`Comparator antigen` %in% c("Delta")))
  
  base_t <- 0
  spacer <- 1
  
  if(which_plot == "titer_drop") {
    plots <- split_up_forest_plot_drops(data, hline_by, row_label = TRUE, show_means = TRUE, axis_text_size)
    if(hline_by[1] == "standardise_encounters") {
      plots_no_lab <- split_up_forest_plot_drops(data, hline_by, row_label = FALSE, show_means = TRUE, axis_text_size)
    }
    
    
    
    to_save <- paste0(path_to_save, to_save, "_fold_drops")
    height_a <- height #  + 4*spacer
    height_bcd <- height # -2 *spacer
    
    height_a <- height
    height_bcd <- height*(nrow_alpha+nrow_beta+nrow_delta)/nrow_wt #-2 *spacer
    
  } else {
    plots <- split_up_forest_plot_titer(data, hline_by, row_label = TRUE, show_means = TRUE, axis_text_size, order_by_omi = order_by_omicron)
    if(hline_by[1] == "standardise_encounters") {
      plots_no_lab <- split_up_forest_plot_titer(data, hline_by, row_label = FALSE, show_means = TRUE, axis_text_size, order_by_omi = order_by_omicron)
    }
    
    to_save <- paste0(path_to_save, to_save, "_titers")
    height_a <- height + 3
    height_bcd <- height + 3
    
    
  }
  
  layout_wt <- c(area(t = base_t, l = 1, b = base_t+nrow_wt, r = 2))
  layout_alpha <- c(area(t = base_t, l = 1, b = base_t+nrow_alpha, r = 2))
  layout_beta <- c(area(t = base_t, l = 1, b = base_t+nrow_beta, r = 2))
  layout_delta <- c(area(t = base_t, l = 1, b = base_t+nrow_delta, r = 2))
  layout_bcd <- c(
    area(t = base_t, l = 1, b = base_t+nrow_alpha +spacer, r = 2),
    area(t = 2*spacer+(base_t+nrow_alpha), l = 1, b = 2*spacer+base_t+nrow_alpha+nrow_beta, r = 2),
    area(t = 3*spacer+(base_t+nrow_beta + nrow_alpha), l = 1, b = 3*spacer+base_t+nrow_alpha+nrow_beta+ nrow_delta, r = 2)
  )
  
  
#  bcd <-plots$Alpha + plots$Beta + plots$Delta + plot_layout(design = layout_bcd) + 
 #   plot_annotation(tag_levels = "A") &
#    theme(plot.tag = element_text(size = 16, face = "bold"))

    ggsave(paste0(to_save, "_WT.",fileext), plots$WT, width = width, height = height_a, dpi = 300)
    ggsave(paste0(to_save, "_Alpha.",fileext), plots$Alpha, width = width, height = height_a*((nrow_alpha*1.2)/nrow_wt), dpi = 300)
    ggsave(paste0(to_save, "_Beta.",fileext), plots$Beta, width = width, height = height_a*(nrow_beta/nrow_wt), dpi = 300)
    ggsave(paste0(to_save, "_Delta.",fileext), plots$Delta, width = width, height = height_a*(nrow_delta/nrow_wt), dpi = 300)
    
    if(hline_by[1] == "standardise_encounters") {
    ggsave(paste0(to_save, "_no_lab_WT.",fileext), plots_no_lab$WT, width = width, height = height_a*0.85, dpi = 300)
  #  ggsave(paste0(to_save, "_no_lab_Alpha.",fileext), plots_no_lab$Alpha, width = width, height = height_a*((nrow_alpha*1.2)/nrow_wt)*0.85, dpi = 300)
  #  ggsave(paste0(to_save, "_no_lab_Beta.",fileext), plots_no_lab$Beta, width = width, height = height_a*(nrow_beta/nrow_wt)*0.85, dpi = 300)
  #  ggsave(paste0(to_save, "_no_lab_Delta.",fileext), plots_no_lab$Delta, width = width, height = height_a*(nrow_delta/nrow_wt)*0.85, dpi = 300)
    }
  
}


#old width and height:  width = 12, height = 18
save_in_scaled_format_encounter <- function(data, which_plot, width = 8, height = 12, to_save= "plots", hline_by = "", single_plots = TRUE, axis_text_size=axis_text_size_p) {
  nrow_wt <- nrow(data %>% filter(`Comparator antigen` %in% c("WT","Wu-1", "WA1","D614G","B.1", "Wu-1?")))
  nrow_alpha <- nrow(data %>% filter(`Comparator antigen` %in% c("Alpha")))
  nrow_beta <- nrow(data %>% filter(`Comparator antigen` %in% c("Beta")))
  nrow_delta <-  nrow(data %>% filter(`Comparator antigen` %in% c("Delta")))
  
  base_t <- 0
  spacer <- 1
  
  if(which_plot == "titer_drop") {
    plots <- split_up_forest_plot_drops(data, hline_by, row_label = TRUE, show_means = TRUE, axis_text_size)
    
    to_save <- paste0(path_to_save, to_save, "_fold_drops")
    height_a <- height #  + 4*spacer
    height_bcd <- height # -2 *spacer
    
    height_a <- height
    height_bcd <- height*(nrow_alpha+nrow_beta+nrow_delta)/nrow_wt #-2 *spacer
    
    # save for each subgroup
    if(single_plots & hline_by[1] == "standardise_encounters") {
      
      groupings <- as.character(data %>% pull(hline_by[1]) %>% unique())
      
      #standard height is 2x vax
      nrows_sub <- data %>% filter(`Comparator antigen` %in% c("WT","Wu-1", "WA1","D614G","B.1", "Wu-1?")) %>% group_by_at(hline_by[[1]]) %>% count()
      
      nrows_sub$height <- nrows_sub$n * 7.3/nrows_sub$n[nrows_sub$standardise_encounters == "2x Vax"]
      nrows_sub$height[nrows_sub$height < 1] <- 1
      lapply(groupings, function(group) {
        
        temp_data <- data[data[,hline_by[1]] == group,]
        save_in_scaled_format_encounter(temp_data, which_plot, width =8, height =nrows_sub$height[nrows_sub$standardise_encounters == group], to_save = paste0(to_save, "_", hline_by[1], "_", group), hline_by, single_plots = FALSE, axis_text_size = axis_text_size_p*1.6)
      })
    }
    
    
  } else {
    plots <- split_up_forest_plot_titer(data, hline_by, row_label = TRUE, show_means = TRUE, axis_text_size)
    
    to_save <- paste0(path_to_save, to_save, "_titers")
    height_a <- height +1
    height_bcd <- height +1
    
    # save for each subgroup
    if(single_plots & hline_by[1] == "standardise_encounters") {
      # height <- height 
      height_a <- height + 2
      height_bcd <- height + 2
      
      groupings <- as.character(data %>% pull(hline_by[1]) %>% unique())
      
      #standard height is 2x vax
      nrows_sub <- data %>% filter(`Comparator antigen` %in% c("WT","Wu-1", "WA1","D614G","B.1", "Wu-1?")) %>% group_by_at(hline_by[[1]]) %>% count()
      
      nrows_sub$height <- nrows_sub$n * 7/nrows_sub$n[nrows_sub$standardise_encounters == "2x Vax"]
      nrows_sub$height[nrows_sub$height < 1] <- 0.7
      lapply(groupings, function(group) {
        
        temp_data <- data[data[,hline_by[1]] == group,]
        save_in_scaled_format_encounter(temp_data, which_plot, width =8, height =nrows_sub$height[nrows_sub$standardise_encounters == group], to_save = paste0(to_save, "_", hline_by[1], "_", group), hline_by, single_plots = FALSE, axis_text_size = axis_text_size_p*1.6)
      })
    }
    
  }
  
  layout_wt <- c(area(t = base_t, l = 1, b = base_t+nrow_wt, r = 2))
  layout_alpha <- c(area(t = base_t, l = 1, b = base_t+nrow_alpha, r = 2))
  layout_beta <- c(area(t = base_t, l = 1, b = base_t+nrow_beta, r = 2))
  layout_delta <- c(area(t = base_t, l = 1, b = base_t+nrow_delta, r = 2))
  layout_bcd <- c(
    area(t = base_t, l = 1, b = base_t+nrow_alpha +spacer, r = 2),
    area(t = 2*spacer+(base_t+nrow_alpha), l = 1, b = 2*spacer+base_t+nrow_alpha+nrow_beta, r = 2),
    area(t = 3*spacer+(base_t+nrow_beta + nrow_alpha), l = 1, b = 3*spacer+base_t+nrow_alpha+nrow_beta+ nrow_delta, r = 2)
  )
  

    ggsave(paste0(to_save, "_WT.",fileext), plots$WT, width = width, height = height_a, dpi = 300)
  
}


reset_rowlabel <- function(data) {
  length_rowlabel <-min(unlist(lapply(data$rowlabel, function(x) nchar(strsplit(x, ")")[[1]][1]))))
  data$rowlabel <- unlist(lapply(data$rowlabel, function(x) substr(x, 1, length_rowlabel+3)))
  
  return(data)
}


# Huxtable format functions
format_huxtable_cols_standard <- function(ht, single_grouping = TRUE, outer_group = "") {
  base <- 1
  if(!single_grouping) {
    base <- 2
    ht %>%
      merge_cells(.,1:2,1) -> ht
    ht[1,1] <- outer_group
  }
  
  ht %>%
    merge_cells(.,1, (base+1):(base+6)) %>%
    merge_cells(.,1, (base+7):(base+11)) %>%
    merge_cells(.,1, (base+12):(base+15)) %>%
    merge_cells(.,1:2,base) -> ht
  
  ht[1,base] <- "Serum Group"
  ht[1,base+1] <- "GMT"
  ht[1,base+7] <- "Mean Omicron fold drop from"
  ht[1,base+12] <- "Mean fold drop to WT"
  ht[2,(base+1):(base+6)] <- c("WT", "Alpha", "Beta", "Gamma", "Delta","Omicron")
  ht[2,(base+7):(base+11)] <- c("WT", "Alpha", "Beta", "Gamma", "Delta")
  ht[2,(base+12):(base+15)] <- c("Alpha", "Beta", "Gamma", "Delta")
  
  ht %>% set_align(1, everywhere, "center") %>%
    set_align(-1, 1, ".") %>%
    set_all_padding(2) %>% 
    set_outer_padding(0.1) %>% 
    set_bold(row = 1, col = everywhere) %>% 
    set_bold(row = everywhere, col = 1) %>% 
    set_bottom_border(row = 1, col = everywhere) %>% 
    set_width(0.9) %>% 
    set_font_size(7) %>% 
    theme_article() %>%
    set_right_border(everywhere, base, 1) %>% 
    set_right_border(everywhere, base+6, 1) %>% 
    set_right_border(everywhere, base+11, 1) %>% 
    set_right_border_color(everywhere, base, "grey")  %>%
    set_right_border_color(everywhere, base+6, "grey") %>%
    set_right_border_color(everywhere, base+11, "grey") -> ht
  
  return(ht)
}

format_huxtable_row <- function(ht, data, column_name, rev = FALSE) {
  
  nr_rows <- table(data[,column_name])
  nr_rows <- nr_rows[nr_rows > 0]
  if(rev) {
    nr_rows <- rev(nr_rows)
  }
  nr_rows <- c(0, nr_rows)

  start <- 3
  for(i in 2:(length(nr_rows))) {
    start_pos <-start+sum(nr_rows[(i-1):1])
    end_pos <-(start-1)+sum(nr_rows[i:1])
 
    ht <- ht %>%
      merge_cells(.,start_pos:end_pos,1) %>%
      set_bottom_border(row =end_pos, col = everywhere, value = 0.2)
  
  }
  
  return(ht)
}
