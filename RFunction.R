library(move2)
library(survival)
library(survminer)
library(ggplot2)   
library(dplyr)
library(lubridate)
library(stringr)
library(sf)
library(forcats)

# logger.fatal(), logger.error(), logger.warn(), logger.info(), logger.debug(), logger.trace()

# Survival Function 
rFunction = function(data, sdk, time_period_start, time_period_end, fix_na_start_times, fix_na_end_times, censor_capture_mortality, group_comparison_individual) {
  
  logger.info(paste("Welcome to the", sdk))
  
  ## Cleaning and cropping ----------------------------------------------------
  
  data <- dplyr::filter(data, !sf::st_is_empty(data))       # Exclude empty locations
  data <- mt_filter_unique(data)                            # Exclude marked outliers 
  data <- data %>% filter_track_data(is_test == FALSE)      # Exclude data marked "test"
  
  
  ## Aggregate across multiple deployments (where present) ---
  
  # Extract event-level data 
  events <- data |>
    as_tibble() |>
    dplyr::select(any_of(c(
      "deployment_id",
      "individual_local_identifier",   
      "timestamp")))
  
  # Extract relevant track-level attributes
  tracks <- mt_track_data(data) |>
    mutate(mortality_location_filled = if_else(
      is.na(mortality_location) | st_is_empty(mortality_location),
      0L, 1L)) |> 
    dplyr::select(any_of(c(
      "individual_local_identifier",
      "deployment_id",
      "individual_id",
      "sex",
      "animal_life_stage",
      "animal_reproductive_condition",
      "attachment_type",
      "study_site",
      "deploy_on_timestamp",
      "deploy_off_timestamp",
      "deployment_end_type",
      "deployment_end_comments",
      "death_comments",
      "mortality_location_filled", 
      "mortality_date",
      "mortality_type", 
      "timestamp_first_deployed_location",
      "timestamp_last_deployed_location")))
  
  # Join track attributes to every event row
  use_deployment_join <- all(c("deployment_id") %in% names(events), 
                             "deployment_id" %in% names(tracks)) &&
    any(!is.na(events$deployment_id)) &&
    any(!is.na(tracks$deployment_id))
  
  if (use_deployment_join) {
    events_with_ind <- events |>
      left_join(tracks, by = "deployment_id")
    
  } else {
    if (!"individual_local_identifier" %in% names(events)) {
      logger.fatal("Cannot join: neither deployment_id nor individual_local_identifier is available in events")
    }
    logger.info("Joining on individual_local_identifier (deployment_id join not possible)")
    
    events_with_ind <- events |>
      left_join(tracks, by = "individual_local_identifier")
  }
  
  events_with_ind <- events_with_ind |>
    relocate(any_of(c("individual_id", "individual_local_identifier", "deployment_id",
                      "timestamp")),
             .before = everything())
  
  # Summarize timestamps and location count per individual
  summary_table <- events_with_ind |>
    group_by(individual_id, individual_local_identifier) |>
    summarise(
      first_timestamp = min(as.Date(timestamp), na.rm = TRUE),
      last_timestamp  = max(as.Date(timestamp), na.rm = TRUE),
      n_locations     = n(),
      n_deployments   = if ("deployment_id" %in% names(events_with_ind)) {
        n_distinct(deployment_id, na.rm = TRUE)
      } else {
        1L  
      },
      
      # Timestamp columns: min / max if present
      timestamp_first_deployed_location = if ("timestamp_first_deployed_location" %in% names(events_with_ind))
        min(timestamp_first_deployed_location, na.rm = TRUE) else NA,
      
      timestamp_last_deployed_location = if ("timestamp_last_deployed_location" %in% names(events_with_ind))
        max(timestamp_last_deployed_location, na.rm = TRUE) else NA,
      
      deploy_on_timestamp = if ("deploy_on_timestamp" %in% names(events_with_ind)) {
        if (all(is.na(deploy_on_timestamp))) as.POSIXct(NA) else min(deploy_on_timestamp, na.rm = TRUE)
      } else as.POSIXct(NA),
      
      deploy_off_timestamp = if ("deploy_off_timestamp" %in% names(events_with_ind)) {
        if (all(is.na(deploy_off_timestamp))) as.POSIXct(NA) else max(deploy_off_timestamp, na.rm = TRUE)
      } else as.POSIXct(NA),
      
      # Mortality column: 1/0 if filled if present 
      mortality_location_filled = if ("mortality_location_filled" %in% names(events_with_ind))
        as.integer(any(mortality_location_filled >= 1, na.rm = TRUE)) else NA_integer_,
      
      # Categorical columns: collapsed unique if present
      sex = if ("sex" %in% names(events_with_ind))
        str_c(unique(sex[!is.na(sex)]), collapse = " | ") else NA_character_,
      
      death_comments = if ("death_comments" %in% names(events_with_ind))
        str_c(unique(death_comments[!is.na(death_comments)]), collapse = " | ") else NA_character_,
      
      deployment_end_comments = if ("deployment_end_comments" %in% names(events_with_ind))
        str_c(unique(deployment_end_comments[!is.na(deployment_end_comments)]), collapse = " | ") else NA_character_,
      
      deployment_end_type = if ("deployment_end_type" %in% names(events_with_ind))
        str_c(unique(deployment_end_type[!is.na(deployment_end_type)]), collapse = " | ") else NA_character_,
      
      animal_life_stage = if ("animal_life_stage" %in% names(events_with_ind))
        str_c(unique(animal_life_stage[!is.na(animal_life_stage)]), collapse = " | ") else NA_character_,
      
      animal_reproductive_condition = if ("animal_reproductive_condition" %in% names(events_with_ind))
        str_c(unique(animal_reproductive_condition[!is.na(animal_reproductive_condition)]), collapse = " | ") else NA_character_,
      
      attachment_type = if ("attachment_type" %in% names(events_with_ind))
        str_c(unique(attachment_type[!is.na(attachment_type)]), collapse = " | ") else NA_character_,
      
      .groups = "drop"
    ) |>
    mutate(
      
      # Clean empty strings (fill NA) for columns that exist
      across(
        any_of(c(
          "death_comments",
          "deployment_end_comments",
          "deployment_end_type",
          "animal_life_stage",
          "animal_reproductive_condition"
        )),
        ~ if_else(. == "", NA_character_, .)
      ),
      
      # Convert deploy timestamps 
      across(
        any_of(c("deploy_on_timestamp", "deploy_off_timestamp")),
        as.Date
      )
    )
  
  
  ## Clean dates ---
  
  # Start times 
  if(fix_na_start_times == "timestamp"){
    summary_table <- summary_table %>% 
      mutate(missing_timestamp_start = is.na(deploy_on_timestamp))
    n_missing <- sum(is.na(summary_table$deploy_on_timestamp), na.rm = TRUE)
    
    summary_table <- summary_table %>%
      mutate(deploy_on_timestamp = if_else(
        is.na(deploy_on_timestamp),
        as.Date(first_timestamp),
        deploy_on_timestamp)) %>% 
      dplyr::select(-missing_timestamp_start)
    
    if (n_missing > 0) {
      log.info(sprintf("Warning: Replaced %d missing deploy_on_timestamp value%s with first_timestamp.",
                       n_missing,
                       if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  }
  
  if(fix_na_start_times == "remove"){
    n_missing <- sum(is.na(summary_table$deploy_on_timestamp))
    summary_table <- summary_table %>% filter(!is.na(deploy_on_timestamp))
    
    if (n_missing > 0) {
      log.info(sprintf("Warning: Removed %d deploy_on_timestamp value%s that were NA.", n_missing,
                       if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  }
  
  # End times 
  if(fix_na_end_times == "timestamp"){
    summary_table <- summary_table %>%
      mutate(missing_timestamp_end = is.na(deploy_off_timestamp))
    n_missing <- sum(is.na(summary_table$deploy_off_timestamp), na.rm = TRUE)
    
    summary_table <- summary_table %>%
      mutate(deploy_off_timestamp = if_else(
        is.na(deploy_off_timestamp),
        as.Date(last_timestamp),
        deploy_off_timestamp)) %>% 
      dplyr::select(-missing_timestamp_end)
    
    if (n_missing > 0) {
      log.info(sprintf("Warning: Replaced %d missing deploy_off_timestamp value%s with last_timestamp.", 
                       n_missing,
                       if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  }
  
  if(fix_na_end_times == "systime"){
    summary_table <- summary_table %>%
      mutate(missing_timestamp_end = is.na(deploy_off_timestamp))
    n_missing <- sum(is.na(summary_table$deploy_off_timestamp), na.rm = TRUE)
    
    summary_table <- summary_table %>%
      mutate(deploy_off_timestamp = if_else(
        is.na(deploy_off_timestamp),
        Sys.Date(), 
        deploy_off_timestamp))%>% 
      dplyr::select(-missing_timestamp_end)
    
    if (n_missing > 0) {
      log.info(sprintf("Warning: Replaced %d missing deploy_off_timestamp value%s with current date.",
                       n_missing,
                       if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  }
  
  if(fix_na_end_times == "remove"){
    n_missing <- sum(is.na(is.na(summary_table$deploy_off_timestamp)))
    summary_table <- summary_table %>% filter(!is.na(deploy_off_timestamp))
    
    if (n_missing > 0) {
      log.info(sprintf("Warning: Removed %d deploy_off_timestamp and/or deploy_on_timestamp value%s that were NA.", n_missing,
                       if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  }
  
  # Remove data for individuals where "deploy_off_timestamp" occurs before "deploy_on_timestamp" 
  n_original <- nrow(summary_table) 
  summary_table <- summary_table %>%
    filter(deploy_off_timestamp >= deploy_on_timestamp)
  n_removed <- n_original - nrow(summary_table)  
  
  if (n_removed > 0) {
    log.info(sprintf("Warning: Removed %d individual%s where deploy_off_timestamp < deploy_on_timestamp.",
                    n_removed, if (n_removed == 1) "" else "s"),
            call. = FALSE, immediate. = TRUE)
  }
  
  
  ## Crop data to user-defined temporal windows ---
  
  # Removed censored data (mortalities within set period of capture)
  if(censor_capture_mortality > 0){
    n_before <- nrow(summary_table)
    
    summary_table <- summary_table %>%
      mutate(raw_deploy_on_timestamp = deploy_on_timestamp) %>%
      mutate(censor_cutoff = deploy_on_timestamp + lubridate::days(censor_capture_mortality)) %>%
      mutate(remove_due_to_early_end = !is.na(deploy_off_timestamp) & deploy_off_timestamp <= censor_cutoff) %>%
      filter(!remove_due_to_early_end) %>%
      mutate(deploy_on_timestamp = censor_cutoff) %>%
      select(-censor_cutoff, -remove_due_to_early_end)
    
    n_after  <- nrow(summary_table)
    n_removed <- n_before - n_after
    
    if (n_removed > 0) {
      log.info(paste0("Warning: Removed ", n_removed, " individual(s) because deploy_off_timestamp occurred within ", censor_capture_mortality, " day(s) after deploy_on_timestamp"),
               call. = FALSE, immediate. = TRUE)
    } 
  }
  
  # Crop to study period of interest 

  # Save original deploy_off_time 
  summary_table <- summary_table %>% mutate(raw_deploy_off_timestamp = deploy_off_timestamp) 
  
  # Define window 
  effective_start <- if (is.na(time_period_start)) {
    min(summary_table$deploy_on_timestamp, na.rm = TRUE)
  } else {
    time_period_start
  }
  
  effective_end <- if (is.na(time_period_end)) {
    max(summary_table$deploy_off_timestamp, na.rm = TRUE)
  } else {
    time_period_end
  }  
  
  # Run updates 
  if(!is.na(time_period_start) | !is.na(time_period_end)){
    
    # Crop to window 
    n_original <- nrow(summary_table) 
    summary_table <- summary_table %>%
      
      # Determine if the deployment overlaps study window 
      mutate(overlaps_study = deploy_on_timestamp <= effective_end & 
               deploy_off_timestamp  >= effective_start) %>%
      filter(overlaps_study | is.na(overlaps_study)) %>%    
      
      # Crop to window 
      mutate(first_timestamp = pmax(deploy_on_timestamp, effective_start, na.rm = TRUE),
             last_timestamp  = pmin(deploy_off_timestamp, effective_end,   na.rm = TRUE)) %>%
      
      # Clean 
      select(-overlaps_study) 
    
    n_removed <- n_original - nrow(summary_table)
    if (n_removed > 0) {
      log.info(sprintf("Warning: %d record%s did not overlap the user-defined study window and were removed.",
                      n_removed, if (n_removed == 1) "" else "s"),
              call. = FALSE, immediate. = TRUE)
    }
  }
  
  
  ## Calculate entry time and exit time (for staggered entry) ---
  origin_date <- if_else(is.na(time_period_start),
                         min(summary_table$deploy_on_timestamp, na.rm = TRUE),
                         time_period_start)
  
  summary_table <- summary_table %>%
    mutate(origin_date = origin_date, 
           entry_time_days  = as.numeric(difftime(deploy_on_timestamp, origin_date, units = "days")),
           exit_time_days   = as.numeric(difftime(deploy_off_timestamp, origin_date, units = "days")))
  
  
  ## Calculate mortality indicator ---
  # Here, event = 1 if observed death, 0 if censored or survived 
  # death comments to flag 
  positive_pattern <- "dead|death|cod|predation|predator|vehicle|collision|killed|poach|shot|hunt|harvest|mortality"
  
  # search in data 
  summary_table <- summary_table %>%
    
    # Initialize mortality event 
    mutate(mortality_event = NA_real_) %>%
    
    # Identify survivors (individuals who last beyond study)
    mutate(
      survived_beyond_study = !is.na(raw_deploy_off_timestamp) &
        raw_deploy_off_timestamp > as.Date(effective_end),
      
      mortality_event = if_else(survived_beyond_study, 0L, mortality_event),
      
      # Update columns to remove ambiguity (e.g., if animal dies after study window)
      death_comments = if ("death_comments" %in% names(.)) {
        if_else(survived_beyond_study, "survived beyond study", death_comments)
      } else death_comments,
      
      deployment_end_comments = if ("deployment_end_comments" %in% names(.)) {
        if_else(survived_beyond_study, "survived beyond study", deployment_end_comments)
      } else deployment_end_comments,
      
      deployment_end_type = if ("deployment_end_type" %in% names(.)) {
        if_else(survived_beyond_study, "survived beyond study", deployment_end_type)
      } else deployment_end_type,
      
      mortality_location_filled = if ("mortality_location_filled" %in% names(.)) {
        if_else(survived_beyond_study, 0L, mortality_location_filled)
      } else mortality_location_filled
    ) %>%
    
    # Search for mortality indicators 
    # A. death_comments keywords
    mutate(
      mortality_event = case_when(
        "death_comments" %in% names(.) & str_detect(tolower(death_comments), "\\bnot\\b") ~ 0L,
        "death_comments" %in% names(.) & str_detect(tolower(death_comments), positive_pattern) ~ 1L,
        mortality_event == 1L ~ 1L,
        TRUE ~ mortality_event
      )
    ) %>% 
    
    # B. mortality_location_filled
    mutate(
      mortality_event = case_when(
        "mortality_location_filled" %in% names(.) &
          mortality_location_filled >= 1 ~ 1L,
        
        mortality_event == 1L ~ 1L,
        TRUE ~ mortality_event
      )
    ) %>%
    
    # C. deployment_end_type  
    mutate(
      mortality_event = case_when(
        mortality_event == 1L ~ 1L,
        
        # Mortality indication
        "deployment_end_type" %in% names(.) &
          str_detect(tolower(deployment_end_type), "\\bdead\\b|\\bdeath\\b") ~ 1L,
        
        # Censoring indication
        "deployment_end_type" %in% names(.) &
          tolower(deployment_end_type) %in% c("removal", "other", "unknown", "survived beyond study") ~ 0L,
        
        # Missing column OR NA value → censored 
        (!"deployment_end_type" %in% names(.) | is.na(deployment_end_type)) &
          is.na(mortality_event) ~ 0L,
        
        TRUE ~ mortality_event
      )
    ) %>%
    
    # Final censoring: remaining NA → 0 
    mutate(
      mortality_event = if_else(
        is.na(mortality_event) & !is.na(deploy_off_timestamp),
        0L,
        mortality_event
      )
    ) %>%
    
    # Clean up & relocate
    select(-survived_beyond_study) %>%
    relocate(mortality_event, .after = deployment_end_type)
  
  # Error out: No deaths 
  n_mort_events <- sum(summary_table$mortality_event == 1, na.rm = TRUE)
  if (n_mort_events == 0) {
    logger.fatal("Cannot run survival analysis: no mortality events detected.",
            call. = FALSE, immediate. = TRUE)
  }
  
  # Warning: Small proportion of deaths 
  if (n_mort_events <= 10) {
    logger.warn(sprintf("Few (%d) deaths detected. Model may have low statistical power, potentially resulting in unreliable estimates and poor predictive power.", n_mort_events),
            call. = FALSE, immediate. = TRUE)
  }
  
  
  ## Clean user-defined attributes for comparison (if applicable) --- 
  if(group_comparsion_individual == "sex"){
    
    # Remove NAs 
    n_original <- nrow(summary_table)
    summary_table <- summary_table[!is.na(summary_table$sex),] 
    
    # Get unique sexes after cleaning
    unique_sexes <- sort(unique(summary_table$sex))
    n_sexes <- length(unique_sexes)
    
    logger.info(sprintf("%d sexes detected after cleaning: %s", n_sexes, 
                    paste(unique_sexes, collapse = ", ")),
            call. = FALSE, immediate. = TRUE)
    
    n_lost <- n_original - nrow(summary_table)
    if (n_lost > 0) {
      logger.warn(sprintf("%d individuals with NA sex removed from study.", n_lost),
              call. = FALSE, immediate. = TRUE)
    }
  }
  
  if(group_comparsion_individual == "lifestage"){
    
    # Clean data, remove NA rows
    n_original <- nrow(summary_table)
    summary_table <- summary_table %>%
      filter(!is.na(animal_life_stage)) %>%
      mutate(
        animal_life_stage = str_trim(animal_life_stage),
        animal_life_stage = str_replace_all(animal_life_stage, "\\s+", ""),
        animal_life_stage = str_extract(animal_life_stage, "^[^|]+"),
        animal_life_stage = str_replace(animal_life_stage, "–", "-"))
    
    # Get unique life-stages after cleaning
    unique_stages <- sort(unique(summary_table$animal_life_stage))
    n_life_stages <- length(unique_stages)
    
    logger.info(sprintf("%d life-stages detected after cleaning: %s", n_life_stages, 
                        paste(unique_stages, collapse = ", ")),
                call. = FALSE, immediate. = TRUE)
    
    n_lost <- n_original - nrow(summary_table)
    if (n_lost > 0) {
      logger.warn(sprintf("%d individuals with NA life stage removed from study.", n_lost),
              call. = FALSE, immediate. = TRUE)
    }
  }
  
  if(group_comparsion_individual == "reproCond"){
  
    # Clean data, remove NAs 
    n_original <- nrow(summary_table)
    summary_table <- summary_table %>%
      filter(!is.na(animal_reproductive_condition)) %>%
      mutate(
        animal_reproductive_condition = str_trim(animal_reproductive_condition),
        animal_reproductive_condition = str_replace_all(animal_reproductive_condition, "\\s+", ""),
        animal_reproductive_condition = str_extract(animal_reproductive_condition, "^[^|]+"),
        animal_reproductive_condition = str_replace(animal_reproductive_condition, "–", "-"))
    
    # Get unique conditions after cleaning
    unique_conditions <- sort(unique(summary_table$animal_reproductive_condition))
    n_conditions <- length(unique_conditions)
    
    logger.info(sprintf("%d reproductive conditions detected after cleaning: %s", n_conditions, 
                        paste(unique_conditions, collapse = ", ")),
                call. = FALSE, immediate. = TRUE)
    
    n_lost <- n_original - nrow(summary_table)
    if (n_lost > 0) {
      logger.warn(sprintf("%d individuals with NA reproductive conditions removed from study.", n_lost),
              call. = FALSE, immediate. = TRUE)
    }
  }
  
  if(group_comparsion_individual == "attachment"){
    
    # Clean data, remove NAs 
    n_original <- nrow(summary_table)
    summary_table <- summary_table %>%
      filter(!is.na(attachment_type)) %>%
      mutate(
        attachment = str_trim(attachment_type),
        attachment = str_replace_all(attachment_type, "\\s+", ""),
        attachment = str_extract(attachment_type, "^[^|]+"),
        attachment = str_replace(attachment_type, "–", "-"))
    
    n_attaches <- length(unique(summary_table$attachment))
    logger.info(sprintf("%d attachment types detected.", n_attaches), call. = FALSE, immediate. = TRUE)
    
    n_lost <- n_original - nrow(summary_table)
    if (n_lost > 0) {
      logger.warn(sprintf("%d individuals with NA attachment type removed from study.", n_lost), 
              call. = FALSE, immediate. = TRUE)
    }
  } 
  
  
  ## Basic summaries of data --------------------------------------------------
  
  # Plot each individual's tracking history (across entire dataset) --- 
  deployment_to_ind <- mt_track_data(data) |>
    dplyr::select(deployment_id, individual_id) |>
    distinct()
  
  data_with_ind <- data |>
    mt_as_event_attribute(c("individual_id", "deployment_id"), .keep = FALSE)
  
  track_times <- data_with_ind |>
    group_by(individual_id) |>
    summarise(
      start     = min(timestamp, na.rm = TRUE),
      end       = max(timestamp, na.rm = TRUE),
      n_locs    = n(),
      n_deploy  = n_distinct(deployment_id),     
      .groups   = "drop"
    ) |>
    left_join(
      mt_track_data(data) |>
        distinct(individual_id, .keep_all = TRUE),
      by = "individual_id"
    ) |>
    mutate(
      duration_days = round(as.numeric(difftime(end, start, units = "days")), 1),
      track_label   = fct_reorder(as.character(individual_id), start)
    ) |>
    arrange(start)
  
  (tracking_history <- ggplot(track_times) +
      geom_segment(aes(x = start, xend = end, y = track_label, yend = track_label),
                   linewidth = 3, color = "steelblue") +
      geom_point(aes(x = start, y = track_label), color = "darkgreen", size = 3.5) +
      geom_point(aes(x = end,   y = track_label), color = "firebrick",  size = 3.5) +
      labs(title    = "Individual Tracking History (Full Data)",
           subtitle = sprintf("%d unique individuals • %d total deployments • %d locations",
                              nrow(track_times),
                              sum(track_times$n_deploy, na.rm = TRUE),
                              sum(track_times$n_locs, na.rm = TRUE)),
           x = "Time",
           y = "Individual ID") +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_text(size = 8),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title   = element_text(face = "bold", size = 14), 
            plot.subtitle = element_text(size = 12, color = "grey50"), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      scale_x_datetime(date_breaks = "1 year",
                       date_labels = "%Y"))
  
  # want to figure out how to add width, height, units, dpi, bg to artifact 
  artifact <- appArtifactPath("tracking_history_full.png")
  logger.info(paste("Saving tracking history (full) plot:", artifact))
  png(artifact)
  tracking_history
  dev.off()

  # Plot each individual's tracking history (during study period) ---
  track_times <- summary_table %>%
    mutate(start         = first_timestamp,
           end           = last_timestamp,
           duration_days = round(as.numeric(difftime(end, start, units = "days")), 1),
           track_label   = fct_reorder(as.character(individual_id), start)) %>%
    arrange(start)
  
  (tracking_history_subset <- ggplot(track_times) +
      geom_segment(aes(x = start, xend = end, y = track_label, yend = track_label),
                   linewidth = 3, color = "steelblue") +
      geom_point(aes(x = start, y = track_label),
                 color = "darkgreen", size = 3.5) +
      geom_point(aes(x = end, y = track_label),
                 color = "firebrick", size = 3.5) +
      labs(title = "Individual Tracking History (Data Subset)",
           subtitle = sprintf("%d unique individuals • %d total deployments • %d locations",
                              nrow(track_times),
                              sum(track_times$n_deploy, na.rm = TRUE),
                              sum(track_times$n_locations, na.rm = TRUE)),  
           x = "Time",
           y = "Individual ID") +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_text(size = 8),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title   = element_text(face = "bold", size = 14), 
            plot.subtitle = element_text(size = 12, color = "grey50"), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      scale_x_datetime(date_breaks = "1 year",
                       date_labels = "%Y")) 
  
  
  # want to figure out how to add width, height, units, dpi, bg to artifact 
  artifact <- appArtifactPath("tracking_history_cropped.png")
  logger.info(paste("Saving tracking history (cropped) plot:", artifact))
  png(artifact)
  tracking_history_subset
  dev.off()
  
  
  ## Survival Analysis ----------------------------------------------------------
  
  ## Fit Kaplan-Meier with staggered entry --- 
  km_fit <- survfit(Surv(entry_time_days, exit_time_days, mortality_event) ~ 1, 
                    data = summary_table)
  
  
  ## Life table --- 
  times <- round(seq(min(summary_table$entry_time_days), max(summary_table$exit_time_days), 
                     length.out = 10))
  s <- summary(km_fit, times = times)
  life_table <- data.frame(time_days        = s$time,
                           n_risk           = s$n.risk,
                           n_event          = s$n.event,
                           survival_prob    = s$surv,
                           std_err          = s$std.err,
                           lower_95         = s$lower,
                           upper_95         = s$upper)
  
  # unsure if the row.names will cause issues - double check 
  write.csv(life_table, file = appArtifactPath("life_table.csv", row.names = F))
  
  
  ## Plotting statistics --- 
  n.ind <- nrow(summary_table)
  n.events <- nrow(summary_table[summary_table$mortality_event == 1,])
  n.days <- as.numeric(summary(km_fit)$table["median"])
  
  ## KM Survival Curve ---
  km_curve <- ggsurvplot(
    km_fit,
    data = summary_table,
    title = "Kaplan-Meier Survival Curve",
    subtitle = paste0("N = ", n.ind, ", Events = ", n.events, ", Median Survival = ", 
                      med$median, " days"),
    xlab = "Time (days)",
    ylab = "Survival Probability",
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.title = "Number at Risk",
    risk.table.y.text = FALSE,     
    risk.table.height = 0.18,
    conf.int = TRUE,
    censor.shape = "|",
    censor.size = 3,
    legend = "none",
    pval = TRUE,
    surv.median.line = "hv",        
    palette = c("#E69F00", "#56B4E9"),
    ggtheme = theme_classic(base_size = 12) + 
      theme(plot.title         = element_text(face = "bold", size = 14), 
            plot.subtitle      = element_text(size = 12, color = "gray50"),
            axis.text          = element_text(color = "black"),
            panel.grid.major.y = element_line(color = "gray90"), 
            panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
            plot.margin        = margin(10, 10, 10, 10)))
          
  # want to figure out how to add width, height, units, dpi, bg to artifact 
  artifact <- appArtifactPath("km_survival_curve.png")
  logger.info(paste("Saving Kaplan-Meier survival curve plot:", artifact))
  png(artifact)
  km_curve
  dev.off()
  
  ## Cumulative hazard plot ----
  cum_hazard <- ggsurvplot(
    km_fit,
    fun = "cumhaz",
    conf.int = TRUE,
    risk.table = TRUE,
    cumevents = TRUE,                 
    tables.height = 0.18,              
    tables.y.text = FALSE,            
    surv.median.line = "hv",         
    pval = TRUE,                     
    xlab = "Time (days)",
    ylab = "Cumulative Hazard",
    title = "Cumulative Hazard",
    subtitle = paste0("N = ", n.ind, ", Events = ", n.events), #update events 
    palette = c("#E69F00", "#56B4E9"),
    ggtheme = theme_classic(base_size = 12) + 
      theme(plot.title   = element_text(face = "bold", size = 14), 
            plot.subtitle = element_text(size = 12, color = "gray50"),
            axis.text    = element_text(color = "black"),
            panel.grid.major.y = element_line(color = "gray90"), 
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            plot.margin  = margin(10, 10, 10, 10)))
  
  # want to figure out how to add width, height, units, dpi, bg to artifact 
  artifact <- appArtifactPath("cumulative_hazard_plot.png")
  logger.info(paste("Saving cumulative hazard plot:", artifact))
  png(artifact)
  cum_hazard
  dev.off()
  
  
  
  
  result <- if (any(lubridate::year(move2::mt_time(data)) == year)) { 
    data[lubridate::year(move2::mt_time(data)) == year,]
  } else {
    NULL
  }
  
  if (!is.null(result)) {
    # Showcase creating an app artifact. 
    # This artifact can be downloaded by the workflow user on Moveapps.
    artifact <- appArtifactPath("plot.png")
    logger.info(paste("plotting to artifact:", artifact))
    png(artifact)
    plot(result[move2::mt_track_id_column(result)], max.plot=1)
    dev.off()
  } else {
    logger.warn("nothing to plot")
  }
  
  # Showcase to access a file ('auxiliary files') that is 
  # a) provided by the app-developer and 
  # b) can be overridden by the workflow user.
  fileName <- getAuxiliaryFilePath("auxiliary-file-a")
  logger.info(readChar(fileName, file.info(fileName)$size))

  # Pass original to the next app in the MoveApps workflow
  return(data)
}
