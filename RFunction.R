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
  
  ## Basic cleaning ---
  
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
        "death_comments" %in% names(.) &
          str_detect(tolower(death_comments),
                     "dead|death|cod|predation|predator|vehicle|collision|killed|poach|shot|hunt|harvest") ~ 1L,
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
  
  ADD ATTRIBUTE CODE 
  
  
  
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
