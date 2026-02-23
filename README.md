# Name of App 

KM Survival App 

Github repository: 

https://github.com/meredithspalmer/MoveApps_Survival

## Description

Perform basic Kaplan-Meier survival analyses. 

## Documentation

This app implements fundamental Kaplan-Meier (KM) survival estimation functions, including producing life tables, survival curves, cumulative hazard curves, and per-group estimation. 

Users can define a study period, censor for post-capture mortality, and indicate how null timestamp data should be handled. 

Users can also indicate whether they want additional analyses comparing key groups of interest, currently, segregated by sex, lifestage, reproductive condition, or attachment style. 

Data are cleaned (removing empty locations, marked outliers, marked test data, start dates occurring after end dates), dates are processed occurring to user inputs (crop to study period and format any missing start/end dates). 

Next, data are summarized into a table containing the duration (start and end dates, to allow for staggered entry) of each individual and a survival event indicator is derived (1 if mortality occurs during the study period, 0 if individual survives or data are censored, e.g., collar or individual is lost). 

Basic summaries are producing, including graphs and tables detailing each individual's survival tenure during the study period. 

A Kaplan-Meier survival estimation is run allowing for staggered entry. Outputs include a life table, a KM survival curve, and a cumulative hazard plot. 

Users have the option to select groups for comparison: if this option is selected, the app will also produce a table of log-rank test outputs along with a plot contrasting each survival curve. 


### Application scope

#### Generality of App usability

The app will produce a warning and terminate if none of the individuals in the study experienced a mortality event during the study period. 

This app allows for staggered entry during the defined study period. 

#### Required data properties

**Events**: This app can only be used if mortality information is captured in one of the following columns: death_comments, deployment_end_comments, deployment_end_type, mortality_location_filled. 

**Sample size**: The required sample size for a KM survival analysis depends primarily on the total number of deaths, rather than the total number of individuals. In general, the lower the event (death) rate, the higher the number of required individuals. This app does NOT perform a power analysis prior to performing the survival analyzes. However, if fewer than 10 mortality events are detected, the app will generate a warning that the model may have low statistical power, potentially resulting in unreliable estimates and poor predictive power. 

Note that a larger sample size is required for comparison across groups, 

### Input type

`move2::move2_loc`

### Output type

`move2::move2_loc`

### Artefacts

*Tracking history:* `tracking_history.png`: Figure (png) detailing the start and end dates of each individual across the study and during the tracking period, along with an indicator of how each individual was terminated (death, censored, survived). 

*Life table:* `survival_statistics.csv`: Output of KM survival analysis; table (csv) with the time, number of individuals at risk, number of events, survival, standard error, and upper & lower 95% confidence intervals. 

*KM survival curve:* `km_survival_curve.png`: KM survival curve plot (png), detailing median survival time. 

*Cumulative hazard plot:* ... includes number of individuals at risk per time period and number of events per time period. 

*Log-rank test:* `logrank_table_statistics.csv`: Output of comparing survival curves between groups; table (csv) with test statistics, degrees of freedom, p-value, and pairwise comparisons. 

*Comparison curves:*:  .... 


### Settings 

*Example:* `Radius of resting site` (radius): Defined radius the animal has to stay in for a given duration of time for it to be considered resting site. Unit: `metres`.

*Please list and define all settings/parameters that the App requires to be set by the App user, if necessary including their unit. Please first state the Setting name the user encounters in the Settings menu defined in the appspecs.json, and between brackets the argument used in the R function to be able to identify it quickly in the code if needed.*

*Example:* `Radius of resting site` (radius): Defined radius the animal has to stay in for a given duration of time for it to be considered resting site. Unit: `metres`.

### Most common errors

Please send errors to mspalmer.zool@gmail.com 

### Null or error handling

*Please indicate for each setting as well as the input data which behaviour the App is supposed to show in case of errors or NULL values/input. Please also add notes of possible errors that can happen if settings/parameters are improperly set and any other important information that you find the user should be aware of.*

**Setting `radius`:** If no radius AND no duration are given, the input data set is returned with a warning. If no radius is given (NULL), but a duration is defined then a default radius of 1000m = 1km is set. 
