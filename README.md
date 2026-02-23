# KM Survival App 

Github repository: 

https://github.com/meredithspalmer/MoveApps_Survival

## Description

Perform basic Kaplan-Meier survival analyses and optional group comparisons via the log-rank test

## Documentation

This app implements fundamental Kaplan-Meier (KM) survival estimation functions. It produces life tables, survival curves, cumulative hazard curves, and, if applicable, statistical comparisons of per-group survival estimation.

**Kaplan-Meier Survival Estimation:** The KM estimator is a non-parametric method used to estimate the survival function, that is, the probability that an individual survives past time t, from lifetime data. This analysis allows for:
- *Right-censoring*, where the exact time of death is unknown for some individuals because they are still alive at the end of the study period, lost to follow-up (e.g., collar failure), or exit the study period alive for other reasons. 
- *Staggered entry* (also called left truncation or delayed entry), where individuals enter the study at different times rather than all starting at the same baseline. 

**Log-rank Test:** Users can request log-rank tests (Mantel–Cox test) to compare survival distributions across groups. The log-rank test assesses whether there are statistically significant differences in survival between two or more independent groups. Currently supported grouping variables are: sex, life stage, reproductive condition, or tag/collar attachment type. 

Users define a study period, can censor data to exclude post-capture mortality events, and specify how to handle missing timestamp information.

Data pre-processing includes:  
- Removing or updating empty or invalid data (according to user specification)
- Flagging and handling marked outliers and test data  
- Checking for logical errors (e.g., start date after end date)  
- Subsetting data to study period or select individuals (according to user specification)

The app then summarizes the data into a per-individual table and figure containing:  
- Duration within the study period (accounting for staggered entry and censoring)  
- Survival event indicator (1 = death occurred during study period, 0 = censored or survived)  

Finally, the app performs Kaplan-Meier survival estimation on the cleaned dataset, generating: 
- Summary statistics 
- A KM survival curve (estimated survival probability over time, with 95% confidence intervals)
- A cumulative hazard plot (total accumulated risk over time, with 95% confidence intervals)

When group comparison is enabled, the app additionally produces:  
- A table of log-rank test results  
- A comparative survival curve plot contrasting the groups


### Application scope

#### Generality of App usability

This app has currently been tested on mammals and birds with N_individuals > 100, N_deaths > 50, and study duration >1 year, but should be applicable to any dataset containing mortality data of sufficient sample size (see below). 

This app allows for staggered entry during the defined study period and for censored data (individuals that are "lost" from the study, e.g., due to equipment failure and individuals that survive the study period). 


#### Required data properties

**Events**: This app can only be used if mortality information (indication of event along with an associated end date) is captured in one of the following columns: death_comments, deployment_end_comments, deployment_end_type, mortality_location_filled. 

The app will produce a warning and terminate if none of the individuals in the study experienced a mortality event during the study period. 

**Sample size**: The required sample size for a KM survival analysis depends primarily on the total number of deaths, rather than the total number of individuals. In general, the lower the event (death) rate, the higher the number of required individuals. 

This app does *not* perform a power analysis prior to performing the survival analyzes. However, if fewer than 10 mortality events are detected, the app will generate a warning that the model may have low statistical power, potentially resulting in unreliable estimates and poor predictive power. 

Note that a larger sample size is required for comparison across groups, 

**Comparison types:** This app enables the user to select one of four grouping options to compare across. These currently include sex, lifestage, reproductive condition, and collar attachment type. These grouping classifications are cleaned and standardized for capitalization and white-space errors (e.g., "male" and "Male" will register as the same class) but note that other mispecifications (e.g., "male" vs. "M") will calculate as separate groups. 


### Input type

`move2::move2_loc`

### Output type

`move2::move2_loc`

### Artefacts

*Tracking history*: Figure (`tracking_history.png`) detailing the start and end dates of each individual across the study and during the tracking period, along with an indicator of how each individual was terminated (death, censored, survived). 

*Life table:* Output of KM survival analysis; table (`survival_statistics.csv`) with the time, number of individuals at risk, number of events, survival, standard error, and upper/lower 95% confidence intervals. 

*KM survival curve:* Plot (`km_survival_curve.png`), depicting survival over time with median survival time indicated. 

*Cumulative hazard plot:* Plot (`cumulative_hazard_plot.png`), depicting total accumulated risk (expected number of accumulated deaths) over time. 

*Log-rank test:* Output of comparing survival curves between groups; table (`logrank_table_statistics.csv`) with test statistics, degrees of freedom, p-value, and pairwise comparisons. 

*Comparison curves:*: Plot (`km_comparison_curves.png`), depicting survival curves by selected group. 


### Settings 

`Start date` and `End date`: Temporal limits of the study period. If left null (default), the analysis will encompass the entire study period. Useful for defining a study year. Unit: `date`. 

`Fix empty start times` and `Fix empty end times`: Defines how the app handles NA dates (`deploy_on_timestamp` and`deploy_off_timestamp`) at the beginning and end of the study. Options include: 
- Replacing with the first/last recorded timestamp (default)
- Replacing with the current date (for end times only)
- Removing missing data 

`Censor capture-related mortality`: If capture-related mortality is a concern, this setting allows users to define a number of days post-capture to exclude from the overall analysis. Default is no censoring. Unit: `days`. 

`Groups for comparison`: If interested in comparing across groups, this identifies the grouping variable. Default is no group comparisons. Options currently include: 
- Sex
- Lifestage 
- Reproductive condition
- Attachment type 


### Most common errors

Please document and send errors to mspalmer.zool@gmail.com. 


### Null or error handling

See null handling outlined in *Settings*. 

The app will log warnings quantifying the reductions in sample size due to data cleaning and censoring, as well as indicating how many null timestamps were updated or removed.  

The app will log warnings when data includes fewer than 10 mortality events and will terminate if there are no mortality events. 

Errors may arise due to how mortality events are recorded (what metadata flag mortality and whether mortality is associated with a specific end date). If mortality is logged in a column not mentioned in this documentation, it will not be recognized. 



*Please indicate for each setting as well as the input data which behaviour the App is supposed to show in case of errors or NULL values/input. Please also add notes of possible errors that can happen if settings/parameters are improperly set and any other important information that you find the user should be aware of.*

**Setting `radius`:** If no radius AND no duration are given, the input data set is returned with a warning. If no radius is given (NULL), but a duration is defined then a default radius of 1000m = 1km is set. 
