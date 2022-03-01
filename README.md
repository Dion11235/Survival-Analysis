# Survival Analysis of Patients with Recurrent Bladder Cancer
In this project, we have tried to estimate the survival and hazard functions of 118 patients with bladder cancer, who have entered a study where they received 3 different kinds of treatments :
- placebo
- pyridoxine
- thiotepa

Further down the study, some of the patients have observed multiple recurrences (upto 9), some have died, and the rest have been censored. In order to deal with this data, we have divided the entire dataset into 2 parts :
- 118 patients with only their first event/censoring times (independent)
- the same 118 patients with all recurrence times (dependent)

## Basic EDA :

<p align="center">
  <img width="200" height="200"src="/plots/pie.png">
  <img width="300" height="200"src="/plots/bar1.png">
  <img width="250" height="200"src="/plots/box1.png"><br>
</p>
<p align="center">
  <img width="400" height="300"src="/plots/bar2.png"> <br>
</p>

Since some of the patients have suffered multiple recurrences, the original data has been divided in 2 parts to deal with single recurrences and multiple recurrences separately :
- 118 patients with only their first event/censoring times (independent time-to-event)
- the same 118 patients with all recurrence times (dependent time-to-event)


## Modelling :
### First Group Analysis :
#### Kaplan-Meier Estimate :
The Kaplan–Meier estimator is a non-parametric statistic used to estimate the survival function from time-to-event data. The following is the Kaplan-Meier estimator of the patients grouped by the 3 treatment methods mentioned before. From the graph, we can conclude that,
1. the group on pyridoxine has the highest median survival time among all patients.
2. patients on pyridoxine have a median survival time approximately 4 times greater than the patients on placebo.

<p align="center">
  <img width="700" height="450"src="/plots/KM1.png"> <br>
</p>


The pairwise log-rank test between these 3 groups yield the following results:

| Log-Rank test between treatments | Test Statistic | p-value |
|----------|---------|-------|
| Placebo & Pyridoxine | 1.326279 | 0.249468 | 
| Placebo & Thiotepa | 1.520945 | 0.217477 |
| Pyridoxine & Thiotepa  | 0.001407 | 0.970077 |

The high p-value from the log-rank test between pyridoxine and thiotepa is indicative of the fact that their survival curves are not quite dissimilar, which is also clear from the graph above. The low p-value between placebo and the other 2 groups supports the fact that their survival curves deviate completely in later time points.

Next we have thresholded on the number and size of tumours respectively, on the threshold for which the pairwise log-rank test gives the lowest p-value, which indicates that the two groups are the most dissimilar in terms of survival and hazard functions. The following are the Kaplan-Meier estimators of the patients grouped by the thresholded number and size of tumours respectively. From the graphs below, we can conclude that 
1. patients with less than 4 tumours have a median survival time approximately 7 times greater than those with greater than or equal to 4 tumours.
2. patients whose size of the largest initial tumour was less than 6 cm have a median survival time approximately 4 times greater than the other group.

<p align="left">
  <img width="400" height="250"src="/plots/KM2.png">
  <img width="400" height="250"src="/plots/KM3.png"> <br>
</p>

#### Cox-Proportional Hazard Model :
In the previous section we have worked with the Kaplan-Meier estimator which does not take any covariates into consideration. In this section, we shall be modelling the survival functions by taking into account the 2 covariates - number and size of tumour into the model itself.

**Cox proportional hazards assumptions:**
Before modelling, we must check whether the covariates follow the Cox proportional hazard assumptions - that the hazards are proportional. This is achieved with the Schoenfeld residuals test, which plots the residual values of the covariates against time to look for time dependency.

<p align="center">
  <img width="550" height="400"src="/plots/coxph_assumption.png"> <br>
</p>

As we can clearly see, the covariate residuals do not exhibit any pattern, and are completely random. Also, the p-values for both the covariates are greater than 0.05, hence both the covariates satisfy the proportional hazard assumption, and we can go ahead with modelling. 


### Second Group Analysis :
Unlike the first group, the second group has got multiple recurrences for individual patients.

<p align="center">
  <img width="600" height="400" src="https://github.com/Dion11235/Survival-Analysis/blob/main/plots/group2_snapshot.png"><br>
</p>

Cox proportional Hazard models are unable to capture the dependencies in between the data points (Here in between two recurrence data which comes from the same patient).

So, the Shared Frailty model has been used to model the hazard function of different treatment groups here. It has been assumed that the random variable Z follows Gamma Distribution with finite variance and mean 1, as we have used the multiplicative frailty modeling technique.

<p align="center">
  <img width="500" height="400" src="https://github.com/Dion11235/Survival-Analysis/blob/main/plots/gamma_frailty_summary.png"><br>
</p>

Here we can see the summary of the fitted gamma shared frailty model. The first part i.e. Regression Coefficients ,  gives us :

- Unit increase in the number of tumors and the maximum size of the tumors in the patient with all other covariates unchanged, increases the hazard 1.3745 times and 1.0456 times respectively.
- A patient with a constant number of tumors with constant maximum size is 1.0965 times as likely to die in Pyridoxine treatment as the same in Placebo treatment.
- A patient with a constant number of tumors with constant maximum size is 0.5115 times as likely to die in Thiotepa treatment as the same in Placebo treatment.


In the Commenges-Andersen Test, we check for heterogeneity in the data. And due to a very low p-value we can reject the null hypothesis that is There is no heterogeneity in the data i.e. the frailty variable Z does not follow a mixed distribution.

![This is an image](https://github.com/Dion11235/Survival-Analysis/blob/main/plots/frailty_wrt_trt_both.png)

These are the conditional and marginal survival curves (i.e. the expected survival over Z), predicted by the fitted gamma shared frailty model,  for an individual with one tumor of size 1 cm with respect to all three different treatments. We can instantly conclude :

- Median Life Time for the patient in Placebo treatment : 20 months (approx)
- Median Life Time for the patient in Pyridoxine treatment : 18 months (approx)
- Median Life Time for the patient in Thiotepa treatment : 40 months (approx)

These values correspond to the exponential of coefficients quite clearly, and this gives us the conclusion that, *“The treatment Thiotepa is almost 2 times more effective for Bladder cancer than Pyridoxine or Placebo.”* 

With that the analysis ends here. Thank you for you patience and valuable time.

