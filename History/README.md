# MSc Thesis project

## Last week's summary
- [x] Fomulated JMs for one longitudinal and time-to-event data as Latent Gaussian Models (LGM).
- [x] Set a simulation study for JMs, in which the baseline hazard of survival part is parametirc model.
- [x] Set a simulation study for JMs with shared random effects association structure.
- [x] Set a simulation study for JMs with multiple longitudinal biomarkers.

## Summary of the weekly meeting

| Timeline | Description | Tasks & Feedback & suggections |
| --- | --- | --- |
| 6 Dec 2022 | First meeting with Dr. Bardia and supervisors | Meeting summary, understand research questions |
| 15 Dec 2022 | Weekly meeting | Read papers, prepare proposal |
| 22 Dec 2022 | Weekly meeting | Read papers, finish proposal, review LMM and survival analysis |
| 12 Jan 2023 | Weekly meeting | Access to LUMC account, try LMM and Cox model with INLA on two case studies, revised proposal based on feedback |
| 19 Jan 2023 | Weekly meeting | Specify prior, check label information in SPSS, apply LMM, Cox model to the real data set, initial data analysis |
| 26 Jan 2023 | Weekly meeting | Find patients whose AR grade dropped, find outliers, Access the HPC environment, and Try to add more random effects to the LMM |
| 2 Feb 2023 | Weekly meeting | Residual analysis (focus on the gradient), Cox regression (Include time-dependent covariates) |
| 7 Feb 2023 | Weekly meeting | Prepare for client meeting: defining the data |
| 10 Feb 2023 | Client meeting | Receive updated data, reconstruct the dataset, exploratory data analysis |
| 17 Feb 2023 | Weekly meeting | Reconstruct the dataset in R, Exploratory data analysis, Exploring joint modeling with INLA (Time-dependent cox regression?) |
| 24 Feb 2023 | Weekly meeting | Currently parametric model for the survival part makes too strong an assumption on baseline hazard; try cox regression in joint model |
| 2 Mar 2023 | Weekly meeting | Try the 'copy' feature in R-INLA to build a joint model, start with a joint of two models with normal outcomes |
| 9 Mar 2023 | Weekly meeting | Try joint of two models with normal outcomes, unequal number of observations, compare the coefficients from different methods |
| 16 Mar 2023 | Weekly meeting | Only include quantitive variables, then compare the coefficients. Thesis writing |
| 23 Mar 2023 | Weekly meeting | Concern about the baseine hazard that we used, try with parametric assumption on the baseline hazard, Thesis writing |
| 30 Mar 2023 | Weekly meeting | Structure of Thesis, How the INLA calculate conditional joint likelihood, check the estimates of coefficients and their sd |
| 13 April 2023 | Weekly meeting | How the INLA calculates conditional joint likelihood, Possible to use change value in a small time interval to approximate slope |
| 20 April 2023 | Weekly meeting | How the INLA calculates conditional joint likelihood, discussion about Thesis defence schedule |
| 26 April 2023 | Weekly meeting | Discussion about contents in Thesis, and Thesis defense schedule |

## History (Finished tasks before the meeting)

26 April 2023
- [x] set a simulation study about JMs with current slope structure
- [x] set a simulation study about JMs with cumulative effects structure

20 April 2023
- [x] Created slides to show how to build JMs using INLA
- [x] Developed JMs with current slope, cumulative effect association structure
- [x] Simulate a dataset for joint modeling with current value association structure
- [ ] Tried to find papers or codes that showed INLA can compute conditional joint likelihood

13 April 2023
- [x] Finished joint model using INLA for one longitudinal biomarker and time-to-event process
- [x] Compared our results with that from JM, JMbayes
- [x] Tried to find papers or codes that showed INLA can compute conditional joint likelihood
- [x] Tired to developed other association structures (current slope, cumulative effect)

30 Mar 2023
- [x] Debug my current solution to build joint model for longitudinal and survival data
- [x] Construct a new dataset that necessary information is included 
- [x] Exploring the concern about baseline hazard
- [x] Thesis writing (Introduction + part of Method)

23 Mar 2023
- [x] Coefficients comparison (only quantitive variables)
- [x] Coefficients of fixed effects computed by INLA are the same as JM, etc.
- [x] Covariance matrix of random effects computed by INLA are the same as JM, etc.
- [ ] Association (scaling factor) is different? Error in computing method?

16 Mar 2023
- [x] Copied coefficients of longitudinal submodel (function of M(t))
- [x] Using part of the measurement points to split the follow-up time
- [x] Using mid-points as the t in M(t) for survival submodel
- [x] Joint models with current value association structure
- [ ] Coefficients in the survival submodel differ from that computed by JM, etc.
- [ ] It could be caused by factor expanding strategy in INLA

9 Mar 2023
- [x] Learned the 'copy' feature in R-INLA
- [x] Built a joint model of longitudinal measurements and time-to-event process (the standard way in INLA to convert Cox to Poisson)
- [x] Copied random effects in longitudinal part to survival part
- [ ] Copied fixed effects in the longitudinal part while failing to incorporate them into the survival part.

2 Mar 2023
- [x] Fined joint models with INLA (Semi-parametric model for survival part)
- [x] Inference based on time-dependent cox model
- [x] Inference based on linear mixed model

24 Feb 2023
- [x] Built several joint models with INLA (only consider reoperation as the event);
- [x] Built several joint models with INLAjoint (reoperation and death, regarded as competing events)
- [x] Built several joint models with JMbayes2.
- [x] Summarise problems during data reconstruction.

17 Feb 2023
- [x] Reconstruct the dataset in R
- [x] Exploratory data analysis
- [x] Exploring joint modeling with INLA (Time-dependent cox regression?)

7 Feb 2023
- [x] Access the Research LUMC platform.
- [x] LMM: a residual analysis based on the full model (including two-way interactions).
- [x] Time-dependent cox regression (survival package).
- [ ] Time-dependent cox regression (INLA)

2 Feb 2023
- [x] Find patients whose AR grade dropped during follow-up.
- [x] Find outliers of each variable.
- [x] Fix errors when using 'xyplot'.
- [x] Residual analysis for LMM.
- [x] Build summary table on GitHub.
- [ ] Tried joint modeling in the JM package. - Failed due to NAs.
- [ ] Tried to build a LMM with ideal marginal residuals. - Linear trend in residuals at this moment.

25 Jan 2023
- [x] Investigate how to specify the prior parameters for the iidkd model.
- [x] Applied LMM, Survival analysis with INLA to the real data set.
- [x] Initial data analysis: check the frequency of each AR grade transition; check label information in SPSS.
- [ ] Access LUMC computing environment (R and R studio)

17 Jan 2023
- [x] Thesis proposal v2.0
