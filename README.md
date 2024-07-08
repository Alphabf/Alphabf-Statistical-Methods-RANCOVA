# Statistical-Methods-RANCOVA
SAS code for repeated measure analysis of covariance (RANCOVA).
The parameter will be analyzed by the analysis phase, with a repeated measure analysis of
covariance (RANCOVA). 
The model will consider baseline, period, animal, treatment, time after dose, and the interaction between treatment group and time after dose.

The analysis uses PROC MIXED in SAS®, with the covariance structure chosen based on AIC comparison.

Dunnett's test will compare treatment levels with the control.

If TRT*TIME interaction is significant, Dunnett's test will be conducted per time interval; otherwise, it will be done across pooled time intervals.
