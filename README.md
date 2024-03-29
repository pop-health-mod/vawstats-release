# vawstats
This repository contains the code used to produce estimates of global, regional, and national violence against women (VAW) statistics is presented here. The methods are described in the preprint titled: [*A framework to model global, regional, and national estimates of intimate partner violence*](https://doi.org/10.1101/2020.11.19.20235101).

(The database has yet been publicly released, however. The clearance process from co-custodians is in progress.)

The first two scripts are used to estimate adjustment factors for the crosswalk for the meta-regression of *lifetime intimate partner violence* (IPV) and *past year IPV*:

`1.0a CrossWalk Ever IPV.R`   

`1.0b CrossWalk Past Year IPV.R`

Once the adjustments factors are estimated, we format the *lifetime IPV* and *past year IPV* databases and calcualte desriptive statistics in the two homonymous scripts below.

`2.0a AgeStn JAGS Ever IPV.R`. 

`2.0b AgeStn JAGS Past Year IPV.R`

Finally, the joint Bayesian meta-regression model can be run, along with in-sample comparisions and out-of-samples predictions.

`2.0c AgeStn JAGS Combined IPV.R`


# Collaborators
Mathieu Maheu-Giroux (McGill University)

Lynnmarie Sardinha (University of Bristol)

Heidi Stöckl (London School of Hygiene & Tropical Medicine)

Sarah Meyer (World Health Organization)

Arnaud Godin (McGill)

Monica Alexander (University of Toronto)

Claudia García-Moreno (World Health Organization)
