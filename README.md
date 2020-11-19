# vawstats
The code use to estimate global, regional, and national violence against women statistics is presented here. (The database has yet been publicly released, however.)

The first two scripts are used to estimate adjustment factors for the crosswalk:
`1.0a CrossWalk Ever IPV.R` 
`1.0b CrossWalk Past Year IPV.R`

Once the adjustments factors are estimated, we format the *lifetime IPV* and *past year IPV* database and calcualte desriptive statistics in the two homonymous scripts below.
`2.0a AgeStn JAGS Ever IPV.R`
`2.0b AgeStn JAGS Past Year IPV.R`

Finally, the joint Bayesian meta-regression model can be run, along with in-sample validationsa and out-of-samples comparisons.
`2.0c AgeStn JAGS Combined IPV.R`
