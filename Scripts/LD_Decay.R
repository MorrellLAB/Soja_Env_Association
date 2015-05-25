#   A script to estimate LD Decay (D' across distance) in the method proposed
#   by Abecasis et al. 2001.
#   Written by Thomas JY Kono

#   The equation put forth by Abecasis et al. 2001 is
#       E(D_t) = D_low + (D_high - D_low)(1 - [theta])^t * D_0
#   Where E(D_t) is the expected value of D'
#         D_low is the minimum expected value of D'
#         D_high is the maximum expected value of D'
#         D_0 is the initial value of D' (taken to be 1, usually)
#         [theta] is the recombination fraction between the two loci
#         t is the number of generations since D'=1
#   In the paper, [theta] of 1cM/Mb was used, but we have better data on this.
#   We will take the average cM/Mb rate in the pericentromeric and euchromatic
#   regions and use those as scaling factors to estimate the recombination
#   fraction between markers from the physical distance.
#   The values of D_low, D_high and t were chosen to minimize the RMSE of the
#   data.

#   Read in the approrpriate files. Col 1 is distance, and Col 2 is LD.
euch <- read.table("Euchromatin_LDDecay.txt.bz2", header=TRUE)
peri <- read.table("Pericentromere_LDDecay.txt.bz2", header=TRUE)
#   The $ doesn't seem to play well with nls()
peri_l <- peri$LD
peri_d <- peri$Distance
euch_l <- euch$LD
euch_d <- euch$Distance

#   Next, we have to set the starting guesses for our parameters to minimize
#   the RMSE. These can't be exactly 0, since that would produce infinity.
#   nls() expects them to be in a list
start <- list(
    dlow=0.05,
    dhigh=0.95,
    gens=5000)
#   Set the lower bounds. These should be in the same order as the starting
#   values.
lower <- c(0.001, 0.001, 1)
upper <- c(0.999, 0.999, 100000)
#   Then, fit the least-squares curve
#   Scale the physical distance to genetic distance. The physical distances are
#   given in bp, so we divide by 1,000,000 to get Mb. In this case, the
#   median cM/Mb in pericentromeric regions is 2.390736 and in
#   euchromatic regions it is 3.589077. cM is the % recombination between two
#   loci, so we divide by 100 again. This equation has some trouble with values
#   for "recombination fraction" greater than 1, so we scale it by 1/10 again.
#   We will undo this scaling later.
ld_peri_model <- nls(
    peri_l ~ dlow + (dhigh - dlow) * (1 - (peri_d / 1000000000) * 2.390736)^gens,
    start=start,
    control=nls.control(
        maxiter=100,
        warnOnly=T),
    lower=lower,
    upper=upper,
    algorithm="port",
    trace=TRUE)

ld_euch_model <- nls(
    euch_l ~ dlow + (dhigh - dlow) * (1 - (euch_d / 1000000000) * 3.589077)^gens,
    start=start,
    control=nls.control(
        maxiter=100,
        warnOnly=T),
    lower=lower,
    upper=upper,
    algorithm="port",
    trace=TRUE)

#   Now that we fit the model, we can run predict() on it.
#   Set up the vector of distances to produce an estimated D' value
xvals <- seq(0, 500000, 500)*10
peri_yvals <- predict(ld_peri_model, list(peri_d=xvals))
euch_yvals <- predict(ld_euch_model, list(euch_d=xvals))

#   Print out the model summaries, save the parameter estimates
peri_summary <- summary(ld_peri_model)
peri_dlow <- peri_summary$parameters["dlow", "Estimate"]
peri_dhigh <- peri_summary$parameters["dhigh", "Estimate"]
peri_gens <- peri_summary$parameters["gens", "Estimate"]
# Formula: peri_l ~ dlow + (dhigh - dlow) * (1 - (peri_d/1e+09) * 2.390736)^gens

# Parameters:
#        Estimate Std. Error t value Pr(>|t|)    
# dlow  2.431e-01  2.549e-04   953.9   <2e-16 ***
# dhigh 6.758e-01  2.041e-03   331.0   <2e-16 ***
# gens  1.268e+03  9.789e+00   129.5   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.2513 on 1139790 degrees of freedom

# Algorithm "port", convergence message: relative convergence (4)
euch_summary <- summary(ld_euch_model)
euch_dlow <- euch_summary$parameters["dlow", "Estimate"]
euch_dhigh <- euch_summary$parameters["dhigh", "Estimate"]
euch_gens <- euch_summary$parameters["gens", "Estimate"]
# Formula: euch_l ~ dlow + (dhigh - dlow) * (1 - (euch_d/1e+09) * 3.589077)^gens

# Parameters:
#        Estimate Std. Error t value Pr(>|t|)    
# dlow  1.761e-01  4.184e-05  4207.7   <2e-16 ***
# dhigh 7.701e-01  1.282e-03   600.7   <2e-16 ***
# gens  8.550e+03  2.819e+01   303.3   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.1777 on 18279904 degrees of freedom

# Algorithm "port", convergence message: relative convergence (4)

#   Use these coefficients to calculate the half life
#   This is just algebra!
peri_half <- peri_yvals[1]/2
peri_half_scaled <- 1 - 10^(log10((peri_half-peri_dlow)/(peri_dhigh - peri_dlow))/peri_gens)
peri_half_bp <- (peri_half_scaled * 1000000000) / 2.390736
print(peri_half_bp)
euch_half <- euch_yvals[1]/2
euch_half_scaled <- 1 - 10^(log10((euch_half-euch_dlow)/(euch_dhigh - euch_dlow))/euch_gens)
euch_half_bp <- (euch_half_scaled * 1000000000) / 3.589077
print(euch_half_bp)

#   And plot it!
pdf(file="LD_Decay_Plot.pdf", 8, 8)
plot(
    c(0, 1000000),
    c(0, 1),
    type="n",
    main="LD Decay (D')",
    xlab="Physical Distance (bp)",
    ylab="D'"
    )
lines(
    xvals,
    peri_yvals,
    lwd=2,
    col="red",
    lty=1
    )
lines(
    xvals,
    euch_yvals,
    lwd=2,
    col="blue",
    lty=2
    )
abline(
    v=peri_half_bp,
    lwd=0.5,
    col="red",
    lty=1
    )
abline(
    v=euch_half_bp,
    lwd=0.5,
    col="blue",
    lty=2
    )
legend(
    "topright",
    c("Pericentromere", "Euchromatin"),
    col=c("red", "blue"),
    lty=c(1, 2)
    )
dev.off()
