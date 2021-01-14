Linkage analysis for H.zea CHR13 marker with larval growth
================
Katherine Taylor and Megan Fritz
written Jan. 8, 2021

## Mapping families on treated and untreated leaf tissue incorporated diet.

One F2 mapping family was generated and progeny split into two groups at 48h after hatching - half were placed on diet with untreated leaf tissue (orange), the other half were placed on diet containing Cry1Ab treated leaf tissue (blue). A second F2 mapping family was generated and progeny were also split into two groups at 48h after hatching - half were placed on diet containing Cry1A.105 + Cry2Ab2 treated leaf tissue (purple). While we put the other half of the family on diet with untreated leaf tissue from a sweet corn isoline with the same genetic background as the two-toxin treated tissue, the data are not shown for simplicity. All larvae were allowed to feed for 7 days and then weighed.

``` r
joined %>%
  filter (genotype_9409b != "") %>%
ggplot(aes(x = genotype_9409b, y = end_weight_mg, fill = TreatByFam, color = TreatByFam )) +
  geom_boxplot(fatten = 1, alpha = 0.8, notch = F) +
  geom_point(position=position_jitterdodge(jitter.width = 0.3)) +
  labs(y = "Weight (mg)", x = "") +
  scale_y_continuous(breaks = c(0, 100, 200, 300)) +
  ylim(0,350) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_x_discrete(breaks  = c("GG", "AG", "AA"),
                   labels = c("A", "A/V", "V")) +
    scale_fill_manual(name="Treatment By Family",
                         breaks=c("Fam 1 CL", "Fam 1 DD", "Fam 2 DD"),
                         values=c("#E69F00", "#56B4E9", "#9932CC")) +
  scale_color_manual(name="Treatment By Family ",
                         breaks=c("Fam 1 CL", "Fam 1 DD", "Fam 2 DD"),
                         values=c("#E69F00", "#56B4E9", "#9932CC")) 
```

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-1-1.png)

## Counts of larvae with each genotype & test for mendelian seg

I tried the chisq.test function with and without simulate.p.value, the latter of which is supposed to be better for smaller sample sizes. It did not change the outcome.

``` r
#test of mendelian segregation

exp_props <- c(0.25, 0.5, 0.25)

geno9409_dist_BCO805_DD <- as.matrix(table(BV_BA52_BZM_P11_A1_DD$genotype_9409b))#allele freqs diag dose BC0805
print(geno9409_dist_BCO805_DD[-1,])#AA is the derived genotype in the field
```

    ## AA AG GG 
    ## 16 36 12

``` r
chisq.test(geno9409_dist_BCO805_DD[-1,], exp_props, simulate.p.value = T)
```

    ## 
    ##  Pearson's Chi-squared test with simulated p-value (based on 2000
    ##  replicates)
    ## 
    ## data:  geno9409_dist_BCO805_DD[-1, ] and exp_props
    ## X-squared = 3, df = NA, p-value = 1

``` r
geno9409_dist_BCO805_CL <- as.matrix(table(BV_BA52_BZM_P11_A1_CL$genotype_9409b))#allele freqs control diet
print(geno9409_dist_BCO805_CL[-1,])#AA is the derived genotype in the field
```

    ## AA AG GG 
    ## 17 23 15

``` r
chisq.test(geno9409_dist_BCO805_CL[-1,], exp_props, simulate.p.value = T)
```

    ## 
    ##  Pearson's Chi-squared test with simulated p-value (based on 2000
    ##  replicates)
    ## 
    ## data:  geno9409_dist_BCO805_CL[-1, ] and exp_props
    ## X-squared = 3, df = NA, p-value = 1

``` r
geno9409_dist_Obs_DD <- as.matrix(table(BV_CV98_03_BZF_I9_DD$genotype_9409b))#allele freqs control diet
print(geno9409_dist_Obs_DD)#AA is the derived genotype in the field
```

    ##    [,1]
    ## AA   15
    ## AG   37
    ## GG   15

``` r
chisq.test(geno9409_dist_Obs_DD, exp_props, simulate.p.value = T)
```

    ## 
    ##  Pearson's Chi-squared test with simulated p-value (based on 2000
    ##  replicates)
    ## 
    ## data:  geno9409_dist_Obs_DD and exp_props
    ## X-squared = 3, df = NA, p-value = 0.3423

## Anova and data transformation for Family 1 - Bv-BA52

We check the assumptions with plots and a shapiro test and then run a boxcox test to determine the appropiate transformation. The lambda value of ~0.5 suggests a square root transformation for individuals on treated diet. Nothing needed for controls.

Progeny on treated diet

``` r
fitDD <- aov(end_weight_mg ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_DD)

shapiro.test(resid(fitDD))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  resid(fitDD)
    ## W = 0.95212, p-value = 0.01457

``` r
plot(fitDD)
```

![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-3-1.png)![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-3-2.png)![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-3-3.png)![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-3-4.png)

``` r
boxcox(fitDD, plotit = TRUE)
```

![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-3-5.png)

Progeny on untreated diet

``` r
fitCL <- aov(end_weight_mg ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_CL)

shapiro.test(resid(fitCL))#no transformation needed.
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  resid(fitCL)
    ## W = 0.97407, p-value = 0.2786

``` r
plot(fitCL)
```

![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-4-1.png)![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-4-2.png)![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-4-3.png)![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-4-4.png)

``` r
boxcox(fitCL, plotit = TRUE)
```

![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-4-5.png)

For progeny on treated diet, trying a non-parametric test because the residuals are not normally distributed.

``` r
kruskal.test(end_weight_mg ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_DD)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  end_weight_mg by genotype_9409b
    ## Kruskal-Wallis chi-squared = 3.586, df = 2, p-value = 0.1665

Then tried a parametric test with sqrt transformed responses.

``` r
fit_BV_BA52_BZM_P11_A1_DD <- aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_DD)
summary(fit_BV_BA52_BZM_P11_A1_DD)
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)  
    ## genotype_9409b  2  22.61   11.30   2.463 0.0936 .
    ## Residuals      61 279.96    4.59                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
shapiro.test(resid(fit_BV_BA52_BZM_P11_A1_DD))#transformation worked
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  resid(fit_BV_BA52_BZM_P11_A1_DD)
    ## W = 0.97051, p-value = 0.1286

Just ran parametric test for individuals raised up on untreated diet.

``` r
fit_BV_BA52_BZM_P11_A1_CL <- aov(end_weight_mg ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_CL)
summary(fit_BV_BA52_BZM_P11_A1_CL)
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)  
    ## genotype_9409b  2  22087   11044   2.432 0.0978 .
    ## Residuals      52 236134    4541                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
shapiro.test(resid(fit_BV_BA52_BZM_P11_A1_CL))#shows no transformation needed.
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  resid(fit_BV_BA52_BZM_P11_A1_CL)
    ## W = 0.97407, p-value = 0.2786

# What about a glm?

Blocking on square because individuals within the same square got diet from the same syringe.

``` r
fit_glmF <- lmer(sqrt(end_weight_mg) ~ 1 + genotype_9409b + (1|square), data = BV_BA52_BZM_P11_A1_DD)
fit_glmR <- lmer(sqrt(end_weight_mg) ~ 1 + (1|square), data = BV_BA52_BZM_P11_A1_DD)

summary(fit_glmF)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: sqrt(end_weight_mg) ~ 1 + genotype_9409b + (1 | square)
    ##    Data: BV_BA52_BZM_P11_A1_DD
    ## 
    ## REML criterion at convergence: 274.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8715 -0.7711 -0.1183  0.8895  1.9493 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  square   (Intercept) 0.1109   0.333   
    ##  Residual             4.4924   2.120   
    ## Number of obs: 64, groups:  square, 7
    ## 
    ## Fixed effects:
    ##                  Estimate Std. Error t value
    ## (Intercept)        8.5444     0.5464  15.638
    ## genotype_9409bAG   0.9909     0.6406   1.547
    ## genotype_9409bGG  -0.3358     0.8137  -0.413
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) g_9409A
    ## gnty_9409AG -0.811        
    ## gnty_9409GG -0.636  0.542

``` r
anova(fit_glmR, fit_glmF)#adding the blocking factor doesn't help.
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: BV_BA52_BZM_P11_A1_DD
    ## Models:
    ## fit_glmR: sqrt(end_weight_mg) ~ 1 + (1 | square)
    ## fit_glmF: sqrt(end_weight_mg) ~ 1 + genotype_9409b + (1 | square)
    ##          Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
    ## fit_glmR  3 286.77 293.25 -140.38   280.77                           
    ## fit_glmF  5 286.06 296.86 -138.03   276.06 4.7036      2     0.0952 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#forcing full model effects without the slope.
fit_glmF2 <- lmer(sqrt(end_weight_mg) ~ 0 + genotype_9409b + (1|square), data = BV_BA52_BZM_P11_A1_DD)
summary(fit_glmF2)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: sqrt(end_weight_mg) ~ 0 + genotype_9409b + (1 | square)
    ##    Data: BV_BA52_BZM_P11_A1_DD
    ## 
    ## REML criterion at convergence: 274.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8715 -0.7711 -0.1183  0.8895  1.9493 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  square   (Intercept) 0.1109   0.333   
    ##  Residual             4.4924   2.120   
    ## Number of obs: 64, groups:  square, 7
    ## 
    ## Fixed effects:
    ##                  Estimate Std. Error t value
    ## genotype_9409bAA   8.5444     0.5464   15.64
    ## genotype_9409bAG   9.5353     0.3762   25.34
    ## genotype_9409bGG   8.2087     0.6284   13.06
    ## 
    ## Correlation of Fixed Effects:
    ##             g_9409AA g_9409AG
    ## gnty_9409AG 0.072            
    ## gnty_9409GG 0.046    0.062

``` r
#now for the controls
fit_glmF <- lmer(end_weight_mg ~ 1 + genotype_9409b + (1|square), data = BV_BA52_BZM_P11_A1_CL)
fit_glmR <- lmer(end_weight_mg ~ 1 + (1|square), data = BV_BA52_BZM_P11_A1_CL)

summary(fit_glmF)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: end_weight_mg ~ 1 + genotype_9409b + (1 | square)
    ##    Data: BV_BA52_BZM_P11_A1_CL
    ## 
    ## REML criterion at convergence: 587.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3194 -0.5340  0.1167  0.4433  1.7867 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  square   (Intercept) 1399     37.40   
    ##  Residual             3400     58.31   
    ## Number of obs: 55, groups:  square, 7
    ## 
    ## Fixed effects:
    ##                  Estimate Std. Error t value
    ## (Intercept)        234.40      20.71  11.319
    ## genotype_9409bAG   -41.78      19.50  -2.143
    ## genotype_9409bGG   -68.84      22.18  -3.103
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) g_9409A
    ## gnty_9409AG -0.551        
    ## gnty_9409GG -0.522  0.534

``` r
anova(fit_glmR, fit_glmF)#but the controls look like there may be statistically significant differences in growth.
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: BV_BA52_BZM_P11_A1_CL
    ## Models:
    ## fit_glmR: end_weight_mg ~ 1 + (1 | square)
    ## fit_glmF: end_weight_mg ~ 1 + genotype_9409b + (1 | square)
    ##          Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
    ## fit_glmR  3 624.76 630.78 -309.38   618.76                           
    ## fit_glmF  5 620.07 630.11 -305.04   610.07 8.6901      2    0.01297 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Excluding small individuals. Some progeny from Family 1 were still alive at the end of the assay but with failure to thrive phenotypes, even in the controls.

Removing small individuals BV\_BA52\_BZM\_P11\_A1\_DD does not change the outcome.

``` r
BV_BA52_BZM_P11_A1_DD_20 <- BV_BA52_BZM_P11_A1_DD %>%
  filter (end_weight_mg > 20)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_DD_20))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)  
    ## genotype_9409b  2  22.61   11.30   2.463 0.0936 .
    ## Residuals      61 279.96    4.59                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
BV_BA52_BZM_P11_A1_DD_25 <- BV_BA52_BZM_P11_A1_DD %>%
  filter (end_weight_mg > 25)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_DD_25))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## genotype_9409b  2  12.75   6.376   1.497  0.232
    ## Residuals      59 251.28   4.259

``` r
BV_BA52_BZM_P11_A1_DD_30 <- BV_BA52_BZM_P11_A1_DD %>%
  filter (end_weight_mg > 30)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_DD_30))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## genotype_9409b  2   5.69   2.844   0.715  0.493
    ## Residuals      57 226.66   3.976

``` r
BV_BA52_BZM_P11_A1_DD_35 <- BV_BA52_BZM_P11_A1_DD %>%
  filter (end_weight_mg > 35)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_DD_35))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## genotype_9409b  2   4.71   2.356   0.656  0.523
    ## Residuals      55 197.64   3.593

``` r
BV_BA52_BZM_P11_A1_DD_40 <- BV_BA52_BZM_P11_A1_DD %>%
  filter (end_weight_mg > 40)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_BA52_BZM_P11_A1_DD_40))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## genotype_9409b  2   2.77   1.386   0.395  0.676
    ## Residuals      54 189.64   3.512

## Anova and data transformation for Family 2 - Bv-CV98

Also check assumptions of normality with plots and a shapiro test and then run a boxcox test to determine the appropiate transformation. The lambda value of ~0.5 suggests a square root transformation for individuals on treated diet. Nothing needed for controls.

Progeny on treated diet

``` r
fitDD <- aov(end_weight_mg ~ genotype_9409b, data = BV_CV98_03_BZF_I9_DD)

shapiro.test(resid(fitDD))#while this is not coming up with a p-value < 0.05, it still seems pretty close to non-normally distributed.
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  resid(fitDD)
    ## W = 0.96808, p-value = 0.08258

``` r
plot(fitDD)
```

![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-10-1.png)![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-10-2.png)![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-10-3.png)![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-10-4.png)

``` r
boxcox(fitDD, plotit = TRUE)
```

![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-10-5.png)

Gonna try a non-parametric test because the residuals questionable.

``` r
kruskal.test(end_weight_mg ~ genotype_9409b, data = BV_CV98_03_BZF_I9_DD)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  end_weight_mg by genotype_9409b
    ## Kruskal-Wallis chi-squared = 2.8318, df = 2, p-value = 0.2427

Then tried a parametric test with sqrt transformed responses.

``` r
fit_BV_CV98_03_BZF_I9_DD <- aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_CV98_03_BZF_I9_DD)
summary(fit_BV_CV98_03_BZF_I9_DD)
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## genotype_9409b  2   15.2   7.599   1.566  0.217
    ## Residuals      64  310.7   4.854

``` r
shapiro.test(resid(fit_BV_CV98_03_BZF_I9_DD))#transformation worked
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  resid(fit_BV_CV98_03_BZF_I9_DD)
    ## W = 0.9688, p-value = 0.09041

# A glm for Family 2?

Blocking on square because individuals within the same square got diet from the same syringe.

``` r
fit_glmF <- lmer(sqrt(end_weight_mg) ~ 1 + genotype_9409b + (1|square), data = BV_CV98_03_BZF_I9_DD)
fit_glmR <- lmer(sqrt(end_weight_mg) ~ 1 + (1|square), data = BV_CV98_03_BZF_I9_DD)

summary(fit_glmF)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: sqrt(end_weight_mg) ~ 1 + genotype_9409b + (1 | square)
    ##    Data: BV_CV98_03_BZF_I9_DD
    ## 
    ## REML criterion at convergence: 291.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9895 -0.7595  0.2044  0.6533  2.0182 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  square   (Intercept) 0.129    0.3591  
    ##  Residual             4.754    2.1803  
    ## Number of obs: 67, groups:  square, 4
    ## 
    ## Fixed effects:
    ##                  Estimate Std. Error t value
    ## (Intercept)        5.8672     0.5918   9.914
    ## genotype_9409bAG   1.0659     0.6688   1.594
    ## genotype_9409bGG   0.1723     0.7976   0.216
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) g_9409A
    ## gnty_9409AG -0.804        
    ## gnty_9409GG -0.673  0.597

``` r
anova(fit_glmR, fit_glmF)#adding the blocking factor doesn't help.
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: BV_CV98_03_BZF_I9_DD
    ## Models:
    ## fit_glmR: sqrt(end_weight_mg) ~ 1 + (1 | square)
    ## fit_glmF: sqrt(end_weight_mg) ~ 1 + genotype_9409b + (1 | square)
    ##          Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## fit_glmR  3 302.12 308.73 -148.06   296.12                         
    ## fit_glmF  5 302.89 313.92 -146.45   292.89 3.2245      2     0.1994

## Excluding small individuals. Some progeny were still alive but with failure to thrive phenotypes.

Removing small individuals only impacts the outcome in one case - removal of progeny with weights of 30mg or less.

``` r
BV_CV98_03_BZF_I9_DD_20 <- BV_CV98_03_BZF_I9_DD %>%
  filter (end_weight_mg > 20)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_CV98_03_BZF_I9_DD_20))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## genotype_9409b  2   3.74   1.870   0.699  0.502
    ## Residuals      51 136.36   2.674

``` r
BV_CV98_03_BZF_I9_DD_25 <- BV_CV98_03_BZF_I9_DD %>%
  filter (end_weight_mg > 25)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_CV98_03_BZF_I9_DD_25))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## genotype_9409b  2   3.45   1.724   0.753  0.477
    ## Residuals      47 107.64   2.290

``` r
BV_CV98_03_BZF_I9_DD_30 <- BV_CV98_03_BZF_I9_DD %>%
  filter (end_weight_mg > 30)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_CV98_03_BZF_I9_DD_30))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)  
    ## genotype_9409b  2  10.95   5.475   3.307 0.0464 *
    ## Residuals      42  69.54   1.656                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
BV_CV98_03_BZF_I9_DD_35 <- BV_CV98_03_BZF_I9_DD %>%
  filter (end_weight_mg > 35)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_CV98_03_BZF_I9_DD_35))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## genotype_9409b  2   5.74   2.871   2.024  0.146
    ## Residuals      38  53.90   1.419

``` r
BV_CV98_03_BZF_I9_DD_40 <- BV_CV98_03_BZF_I9_DD %>%
  filter (end_weight_mg > 40)

summary(aov(sqrt(end_weight_mg) ~ genotype_9409b, data = BV_CV98_03_BZF_I9_DD_40))
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## genotype_9409b  2   5.50   2.752   2.122  0.134
    ## Residuals      36  46.67   1.296

## Final thing. Since filtering out an end weight of 30 or less made a difference for the anova in Family 2 on ObsII, gonna make a plot with those filtered to see what it looks like. The figure clearly shows that this randomly chosen cutoff value produces a trend in the data that is inconsistent with our field study. Furthermore when we removed individuals with weights slightly higher (up to 40 mg or less), the model returns to no statistical significance.

``` r
joined_30 <- joined %>%
  filter (end_weight_mg > 30)

joined_30 %>%
  filter (genotype_9409b != "") %>%
ggplot(aes(x = genotype_9409b, y = end_weight_mg, fill = TreatByFam, color = TreatByFam )) +
  geom_boxplot(fatten = 1, alpha = 0.8, notch = F) +
  geom_point(position=position_jitterdodge(jitter.width = 0.3)) +
  labs(y = "Weight (mg)", x = "") +
  scale_y_continuous(breaks = c(0, 100, 200, 300)) +
  ylim(0,350) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_x_discrete(breaks  = c("GG", "AG", "AA"),
                   labels = c("A", "A/V", "V")) +
  scale_fill_manual(name="Treatment By Family",
                         breaks=c("Fam 1 CL", "Fam 1 DD", "Fam 2 DD"),
                         values=c("#E69F00", "#56B4E9", "#9932CC")) +
  scale_color_manual(name="Treatment By Family ",
                         breaks=c("Fam 1 CL", "Fam 1 DD", "Fam 2 DD"),
                         values=c("#E69F00", "#56B4E9", "#9932CC")) 
```

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

![](Reanalysis_01.11.2020_files/figure-markdown_github/unnamed-chunk-15-1.png)

## Conclusion

All we can really say is that the AA genotype (which produces the V) has significantly faster larval growth on the diet containing untreated leaf tissue. BUT there is a trend that heterozygotes look to growth larger on the diet containing treated leaf tissue for both one and two toxin cultivars.
