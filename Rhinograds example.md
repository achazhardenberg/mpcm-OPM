<h3><span style="text-decoration: underline;"><strong>Sources</strong></span></h3>
  
  <h4>R packages</h4>
		<p style="padding-left: 30px;"><code>"ape"</code> (Paradis et al 2004) </p>
    <p style="padding-left: 30px;"><code>"caper"</code> (Orme et al 2013) </p>
    <p style="padding-left: 30px;"><code>"nlme"</code> (Pinheiro et al 2014) </p>
    <p style="padding-left: 30px;"><code>"geiger"</code> (Harmon et al 2008) </p>

  <h4>Data</h4>
		<p style="padding-left: 30px;"><strong>Rhinograds life history data</strong> (<code>"rhino.csv"</code>), comma separated file of life history traits for 100 species of Rhinogradentia. SP = Code of Species Name; BM = Average Body Mass; NL = Average Nose Length; LS = Average Litter Size; DD = Average Dispersal Distance; RS = Range Size</p>
    <p style="padding-left: 30px;"><strong>Rhinograds phylogenetic tree</strong> (<code>"rhino.tree"</code>), Newick file of the Rhinogradentia phylogeny</p>

<h3><span style="text-decoration: underline;"><strong>Codes</strong></span></h3>

<h4>A Step-by-Step Guide to Phylogenetic Path Analysis Using the d-Sep Method</h4>

To follow the example in Section 8.4, first of all load the needed packages (including all dependencies):


```r
library(ape)
library(caper)
library(nlme)
```

Now, after downloading the ```rhino.csv``` file to your R working directory, create a dataframe from it called ```rhino.dat```:

```r
rhino.dat <- read.csv("rhino.csv")
```

You can look at the structure of the dataframe object you just created:

```r
str(rhino.dat)
```

```
'data.frame':	100 obs. of  6 variables:
 $ SP: Factor w/ 100 levels "s1","s10","s100",..: 1 13 24 35 46 57 68 79 90 2 ...
 $ BM: num  -0.767 -1.01 -1.225 1.331 2.492 ...
 $ NL: num  -2.018 -1.742 -2.469 0.341 1.224 ...
 $ LS: num  -0.0948 -2.1463 -0.6321 -0.1343 0.6082 ...
 $ DD: num  -0.792 1.764 -1.806 0.949 1.624 ...
 $ RS: num  -3.393 -1.047 -2.86 -0.279 0.221 ...
```

Now load the phylogenetic tree:


```r
rhino.tree <- read.tree("rhino.tree")
```

You can also plot it if you like:


```r
plot(rhino.tree, cex = 0.5, no.margin = T)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

The first step in any comparative analysis is to ensure that the species in the dataframe match the species that are present in the phylogenetic tree. This step is not at all trivial as careful comparison of the two lists of species (those in the dataframe and those in the tree) can allow us to detect possible omissions, mismatches or typos in the species names which could lead to problems later on. To compare the species in the tree to the species in the data base you can use this handy script (the code is commented with a hashtag before each comment):


```r
list1 <- as.character(rhino.dat$SP)  #Extracts a vector of species names from the data
list2 <- rhino.tree$tip.label  #Extracts a vector of species names from the tree
list1[which((list1 %in% list2) == FALSE)]  #Species in list1 missing in list2
list2[which((list2 %in% list1) == FALSE)]  #Species in list2 missing in list1
```

In this example, for both of the last two lines of codes you should get in the console the output:    ```character(0)```  
which tells you that there is no species name in one list missing in the other. Otherwise you would have got a vector of the species names present in one list but missing in the other one.

To be able to run the PGLS analyses using the package caper we first need to create a comparative 
object which combines information about the data and the tree in a single object.
Note that for the comparative object the species names in the tree have to be included as an additional column of the data frame, luckily the data base we provided has the species in the column SP.
  

```r
com.dat <- comparative.data(rhino.tree, rhino.dat, SP, vcv = TRUE, vcv.dim = 3, 
    warn.dropped = TRUE)
```

The first model (Model 1 in Fig 8.7) has six conditional indipendencies to be tested. The model is based on 5 vertices (V) and 4 arrows (A). Using this information, with a call to the function ```condNum``` (which we have defined in the "Storks Example" OPM), we can verify if indeed the number of implied conditionial indipendencies is correct:


```r
condNum <- function(V, A) {
    (factorial(V)/(2 * factorial(V - 2))) - A
}
condNum(5, 4)
```

```
[1] 6
```

The 6 d-separation statements which define this causal model are the following (see sections 8.3 and 8.4 for details):  
   
(BM, DD){NL}    
(BM, RS){DD}  
(NL, LS){BM}  
(DD, LS){BM, NL}    
(LS, RS){BM, DD}  
(NL, RS){BM, DD}  
  
Which we can translate into the following PGLS models (Please note that for each model we extract the p-value of the conditional independency we are testing): 


```r
m1.1 <- pgls(DD ~ NL + BM, data = com.dat, lambda = "ML")
m1.1p <- summary(m1.1)$coefficients["BM", 4]

m1.2 <- pgls(RS ~ DD + BM, data = com.dat, lambda = "ML")
m1.2p <- summary(m1.2)$coefficients["BM", 4]

m1.3 <- pgls(LS ~ BM + NL, data = com.dat, lambda = "ML")
m1.3p <- summary(m1.3)$coefficients["NL", 4]

m1.4 <- pgls(LS ~ BM + NL + DD, data = com.dat, lambda = "ML")
m1.4p <- summary(m1.4)$coefficients["DD", 4]

m1.5 <- pgls(RS ~ BM + DD + LS, data = com.dat, lambda = "ML")
m1.5p <- summary(m1.5)$coefficients["LS", 4]

m1.6 <- pgls(RS ~ BM + DD + NL, data = com.dat, lambda = "ML")
m1.6p <- summary(m1.6)$coefficients["NL", 4]
```


We can now use the p-values of each of the conditional independencies we tested to calculate the value of the C-statistic:


```r
C1 <- -2 * (log(m1.1p) + log(m1.2p) + log(m1.3p) + log(m1.4p) + log(m1.5p) + 
    log(m1.6p))
```

And the p-value of the C-statistic:


```r
C1.pval <- 1 - pchisq(C1, 2 * 6)
```

To be able to compare the fit of this path model to all the other models we are testing, we calculate the
value of the CICc (the C-statistic information criterion). We can make our lives easier defining a simple function for this which we will call ```CICc()```:


```r
CICc <- function(C, q, n) {
    C + 2 * q * (n/(n - 1 - q))
}
```

In this function ```C``` is the value of the C-statistics (which, for model 1, we have just assigned to ```C1```), ```q``` is the number of parameters (which can easily be calculated just summing up the number of vertices and arrows in the DAG; in the case of model 1, q is thus equal 9) and ```n``` is the sample size (i.e. the number of species in our dataframe, 100 in this example). Putting these arguments in the just created function we will get the CICc value for model 1, which we assign to ```CICc1```: 


```r
CICc1 <- CICc(C1, 9, 100)
```

The second causal model (Model 2 in Fig. 8.7) has 5 vertices and 5 arrows, it thus requires testing 5 conditional independencies:


```r
condNum(5, 5)
```

```
[1] 5
```

The 5 d-separation statements which define this causal model are the following:  
   
(BM, DD){NL}  
(BM, RS){DD, LS}    
(NL, RS){DD, LS, BM}    
(NL, LS){BM}   
(DD, LS){BM, NL}    
  
As we did for model 1 we translate these conditional indipendence statements into PGLS models and extract the p-values of the conditional indipendencies we are testing: 


```r
m2.1 <- pgls(DD ~ NL + BM, data = com.dat, lambda = "ML")
m2.1p <- summary(m2.1)$coefficients["BM", 4]

m2.2 <- pgls(RS ~ DD + LS + BM, data = com.dat, lambda = "ML")
m2.2p <- summary(m2.2)$coefficients["BM", 4]

m2.3 <- pgls(RS ~ DD + LS + BM + NL, data = com.dat, lambda = "ML")
m2.3p <- summary(m2.3)$coefficients["NL", 4]

m2.4 <- pgls(LS ~ BM + NL, data = com.dat, lambda = "ML")
m2.4p <- summary(m2.4)$coefficients["NL", 4]

m2.5 <- pgls(LS ~ BM + NL + DD, data = com.dat, lambda = "ML")
m2.5p <- summary(m2.5)$coefficients["DD", 4]
```

Now we can calculate the C-statistic and its associated p-value for model 2:

```r
C2 <- -2 * (log(m2.1p) + log(m2.2p) + log(m2.3p) + log(m2.4p) + log(m2.5p))
C2.pval <- 1 - pchisq(C2, 10)
```

And finally the CICc value for this model, which we get using the function defined above:

```r
CICc2 <- CICc(C2, 10, 100)
```


The third model (model 3 in Fig. 8.7) requires testing 6 conditional independencies (you can verify this by yourself now, counting the number of vertices and arrows and using the above defined ```condNum``` function): 
  
(BM, DD){NL}   
(BM, RS){NL}  
(DD, LS){BM, NL}  
(LS, RS){BM, NL}  
(NL, LS){BM}  
(DD, RS){NL}  
   

Note that we have already tested some of these conditional independencies for previous causal models, so there is no need to run these PGLS models once again. Given we created objects containing the p-values of these conditional independencies, we can re-use them below to calculate the value of the C-statistic for this model. We thus only need to perform the following 3 independence tests: 


```r
m3.2 <- pgls(RS ~ NL + BM, data = com.dat, lambda = "ML")
m3.2p <- summary(m3.2)$coefficients["BM", 4]

m3.4 <- pgls(RS ~ BM + NL + LS, data = com.dat, lambda = "ML")
m3.4p <- summary(m3.4)$coefficients["LS", 4]

m3.6 <- pgls(RS ~ NL + DD, data = com.dat, lambda = "ML")
m3.6p <- summary(m3.6)$coefficients["DD", 4]
```

Now we calculate the value of the C-statistic adding the results of tests of conditional independencies we previously used (these are PGLS models m1.1, m1.3 and m1.4 above):


```r
C3 <- -2 * (log(m1.1p) + log(m3.2p) + log(m1.4p) + log(m3.4p) + log(m1.3p) + 
    log(m3.6p))
```

The p-value of the C-statistic for model 3 is:

```r
C3.pval <- 1 - pchisq(C3, 12)
```
And finally the CICc value:


```r
CICc3 <- CICc(C3, 9, 100)
```

Model 4 requires testing 5 conditional independencies:
    
(BM, DD){NL}   
(DD, RS){BM, NL}  
(LS, RS){BM, NL}  
(NL, LS){BM}  
(DD, LS){BM, NL}   
  
2 of which (the second and the third) we have never tested before:


```r
m4.2 <- pgls(RS ~ BM + NL + DD, data = com.dat, lambda = "ML")
m4.2p <- summary(m4.2)$coefficients["DD", 4]

m4.3 <- pgls(RS ~ BM + NL + LS, data = com.dat, lambda = "ML")
m4.3p <- summary(m4.3)$coefficients["LS", 4]
```


Once again we can use results of previous models to calculate the C-statistic:


```r
C4 <- -2 * (log(m1.1p) + log(m4.2p) + log(m4.3p) + log(m1.3p) + log(m1.4p))
C4.pval <- 1 - pchisq(C4, 10)
CICc4 <- CICc(C4, 10, 100)
```

Model 5 has 4 conditional independencies in the basis set: 
  
(BM, DD){NL}    
(LS, RS){BM, NL, DD}  
(NL, LS){BM}    
(DD, LS){BM, NL}  
  
Note that we have already tested 3 of them, so we need to test only the second one:

```r
m5.2 <- pgls(RS ~ BM + NL + DD + LS, data = com.dat, lambda = "ML")
m5.2p <- summary(m5.2)$coefficients["LS", 4]

C5 <- -2 * (log(m1.1p) + log(m5.2p) + log(m1.3p) + log(m1.4p))
C5.pval <- 1 - pchisq(C5, 8)

CICc5 <- CICc(C5, 11, 100)
```

Model 6 has 5 conditional independencies in the basis set:
  
(BM, DD){NL}  
(DD, RS){BM, NL}  
(NL, LS){BM, RS}  
(DD, LS){BM, NL}   
(LS, RS){BM}    
  
Two of these conditional indipendencies (the third and the fifth) have never been tested before: 

```r
m6.3 <- pgls(LS ~ BM + RS + NL, data = com.dat, lambda = "ML")
m6.3p <- summary(m6.3)$coefficients["NL", 4]

m6.5 <- pgls(RS ~ BM + LS, data = com.dat, lambda = "ML")
m6.5p <- summary(m6.5)$coefficients["LS", 4]

C6 <- -2 * (log(m1.1p) + log(m4.2p) + log(m6.3p) + log(m1.4p) + log(m6.5p))
C6.pval <- 1 - pchisq(C6, 10)
CICc6 <- CICc(C6, 10, 100)
```

Model 7 has 4 conditional independencies in the basis set: 
  
(BM, DD){NL}  
(DD, RS){BM, LS, NL}  
(NL, LS){BM}  
(DD, LS){BM, NL}   
  
We need only test one of them (the second) as we have already tested the others: 

```r
m7.2 <- pgls(RS ~ BM + LS + NL + DD, data = com.dat, lambda = "ML")
m7.2p <- summary(m7.2)$coefficients["DD", 4]

C7 <- -2 * (log(m1.1p) + log(m7.2p) + log(m1.3p) + log(m1.4p))
C7.pval <- 1 - pchisq(C7, 8)
CICc7 <- CICc(C7, 11, 100)
```

Model 8 has 6 conditional independencies in the basis set:
  
(BM, RS){$\varnothing$}   
(DD, RS){NL}  
(LS, RS){BM}  
(NL, LS){BM}   
(DD, LS){BM, NL}  
(BM, DD){NL}  

Three of these conditional independencies (the first,the third and the fifth) have never been tested before:

```r
m8.1 <- pgls(RS ~ BM, data = com.dat, lambda = "ML")
m8.1p <- summary(m8.1)$coefficients["BM", 4]

m8.3 <- pgls(RS ~ BM + LS, data = com.dat, lambda = "ML")
m8.3p <- summary(m8.3)$coefficients["LS", 4]

m8.5 <- pgls(LS ~ BM + NL + DD, data = com.dat, lambda = "ML")
m8.5p <- summary(m8.5)$coefficients["DD", 4]

C8 <- -2 * (log(m8.1p) + log(m3.6p) + log(m8.3p) + log(m1.3p) + log(m8.5p) + 
    log(m1.1p))
C8.pval <- 1 - pchisq(C8, 12)
CICc8 <- CICc(C8, 9, 100)
```

Finally, model 9 has 5 conditional independencies in the basis set:
  
(BM,RS){LS}  
(DD, RS){LS, NL}  
(NL, LS){BM}    
(DD, LS){BM, NL}  
(BM, DD){NL}  

The first two of these conditional independencies have never been tested before:

```r
m9.1 <- pgls(RS ~ LS + BM, data = com.dat, lambda = "ML")
m9.1p <- summary(m9.1)$coefficients["BM", 4]

m9.2 <- pgls(RS ~ LS + NL + DD, data = com.dat, lambda = "ML")
m9.2p <- summary(m9.2)$coefficients["DD", 4]

C9 <- -2 * (log(m9.1p) + log(m9.2p) + log(m1.3p) + log(m8.5p) + log(m1.1p))
C9.pval <- 1 - pchisq(C9, 10)
CICc9 <- CICc(C9, 10, 100)
```

You now have all the necessary information to determine which is the causal model which best fits the data. You can reproduce the results in Table 2 and Table 3 of the book chapter, printing the C, C.pval and CICc statistics generated for each causal model above. You can create a summary table of all these statistics as follows (ordered by increasing CICc):  
   
     

```r
model.all <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", 
    "Model 7", "Model 8", "Model 9")
C.all <- c(C1, C2, C3, C4, C5, C6, C7, C8, C9)
C.pval.all <- c(C1.pval, C2.pval, C3.pval, C4.pval, C5.pval, C6.pval, C7.pval, 
    C8.pval, C9.pval)
CICc.all <- c(CICc1, CICc2, CICc3, CICc4, CICc5, CICc6, CICc7, CICc8, CICc9)
C.all <- round(C.all, digits = 3)
C.pval.all <- signif(C.pval.all, digits = 3)
CICc.all <- round(CICc.all, digits = 3)
results.dat <- data.frame(model.all, C.all, C.pval.all, CICc.all)
names(results.dat) <- c("Model", "C statistic", "p-value", "CICc")
results.dat <- results.dat[order(results.dat$CICc), ]
print(results.dat)
```

```
    Model C statistic  p-value  CICc
8 Model 8       7.700 8.08e-01 27.70
6 Model 6       6.439 7.77e-01 28.91
4 Model 4       6.582 7.64e-01 29.05
9 Model 9       7.362 6.91e-01 29.83
5 Model 5       5.258 7.30e-01 30.26
7 Model 7       6.018 6.45e-01 31.02
3 Model 3      28.973 3.98e-03 48.97
1 Model 1      63.809 4.52e-09 83.81
2 Model 2      62.769 1.08e-09 85.24
```
  
  
   
Although as mentioned in the book chapter, the data are simulated and drawn randomly from a multivariate normal distribution with mean 0 and standard deviation of 1, we nonetheless provide example of code to standardize data in R for it may be of use to readers when working on their own data:


```r
rhino.scaled <- apply(rhino.dat[, 2:6], 2, scale)
```

We now need to add the species names again:


```r
SP <- rhino.dat$SP
rhino.scaled <- cbind(SP, as.data.frame(rhino.scaled))
```

To use the standardized data in new PGLS models, we need to create a new comparative.data file, combining this data file with the tree: 


```r
scaled.comp <- comparative.data(rhino.tree, rhino.scaled, SP, vcv = TRUE, vcv.dim = 3)
```

We will use the standardized data to obtain the slope and standard error of the following causal links of the model:

```r
mS8.1 <- pgls(LS ~ BM, data = scaled.comp, lambda = "ML")

mS8.2 <- pgls(NL ~ RS + BM, data = scaled.comp, lambda = "ML")
# Note that here nose length has two causal links: from body mass and from
# range size

mS8.3 <- pgls(DD ~ NL, data = scaled.comp, lambda = "ML")
```
  
  
<h4>Simulation of tree and data for the Rhinograds example</h4>

In this tutorial we will show how we simulated the Rhinogradentia life history data and the related phylogenetic tree. You can use this code to generate your own data with an underlying phylogenetic signal to play around and test PPA and other phylogenetic comparative methods. Please be aware that if you generate new rhino.csv and a new rhino.tree files and then use them instead of the provided files to follow the above Rhinograds example, you will get different estimates for statistics and parameters as those reported in the book chapter, as the underlying data and tree will be different. Although, you still should get qualitativelly similar results, as the underlying correlation matrix among the variables and Pagel's lambda are fixed. You can of course also play with these parameters and change completely the underlying "true" causal structure (changing the correlation matrix) and the degree of phylogenetic signal (changing Pagel's lambda value).  

You will need  the packages: ```ape``` and  ```geiger``` (including all dependencies):


```r
library(ape)
library(geiger)
```

Now let's build the correlation matrix among the five variables according to the true causal structure implied by model 8 (Fig. 8.7 in the book). For simplicity, we fix the correlation coefficient of all direct causal links to 0.5. The correlation coefficient of two variables not directly linked, but for which there is a causal pathway passing through a third variable, will be given by the product of the correlation coefficients of the directly linked variables. For example, for the indirect causal link between BM and DD, which passes through NL (BM->NL->DD), the correlation coefficient will be 0.5x0.5 = 0.25. Following the same reasoning, the correlation coefficient between LS and DD, will be given by the product of the correlations between BM-LS, BM-NL and NL-DD and thus 0.5x0.5x0.5 = 0.125. Variables without any causal link between them have a correlation of 0. The order of the variables in the matrix is the following: BM, NL, LS, DD, RS:


```r
mat <- matrix(c(1, 0.5, 0.5, 0.25, 0, 0.5, 1, 0.25, 0.5, 0.5, 0.5, 0.25, 1, 
    0.125, 0, 0.25, 0.5, 0.125, 1, 0.25, 0, 0.5, 0, 0.25, 1), 5, 5)
rownames(mat) <- c("BM", "NL", "LS", "DD", "RS")
colnames(mat) <- c("BM", "NL", "LS", "DD", "RS")
print(mat)
```

```
     BM   NL    LS    DD   RS
BM 1.00 0.50 0.500 0.250 0.00
NL 0.50 1.00 0.250 0.500 0.50
LS 0.50 0.25 1.000 0.125 0.00
DD 0.25 0.50 0.125 1.000 0.25
RS 0.00 0.50 0.000 0.250 1.00
```

Now we simulate a phylogenetic tree for 100 species. We can plot this tree to see what it looks like, and then we save it in the Newick format for future use (it will be saved in your working directory). We give it a different name (```rhino2.tree```), than the phylogenetic tree file used for the example in the book chapter (```rhino.tree```), to avoid confusion:



```r
tree <- drop.tip(sim.bdtree(b = 1, d = 0, stop = "taxa", n = 101), "s101")
plot(tree, cex = 0.5, no.margin = T)
```

![plot of chunk unnamed-chunk-36](figure/unnamed-chunk-36.png) 

```r
write.tree(tree, file = "rhino2.tree")
```

And now we finally simulate the data providing the phylogenetic tree we just created, the correlation matrix (```mat```) defining the underlying causal structure and fixing Pagel's lambda to 0.8. We save this data file as csv file in your working directory calling it ```rhino2.csv```


```r
data <- data.frame(sim.char(rescale(tree, model = "lambda", 0.8), mat, nsim = 1, 
    model = "BM"))
names(data) <- c("BM", "NL", "LS", "DD", "RS")
SP <- row.names(data)
data <- cbind(SP, data)
write.csv(data, file = "rhino2.csv")
```

<h3><span style="text-decoration: underline;"><strong>References</strong></span></h3>
  <ul>
  <li>Harmon L.J., Weir J.T., Brock C.D., Glor R.D. and   Challenger W. (2008) geiger: investigating evolutionary radiations. Bioinformatics 24:129-131. R package version 2.0.3</li>
  <li>Orme D., Freckleton R., Thomas G., Petzoldt T., Fritz S., Isaac N. and Pearse W. (2013). caper: Comparative Analyses of Phylogenetics and Evolution in R. R package version 0.5.2.</li>
  <li>Paradis E., Claude J. and Strimmer K. 2004. ape:
  analyses of phylogenetics and evolution in R
  language. Bioinformatics 20: 289-290.R package version 3.1-4.</li>
  <li>Pinheiro J., Bates D., DebRoy S., Sarkar D. and R Core Team (2014). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-117.</li>
	</ul>
  


