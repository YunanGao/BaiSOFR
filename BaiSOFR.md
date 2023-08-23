The `BaiSOFR` package is a tool for fitting and interpreting Bayesian
scalar-on-functional regression (SOFR) models, with the acronym
`BaiSOFR` representing **Bayesian Adaptive and Interpretable
Scalar-on-function regression**. The methods implemented in this package
are derived from ([Gao and Kowal
([2022](#ref-gao2022bayesian))](https://arxiv.org/abs/2203.00784)). This
document is designed to help new users get started with the main
functionalities of `BaiSOFR` in R.

# Introduction

Scalar-on-function regression (SOFR) is a type of regression where a
scalar response variable is predicted based on one or several functional
predictor(s). An example of SOFR model is as follows:

$$
y\_i =  \mu + {\bf z\_i}'{\bf\alpha} +  \sum\_{j=1}^p\int\_{\mathcal{T}\_i} X\_{j,i}( t)\beta\_j (t)\\ dt +\epsilon\_i, \quad \epsilon\_i \stackrel{iid}{\sim} \mathcal{N}(0, \sigma^2), \quad  i=1,...,n. \tag{1}
$$

In this model, the scalar response *y*<sub>*i*</sub> ∈ ℝ is linked to
*p* functional covariates
{*X*<sub>*j*, *i*</sub> : 𝒯<sub>*i*</sub> → ℝ},  *j* = 1, .., *p* and
other scalar covariates $\bf z\_i \in \mathbb{R}^{z\_p}$. The compact
domains 𝒯<sub>*i*</sub> ⊆ 𝒯 ⊂ ℝ are subject-specific, and 𝒯 denotes the
maximal domain over which the functional covariates *X*<sub>*i*</sub>
may be observed. The integral term
∫<sub>𝒯<sub>*i*</sub></sub>*X*<sub>*j*, *i*</sub>(*t*)*β*<sub>*j*</sub>(*t*) *d**t*
represents the cumulative effect of *X*<sub>*j*, *i*</sub> over its
domain 𝒯<sub>*i*</sub>, as captured by the unknown regression
coefficient function *β*<sub>*j*</sub>. The linear term
$\bf z\_i'\bf \alpha$ accounts for other scalar covariates that may
influence the response *y*<sub>*i*</sub>.

A closely related model to the SOFR model is the distributed lag model
(DLM) ([Dziak et al. 2019](#ref-dziak2019scalar)). In the DLM, the
integral in (1) is replaced with multiple regression terms that utilize
{*X*<sub>*j*, *i*</sub>(*t*)} as covariates. However, unlike the SOFR
model, the DLM requires the functional covariates *X*<sub>*j*, *i*</sub>
to be observed at the same discrete points for all observations. This
requirement could be challenging to meet in practical applications.
Therefore, we prefer the usage of the SOFR model.

This guide will walk you through two main functionalities of the
`BaiSOFR` package using a simulated dataset:

1.  **Bayesian Scalar-on-Function Regression**: One of the key
    objectives in fitting a SOFR model is to estimate and characterize
    the unknown regression coefficient functions *β*<sub>*j*</sub>, as
    these functions illustrate the influence of functional covariates
    *X*<sub>*j*, *i*</sub> on the scalar response *y*<sub>*i*</sub>. The
    `BaiSOFR` package enables users to fit a Bayesian SOFR model to
    their data using the `fitBASOFR()` function. This Bayesian SOFR
    model is specifically designed to **capture both smooth and abrupt
    variations** in the association function {*β*<sub>*j*</sub>}. It
    provides full posterior uncertainty quantification while maintaining
    computational scalability, concerning both the sample size and the
    number of observation points for each functional covariate.

2.  **Critical window selection for Bayesian SOFR Models**: The other
    main goal of fitting a SOFR model is to identify *critical windows*.
    These are regions within the domain 𝒯 where *X*<sub>*j*, *i*</sub>
    are predictive of *y*<sub>*i*</sub>, adjusted for other functional
    covariates and key confounding variables. Due to the
    high-dimensional and highly correlated nature of
    *X*<sub>*j*, *i*</sub>, this task regarding SOFR interpretation is
    not trivial ([Dziak et al. 2019](#ref-dziak2019scalar)). To improve
    interpretatbility, the `BaiSOFR` package offers **decision
    analysis** tools to extract more *interpretable* **locally constant
    point estimates** from **any** given Bayesian SOFR models using the
    `Interpret_SOFR()` function.

We start from installing and loading the `BaiSOFR` package from
[GitHub](https://github.com/) with:

    # devtools::install_github("YunanGao/BaiSOFR")
    library(BaiSOFR)

# Simulate the SOFR Data

In order to illustrate the usage of the `BaiSOFR` package, we will
simulate data according to the SOFR model (1) using the
`simulate_SOFR()` function. This data includes the response variable
*y*<sub>*i*</sub>, the functional covariates
{*X*<sub>1, *i*</sub>, *X*<sub>2, *i*</sub>, *X*<sub>3, *i*</sub>}, and
the scalar covariates **z**<sub>**i**</sub> (comprising 5 continuous
covariates and 2 dummy covariates). We generate this data for
*i* = 1, ..., *n* = 1000 subjects:

    set.seed(1234) # set random seed to ensure duplication of the results

    # simulate the SOFR data
    data = simulate_SOFR(n = 1000, # number of observations
                         p = 3, # number of functional covariates
                         SNR = 5, # Signal-to-noise ratio
                         num_z_continuous = 5, # number of continuous covariates in z
                         num_z_dummy = 2, # number of dummy covariates in z
                         beta_types = c('spiky','smooth','stepwise'), # the ground truth shape of beta
                         domain_X = c(0, 1), # maximal domain over which X are being observed
                         max_obs = 200, # maximal number of discrete points to observe X
                         subject_specific_domain = TRUE # whether to generate subject-specific domain
                         )

    # View the simulated dataset
    names(data)
    #> [1] "t.grid"                  "X"                      
    #> [3] "y"                       "Z_continuous"           
    #> [5] "Z_dummy"                 "groundTruth"            
    #> [7] "subject_specific_domain"

Here:

-   The `SNR` argument stands for signal-to-noise ratio, defined as
    sd(*y*<sub>*i*</sub>−*ϵ*<sub>*i*</sub>)/sd(*y*<sub>*i*</sub>). Lower
    values generally make the modeling task more challenging.

-   The `domain_X` and `max_obs` arguments define 𝒯, and the discrete
    observation points for all *X*. With `domain_X = c(0, 1)` and
    `max_obs = 200` in our case, we have 𝒯 = \[0,1\] and the discrete
    observation points are `seq(0,1,length.out=200)`. Setting
    `subject_specific_domain=TRUE` means that each
    *X*<sub>*j*, *i*</sub> is observed only on its subject-specific
    domain \[0,*T*<sub>*i*</sub>\], where *T*<sub>*i*</sub> ≤ 1.

-   The `beta_types` argument lets you specify the ground truth shape of
    each *β*<sub>*j*</sub> function. Options include `'spiky'`,
    `'smooth'`, and `'stepwise'`. We can visualize the ground truth
    functions using the `plot_beta()` function:

<!-- -->


    # Visualize the ground truth beta functions
    plot_beta(beta_list = data$groundTruth$beta, x.grids = data$t.grid)

![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-4-1.png)

Given that we specified `beta_types = c('spiky','smooth','stepwise')` in
the simulation code, the curves *β*<sub>1</sub>, *β*<sub>2</sub>, and
*β*<sub>3</sub> correspond respectively to `'spiky'`, `'smooth'`, and
`'stepwise'` in this plot.

# Bayesian SOFR Model

Once the data is prepared, we are interested in fitting the SOFR model
to obtain posterior inference for
$\\{\beta\_1,\beta\_2, \beta\_3, \mu, {\bf\alpha}, \sigma^2 | y, X, \bf z\\}$.
It’s important to note that the true regression function’s shape, such
as *β*<sub>1</sub>(*t*) in the plot above, may either vary smoothly or
exhibit abrupt changes. If we fail to accurately capture the shape of
the true regression function, we can end up with biased estimates and
poorly calibrated uncertainty quantification for *β*<sub>*j*</sub>, as
well as suboptimal predictions for *y*. To adapt to both smooth and
abrupt changes in *β*<sub>*j*</sub>, we developed our method. A brief
technical overview of our approach is provided in the [Bayesian SOFR
Model Overview](#BASOFR_overview) section, and its application is
demonstrated in the [Fitting a Bayesian SOFR Model](#BASOFR_fit)
section.

## Bayesian SOFR Model Overview

The Bayesian Scalar-on-Function Regression method serves as a useful
tool for estimating and making inferences for *β*<sub>*j*</sub>.
Specifically designed to adapt to both smooth and rapid changes, the
BASOFR model applies a B-spline basis expansion for each
*β*<sub>*j*</sub> and employs a dynamic horseshoe prior ([Kowal,
Matteson, and Ruppert 2019](#ref-kowal2019dynamic)) on the second
difference of the basis coefficients. This is defined for
*j* = 1, ..., *p* as:

$$
\begin{aligned}
\beta\_j(t) &= \sum\_{k=1}^{K\_B} B\_{j,k}^\*\psi\_k(t)\\\\
\Delta^2  B\_{j,k}^\* |\lambda\_{j,k} &\stackrel{\text{indep}}{\sim} \mathcal{N}(0, \lambda\_{j,k}^2), \quad \\{\lambda\_{j,k} \\}\sim \text{DHS}
\end{aligned}
$$

Here, {*ψ*<sub>*k*</sub>}<sub>*k* = 1</sub><sup>*K*<sub>*B*</sub></sup>
represents a collection of equally-spaced B-spline bases, and DHS refers
to the dynamic horseshoe prior ([Kowal, Matteson, and Ruppert
2019](#ref-kowal2019dynamic)). This *local* and *adaptive* shrinkage
prior is beneficial for function estimation as it encourages smoothness
but can also capture rapid changes in *β*<sub>*j*</sub>. As a result,
the BASOFR method has shown substantially superior performance in
scenarios that features both smooth and abrupt changes, as evidenced in
various simulation studies in Gao and Kowal
([2022](#ref-gao2022bayesian)).

## Fiting a Bayesian SOFR Model

The `fitBASOFR()` function serves as a straightforward tool to implement
the proposed SOFR method. This function only requires the inputs `y`,
`X`, `z`, and `t.grid` from the SOFR model (1). In our example, `t.grid`
corresponds to `seq(0,1,length.out=200)`, which are the points at which
the functions {*X*<sub>*j**i*</sub>}<sub>*j* = 1</sub><sup>*p*</sup> are
observed, and also the points where *β*<sub>*j*</sub> will be estimated.
Since the `fitBASOFR()` function does not perform any scaling or
centering for pre-processing the `z` values, we recommend users to
conduct any necessary pre-processing, such as normalization or
standardization, before passing them into the `z` argument. The function
uses a Gibbs sampler to draw from the posterior distribution:

    # Create the z-variables matrix
    z = cbind(scale(data$Z_continuous), data$Z_dummy)

    # Fit the BASOFR model
    BASOFR_fit <- fitBASOFR(y = data$y, X = data$X, 
                            t.grid = data$t.grid,  # seq(0, 1, length.out=200)
                            z = z,
                            burnin = 1000,
                            S = 1000)
    #> Initializing the MCMC...
    #> Starting the Gibbs sampler...
    #> Gibbs sampler  10 % complete
    #> [1] "2023-08-23 17:45:32 CDT"
    #> Gibbs sampler  20 % complete
    #> [1] "2023-08-23 17:45:35 CDT"
    #> Gibbs sampler  30 % complete
    #> [1] "2023-08-23 17:45:38 CDT"
    #> Gibbs sampler  40 % complete
    #> [1] "2023-08-23 17:45:42 CDT"
    #> Gibbs sampler  50 % complete
    #> [1] "2023-08-23 17:45:45 CDT"
    #> Gibbs sampler  60 % complete
    #> [1] "2023-08-23 17:45:49 CDT"
    #> Gibbs sampler  70 % complete
    #> [1] "2023-08-23 17:45:52 CDT"
    #> Gibbs sampler  80 % complete
    #> [1] "2023-08-23 17:45:55 CDT"
    #> Gibbs sampler  90 % complete
    #> [1] "2023-08-23 17:45:59 CDT"
    #> Gibbs sampler is done!

Upon completion, `fitBASOFR()` returns an object containing the
posterior draws of
$\\{\beta\_1,\beta\_2, \beta\_3, \mu, {\bf\alpha}, \sigma^2 | y, X, \bf z\\}$,
as well as the posterior expectations and credible intervals for
{*β*<sub>1</sub>, *β*<sub>2</sub>, *β*<sub>3</sub>}:

    names(BASOFR_fit)
    #> [1] "beta.post"        "intercept.post"   "Alpha.post"       "sigma.2.post"    
    #> [5] "beta.post.mean"   "beta.95CI"        "beta.50CI"        "t.grid"          
    #> [9] "other_parameters"

To visualize the fit, we can use the `plot_BASOFR_fit()` function:


    # Obtain the ground truth beta
    groundTruth.beta = matrix(unlist(data$groundTruth$beta), byrow=TRUE, nrow=3)

    # Visualize the fit of Bayesian SOFR
    plot_BASOFR_fit(BASOFR_fit = BASOFR_fit,
                    groundTruth = groundTruth.beta # optional
                    )

![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-7-1.png)![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-7-2.png)![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-7-3.png)

As shown from the plots, the proposed Bayesian SOFR model manages to
effectively capture the ground truth, even when faced with the task of
modelling complex shapes!

# Identifying Critical Windows for Bayesian SOFR Models

We employ a **decision analysis strategy** to identify critical windows
and provide more interpretable model summaries for Bayesian SOFR models
through the `Interpret_SOFR()` function.

For *any* Bayesian SOFR model ℳ, the `Interpret_SOFR()` function offers:

1.  Optimal **locally constant** point estimates of the regression
    functions (with uncertainty quantification of their predictive
    accuracy) of different complexity (i.e. different numbers of
    changing points).

2.  The **acceptable family** of the locally constant point estimates,
    which are *near-optimal* concerning predictive accuracy compared to
    the “best” locally constant estimator.

3.  The **simplest acceptable family member** of the locally constant
    estimates, prioritizing simplicity while maintaining predictive
    accuracy.

For demonstration, we use the previously fitted Bayesian Adaptive SOFR
model to illustrate the use of the `Interpret_SOFR()` function. This
function requires only *sampling* from the posterior as well as the data
observation:

    BayesSOFR <- BASOFR_fit # It can be any Bayesian SOFR model, does not limit to our BASOFR model.

    SOFR_Interpret <- Interpret_SOFR(# the posterior samples
                                     beta.post = BayesSOFR$beta.post, 
                                     Alpha.post = BayesSOFR$Alpha.post,
                                     intercept.post = BayesSOFR$intercept.post,
                                     sigma2.post = BayesSOFR$sigma.2.post,
                                     # the data observations
                                     y = data$y, 
                                     X = data$X, 
                                     z = z,
                                     t.grid = data$t.grid)
    #> Start calculating the aggregated trajectories of the functional covariates...Start to solve the 'fit-to-the-fit' to obtain the locally constant estimates...Start of 'fit-of-the-fit' process for beta 1  to obtain locally constant estimates.
    #> End of 'fit-of-the-fit' process for beta 1  to obtain locally constant estimates.
    #> Start of 'fit-of-the-fit' process for beta 2  to obtain locally constant estimates.
    #> End of 'fit-of-the-fit' process for beta 2  to obtain locally constant estimates.
    #> Start of 'fit-of-the-fit' process for beta 3  to obtain locally constant estimates.
    #> End of 'fit-of-the-fit' process for beta 3  to obtain locally constant estimates.
    #> Start to extract the acceptable family... the simplest member in the acceptable family will be plotted automatically.

![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-8-1.png)![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-8-2.png)![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-8-3.png)![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-8-4.png)![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-8-5.png)![](/private/var/folders/ls/vvqf5_wj6gb3rzvg8j_7r3vw0000gn/T/Rtmpe5jQ1S/preview-ef9e68f0ad3.dir/BaiSOFR_files/figure-markdown_strict/unnamed-chunk-8-6.png)

`Interpret_SOFR()` automatically generates two plots for each regression
function *β*<sub>*j*</sub>:

The first plot shows how the predictive performance varies across
locally constant point estimates with different complexity, specifically
relative (% change) to the “best” model (dashed grey vertical line). For
each regression function, while the predictive performance generally
improves as complexity increases, it is clear that several simpler
models are highly competitive. This is observation holds particularly
true when considering predictive uncertainty.

If we wish to select a single locally constant estimate, a compelling
choice would be the simplest member in the acceptable family (solid grey
vertical line). This choice represents the locally constant model with
the simplest shape while still upholding the high standard for
predictive accuracy. Subsequently, this simplest member in the
acceptable family is presented in the second plot for each regression
function.

Upon execution, the `Interpret_SOFR()` function returns a list, which
includes the simplest member of the acceptable family
(`simplest_acc_family`), all the locally constant estimates with
different complexity (`locally_constant`), as well as their
corresponding empirical (`empirical_MSE`) and uncertainty quantification
of mean squared errors (`predictive_MSE`):

    names(SOFR_Interpret)
    #> [1] "simplest_acc_family" "locally_constant"    "empirical_MSE"      
    #> [4] "predictive_MSE"

An interesting fact to note is that the simplest member of acceptable
family often exhibits a stronger resemblance to the ground truth,
especially when the ground truth takes a `stepwise` pattern (see the
results above for *β*<sub>3</sub>). Conversely, *β*<sub>2</sub>, which
exhibits a more gradually varying ground truth (`smooth`), has the
simplest acceptable member encompassing a greater number of changing
points. This occurs as the model strives to represent the smooth ground
truth using locally constant segments.

# References

Dziak, John J, Donna L Coffman, Matthew Reimherr, Justin Petrovich,
Runze Li, Saul Shiffman, and Mariya P Shiyko. 2019. “Scalar-on-Function
Regression for Predicting Distal Outcomes from Intensively Gathered
Longitudinal Data: Interpretability for Applied Scientists.” *Statistics
Surveys* 13: 150.

Gao, Yunan, and Daniel R Kowal. 2022. “Bayesian Adaptive and
Interpretable Functional Regression for Exposure Profiles.” *arXiv
Preprint arXiv:2203.00784*.

Kowal, Daniel R, David S Matteson, and David Ruppert. 2019. “Dynamic
Shrinkage Processes.” *Journal of the Royal Statistical Society Series
B: Statistical Methodology* 81 (4): 781–804.
