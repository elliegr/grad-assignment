Blog Post: Violating the Probabilistic Assumption
================
Ellie Grace Moore
4/25/2022

      Hey causal cats! Welcome back to my blog. Today we will choose a
specific causal assumption and then generate data that *violates* this
assumption. As we know, the assumptions made in causal inference are key
for the entire process of determining a causal estimate! Thus, in order
to understand these assumptions better, we will see what happens during
the causal process and the final effect violating this assumption has on
our causal estimand.

## The Assumption

      Although there are multiple assumptions in causal inference, today
we will look further into the **probabalistic assumption**. In words,
this means that all units must have non-zero probability for all
exposures. This is also commonly referred to the positivity assumption.
Practically thinking, if we had observations that had a probability of 0
of being assigned to a given exposure, then the causal study is rigged
from the beginning! It is almost as if we are picking and choosing where
we want certain observations to go, thus straying us further away from
the “ideal” randomized trial environment.

### What it Means Mathematically

      Mathematically speaking, the exposure assignment probabilities are
denoted *P*<sub>*i*</sub>(*X*\|*Z*, *Y*(0), *Y*(1)). For the
probabilistic assumption to be met,
*P*<sub>*i*</sub>(*X*\|*Z*, *Y*(0), *Y*(1)) ∈ (0, 1). However, here we
will assign the probability for every fourth observation (arbitrarily
picked) to be zero. So in our case,
*P*<sub>*i*</sub>(*X*\|*Z*, *Y*(0), *Y*(1)) ∈ \[0, 1).

### An Example

       A brief example of a scenario in which the probabilistic
assumption is violated is if we had a study involving both males and
females, but the exposure would only pertain to females. Then the males
would have a zero probability of being exposed.

## Data Violating This Assumption

      To get a better idea of this assumption–and more specifically what
happens when it is violated–we will generate a simple dataset that
violates this assumption. In this first chunk, we generate our
confounder, *Z*, based off of a random standard normal. Then we begin to
compute the probabilities for our assignment mechanism. The first step
is to generate what we call *x\_lin*, which denotes values that are
linearly related to *Z*, then we exponentiate these values in order to
convert these values to probabilities, which is denoted by *x\_p*.

``` r
set.seed(1000)
weefle <- tibble(
    Z = rnorm(1000, 0, 1), #Confounder
    x_lin = 1 + Z + rnorm(1000, 0, 1), 
    x_p = exp(x_lin)/(1+exp(x_lin)))  #Exposure Probabilities
weefle$x_p[seq(1,1000, by = 4)] = 0 #Assigning some probabilities as 0
```

      It is then this *x*<sub>*p*</sub> that we generate our exposure
assignment from. More specifically, we assign each observation to the
exposure using a random binomial with each observation having its
probability from *x*<sub>*p*</sub>. After we have our exposure
assignments all laid out, we will generate our outcomes. First we
calculate *y*0 and *y*1. In this case we will generate them to be the
same in order for our true causal effect to be zero. We then assign
*y\_obs* as *y*1 when the subject is exposed, or *y*0 if not. Here let
us note that in a true causal setting, we would only know the *y\_obs*
value based on the subjects exposure (i.e. we wouldn’t know what
would’ve happened otherwise). We now will have our data frame! We can
see the first few rows below.

``` r
set.seed(10)
weefle <- weefle %>% mutate(
  x = rbinom(1000, 1, x_p) #Exposure assignments
)

d_ellie <- weefle %>%
  mutate(y0 = Z + x, #Outcomes 
          y1 = Z + x,
          y_obs = ifelse(x == 1, y1, y0))

head(d_ellie)
```

    ## # A tibble: 6 × 7
    ##         Z  x_lin   x_p     x     y0     y1  y_obs
    ##     <dbl>  <dbl> <dbl> <int>  <dbl>  <dbl>  <dbl>
    ## 1 -0.446   2.63  0         0 -0.446 -0.446 -0.446
    ## 2 -1.21   -1.43  0.194     0 -1.21  -1.21  -1.21 
    ## 3  0.0411  1.90  0.870     1  1.04   1.04   1.04 
    ## 4  0.639   2.28  0.907     1  1.64   1.64   1.64 
    ## 5 -0.787  -1.08  0         0 -0.787 -0.787 -0.787
    ## 6 -0.385   0.568 0.638     0 -0.385 -0.385 -0.385

## The Causal Effects

      Since we have *y*0 and *y*1, we are able to calculate the true
causal effect (moreso verify that it is zero since we set *y*0 and *y*1
to be equal). Then we calculate the observed causal effect with our data
that violates the probabilistic assumption. From the calculations below,
we see that our observed causal effect is approximately 1.5. With the
frame of reference that our *y’s* were generated based off of a standard
normal, this observed causal effect is rather far away from our true
effect.

``` r
d_ellie %>%
  summarise(true_causal_effect = mean(y1) - mean(y0))
```

    ## # A tibble: 1 × 1
    ##   true_causal_effect
    ##                <dbl>
    ## 1                  0

``` r
d_ellie %>%
  summarise(observed_causal_effect = 
              sum(y_obs * x) / sum(x) -
              sum(y_obs * (1 - x)) / sum(1 - x))
```

    ## # A tibble: 1 × 1
    ##   observed_causal_effect
    ##                    <dbl>
    ## 1                   1.51

### A Weighted Causal Effect

      So our unweighted causal estimate wasn’t great. We will now create
a propensity score model and use an average treatment effect weight in
order to see if it will yield a causal estimate closer to our true
value. After weighting, we see that our observed causal estimate is
slightly closer to zero, but not too much of a change.

``` r
p <- glm(x ~ Z, data = d_ellie, family = binomial()) %>%
  predict(type = "response")
d_ellie <- d_ellie %>%
  mutate(prob_of_exposure = p)

d_ellie %>%
  mutate(w_ate = x / p + (1 - x) / (1 - p)) %>%
  summarise(observed_causal_effect = 
              sum(y_obs * x * w_ate) / sum(x * w_ate) -
              sum(y_obs * (1 - x) * w_ate) / sum((1 - x) * w_ate))
```

    ## # A tibble: 1 × 1
    ##   observed_causal_effect
    ##                    <dbl>
    ## 1                   1.01

## Discussion

      Overall, the biggest thing to note here is that when we violate
the probabilistic assumption, we are no longer conducting an accurate
causal study. Although the final estimates may not be *obviously* wrong,
we know that it would be invalid to use them to make any sort of causal
conclusion. However, numerically speaking, there may not appear to be
anything wrong with it.
