Grad Assignment
================
Ellie Grace Moore
4/18/2022

Hey causal cats! Welcome back to my blog. Today we will choose a
specific causal assumption and then generate data that *violates* this
assumption.

Choose one “Causal Inference Assumption” and write a blog post
describing:

What the assumption is in words What the assumption means mathematically
An example where the assumption is violated Using the simulation
techniques we’ve used in class (for example the “Lucy Land”
simulations") Come up with a data set that violates the assumption
you’ve chosen. Attempt to estimate a causal effect using this dataset
and explain how wrong (or not) the final result is.

``` r
set.seed(1)
weefle <- tibble(
    Z = rnorm(1000, 0, 1), #Confounder
    x_lin = 3 + Z + rnorm(1000, 0, 1), 
    x_p = exp(x_lin)/(1+exp(x_lin)))  #Exposure Probabilities
weefle$x_p[seq(1,1000, by = 4)] = 0 #Assigning some probabilities as 0

weefle <- weefle %>% mutate(
  x = rbinom(1000, 1, x_p) #Exposure assignments
)

d_ellie <- weefle %>%
  mutate(y0 = Z + x,
          y1 = Z + x,
          y_obs = ifelse(x == 1, y1, y0))
```

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
    ## 1                   1.29

``` r
p <- glm(x ~ Z, data = d_ellie, family = binomial()) %>%
  predict(type = "response")
d_ellie <- d_ellie %>%
  mutate(prob_of_exposure = p)
```

``` r
d_ellie %>%
  mutate(w_ate = x / p + (1 - x) / (1 - p)) %>%
  summarise(observed_causal_effect = 
              sum(y_obs * x * w_ate) / sum(x * w_ate) -
              sum(y_obs * (1 - x) * w_ate) / sum((1 - x) * w_ate))
```

    ## # A tibble: 1 × 1
    ##   observed_causal_effect
    ##                    <dbl>
    ## 1                  0.999
