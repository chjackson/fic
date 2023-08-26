test_that("example code runs",{
  expect_error({
    wide.glm <- glm(low ~ lwtkg + age + smoke + ht + ui + smokeage + smokeui,
                    data=birthwt, family=binomial)
    vals.smoke <-    c(1, 58.24, 22.95, 1, 0, 0, 22.95, 0)
    vals.nonsmoke <- c(1, 59.50, 23.43, 0, 0, 0, 0, 0)
    X <- rbind("Smokers" = vals.smoke, "Non-smokers" = vals.nonsmoke)
    inds0 <- c(1,1,0,0,0,0,0,0)
    combs <- all_inds(wide.glm, inds0)
    ficres <- fic(wide = wide.glm, inds = combs, inds0 = inds0,
                  focus = prob_logistic, X = X)
    plot(ficres)
    ggplot_fic(ficres)
  },NA)
})
