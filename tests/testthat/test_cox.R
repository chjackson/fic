
library(survival)
wide <- coxph(Surv(years, death==1) ~ sex + thick_centred + infilt + epith + ulcer + depth + age, data=melanoma)

inds <- rbind(c(1,1,1,1,1,1,1,1,1,1),
              c(1,1,1,1,1,1,1,1,1,0),
              c(1,1,1,1,1,1,1,0,0,0),
              c(1,1,1,1,1,1,0,0,0,0),
              c(1,1,1,1,1,0,0,0,0,0),
              c(1,1,0,0,0,0,0,0,0,0),
              c(1,0,0,0,0,0,0,0,0,0))
inds0 <- c(1,0,0,0,0,0,0,0,0,0)

newdata <- with(melanoma,
                data.frame(sex = c("female","male"),
                           thick_centred = tapply(thick_centred, sex, mean),
                           infilt=4, epith=1, ulcer=1, depth=2,
                           age = tapply(age, sex, mean)))
X <- newdata_to_X(newdata, wide)

ficall <- fic(wide, inds=inds, inds0=inds0, focus="survival", X=X, t=5)

test_that("Built-in focus matches manual focus",{
focus_list <- list(focus = cox_survival,
                   focus_deriv = cox_survival_deriv,
                   focus_dH = cox_survival_dH)
ficall2 <- fic(wide, inds=inds, inds0=inds0, focus=focus_list, X=X, t=5)
expect_equal(ficall, ficall2)
})
