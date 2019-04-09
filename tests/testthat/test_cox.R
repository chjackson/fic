
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
X <- newdata_to_X(newdata, wide, intercept=FALSE)

ficall <- fic(wide, inds=inds, inds0=inds0, focus="survival", X=X, t=5)

expect_error(fic(wide, inds=inds, inds0=inds0, focus="survival", X=X, t=-1), "all > 0")
expect_error(fic(wide, inds=inds, inds0=inds0, focus="survival", X=X, t="foo"), "must be numeric")

test_that("Cox FIC results",{
    expect_equal(ficall$FIC[1:8],
                 c(9.56963761576712, 9.62619861138848, 8.60796715809067, 16.3132883942816, 
                   25.4775188260772, 25.2462463739319, 15.136514133193, 9.45357761358747), tol=1e-06)
    expect_equal(ficall$rmse[1:8],
                 c(0.170905480487675, 0.171710775889101, 0.156580994560079, 0.249207826355617, 
0.326815091143021, 0.325084520541464, 0.282320311026908, 0.171966784884666), tol=1e-06)
})

test_that("Built-in focus matches manual focus",{
focus_list <- list(focus = cox_survival,
                   focus_deriv = cox_survival_deriv,
                   focus_dH = cox_survival_dH)
ficall2 <- fic(wide, inds=inds, inds0=inds0, focus=focus_list, X=X, t=5)
expect_equal(ficall, ficall2)
})
