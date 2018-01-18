if (requireNamespace("MASS", quietly = TRUE)){

birthwt <- within(MASS::birthwt, { 
    lwtkg = lwt*0.45359237 # weight measured in lb in MASS data, kg in FIC book
    race = factor(race, labels = c("white", "black", "other"))
    ptd = factor(ptl > 0)
    ftv = factor(ftv) # keeping all levels, not transforming to 2+
    levels(ftv)[-(1:2)] <- "2+" # doing the transformation
    ftv1=as.numeric(ftv=="1")
    ftv2p=as.numeric(ftv=="2+")
    low = factor(low)
    smoke = (smoke > 0)
    ht = (ht > 0)
    ui = (ui > 0)
    raceblack <- 1*(race=="black")
    raceother <- 1*(race=="other")
    intercpt = 1 + 0*seq(1:length(low)) # intercept
    smokeage=age*smoke
    smokeui=smoke*ui
}
)

use_data(birthwt, pkg="..")
    
}

