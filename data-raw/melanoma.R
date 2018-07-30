
melanoma <- read.table("https://feb.kuleuven.be/public/u0043181/modelselection/datasets/melanoma_data.txt",
                       col.names=c("ptno","death","days","depth","infilt","mucous","epith","ulcer","thick","sex","age","qq"),
                       colClasses=c("numeric","factor","numeric","factor","factor","factor","factor","factor","numeric","factor","numeric","factor"))

melanoma$years <- melanoma$days/365.25
melanoma$thick_centred <- (melanoma$thick - 292)/100
melanoma$infilt <- ordered(melanoma$infilt)
melanoma$depth <- ordered(melanoma$depth)
melanoma$epith <- as.numeric(melanoma$epith==1)
melanoma$ulcer <- as.numeric(melanoma$ulcer==1)

levels(melanoma$sex) <- c("female","male")

use_data(melanoma, pkg="..")
