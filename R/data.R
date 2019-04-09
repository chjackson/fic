#' Risk factors associated with low infant birth weight
#' 
#' @format A data frame with 189 rows and 19 variables.  The first 10 columns are included in the dataset of the same name in the \pkg{MASS} package.  The remaining 9 columns are defined in Claeskens and Hjort (2008), and are included in this dataset for convenience.
#' \describe{
#'     \item{low}{indicator of birth weight less than 2.5 kg}
#'     \item{age}{mother's age in years}
#'     \item{lwt}{mother's weight in pounds at last menstrual period}
#'     \item{race}{mother's race (‘1’ = white, ‘2’ = black, ‘3’ = other)}
#'     \item{smoke}{smoking status during pregnancy}
#'     \item{ptl}{number of previous premature labours}
#'     \item{ht}{history of hypertension}
#'     \item{ui}{presence of uterine irritability}
#'     \item{ftv}{number of physician visits during the first trimester}
#'     \item{bwt}{birth weight in grams}
#'   \item{smokeui}{Binary indicator for both smoking and uterine irritation}
#'   \item{smokeage}{Interaction between age and binary (0/1) smoking status, that is, \code{age} for smokers and zero for non-smokers.}
#'   \item{intercpt}{Intercept term (all 1)}
#'   \item{raceother}{Binary indicator for \code{race} "other"}
#'   \item{raceblack}{Binary indicator for \code{race} "black"}
#'   \item{ftv2p}{Binary indicator for \code{ftv}, number of physician visits during the first trimester, 2 or more}
#'   \item{ftv1}{Binary indicator for \code{ftv} 1}
#'   \item{ptd}{Binary indicator for \code{ptl}, number of previous premature labours, 1 or more}
#'   \item{lwtkg}{Weight measured in kg, as used in Claeskens and Hjort.  Note \code{lwt}, as used in \pkg{MASS}, is in pounds.}
#' }
#' @source \pkg{MASS} package (Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer); originally from Hosmer, D.W. and Lemeshow, S. (1989) Applied Logistic Regression. New York: Wiley
#'
#' @references Claeskens, G., & Hjort, N. L. (2008). Model selection and model averaging (Vol. 330). Cambridge: Cambridge University Press.
"birthwt"

#' Malignant melanoma survival data
#'
#' Data originally analysed by Andersen et al (1993) on the survival of 205 patients in Denmark with malignant melanoma, and used by Claeskens and Hjort (2008) to illustrate focused model selection in Cox regression.
#' 
#' @format A data frame with 205 rows and the following columns:
#' \describe{
#'   \item{ptno}{Patient identification number}
#'   \item{death}{Survival status: 1 = dead from illness, 2 = censored, 4 = dead from other causes}
#'   \item{days}{Survival time in days}
#'   \item{depth}{Invasion depth: factor with levels 1, 2, 3}
#'   \item{infilt}{Infection infiltration level, a measure of resistance to the tumour: factor with levels 1 (high resistance), 2, 3, 4 (low resistance)}
#'   \item{epith}{Indicator for epithelioid cells present}
#'   \item{ulcer}{Indicator for ulceration}
#'   \item{thick}{Thickness of the tumour in 1/100 mm}
#'   \item{sex}{Sex. Factor with levels \code{"female","male"} and reference level \code{"female"}}
#'   \item{age}{Age in years}
#'   \item{years}{Survival time in years (instead of days)}
#'   \item{thick_centred}{Version of \code{thick} centred around its mean and rescaled, defined as (\code{thick} - 292)/100. }
#' }
#' @source The supporting material from Claeskens and Hjort (2008), at \url{https://feb.kuleuven.be/public/u0043181/modelselection/datasets/melanoma_data.txt}.    Versions of this dataset are also given in the \pkg{MASS} and \pkg{boot} packages.
#'
#' @references Claeskens, G., & Hjort, N. L. (2008). Model selection and model averaging (Vol. 330). Cambridge: Cambridge University Press.
#'
#' Andersen, P. K., Borgan, O., Gill, R. D., & Keiding, N. (2012). Statistical models based on counting processes. Springer.
#' 
"melanoma"
