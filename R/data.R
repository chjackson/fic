#' Risk Factors Associated with Low Infant Birth Weight 
#' 
#' @format A data frame with 189 rows and 19 variables.  First 10 columns documented in \code{\link[MASS]{birthwt}}.  Additional 9 columns defined in Claeskens and Hjort (2008).
#' \describe{
#'   \item{}{}
#'   \item{}{}
#'   ...
#' }
#' @source \pkg{MASS} package (Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer); originally from Hosmer, D.W. and Lemeshow, S. (1989) Applied Logistic Regression. New York: Wiley
#'
#' @references Claeskens, G., & Hjort, N. L. (2008). Model selection and model averaging (Vol. 330). Cambridge: Cambridge University Press.
"birthwt"

#' Melanoma
#' 
#' @format A data frame with 205 rows and the following columns:
#' \describe{
#'   \item{ptno}{Patient identification number}
#'   \item{death}{Survival status: 1 2 3 4 }
#'   \item{days}{Survival or censoring time}
#'   \item{depth}{}
#'   \item{infilt}{}
#'   \item{mucous}{}
#'   \item{epith}{}
#'   \item{ulcer}{}
#'   \item{thick}{}
#'   \item{sex}{Sex. Factor with levels \code{"female","male"} and reference level \code{"female"}}
#'   \item{age}{Age in years}
#'   \item{qq}{}
#'   \item{years}{Survival time in years (instead of days)}
#'   \item{thick_centred}{Centred version of \code{thick}, defined as }
#'   ...
#' }
#' @source The supporting material from Claeskens and Hjort (2008), at \url{https://feb.kuleuven.be/public/u0043181/modelselection/datasets/melanoma_data.txt}.
#'
#' @references Claeskens, G., & Hjort, N. L. (2008). Model selection and model averaging (Vol. 330). Cambridge: Cambridge University Press.
"melanoma"
