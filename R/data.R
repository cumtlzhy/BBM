#' ip
#'
#' ip and input are simulated data, simulating reads count data from a 100-site MeRIP-Seq experiment with 15 sample conditions..
#' Three biclusters are implanted in the data, the sizes of which are 15×5,10×5 and 15×7, respectively.
#' The number of reads for each site under each condition in the IP matrix is generated following (2) in the paper BBM.
#' The parameters p of the binomial distribution, which the three biclusters and background follow, are set to 0.3, 0.98, 0.66, and 0.5 respectively,
#' to make the whole distri-bution statistic of the simulated data similar to the real data.The three biclusters did not overlap in sites.
#' In con-ditions, there are 2 overlaps between the first two biclus-ters, and 1 overlaps between the last two biclusters
#'
#'
#' @examples
#'   head(ip)
"ip"

#' input
#'
#' ip and input are simulated data, simulating reads count data from a 100-site MeRIP-Seq experiment with 15 sample conditions..
#' Three biclusters are implanted in the data, the sizes of which are 15×5,10×5 and 15×7, respectively.
#' The number of reads for each site under each condition in the IP matrix is generated following (2) in the paper BBM.
#' The parameters p of the binomial distribution, which the three biclusters and background follow, are set to 0.3, 0.98, 0.66, and 0.5 respectively,
#' to make the whole distri-bution statistic of the simulated data similar to the real data.The three biclusters did not overlap in sites.
#' In con-ditions, there are 2 overlaps between the first two biclus-ters, and 1 overlaps between the last two biclusters
#'
#'
#' @examples
#'   head(input)
"input"
