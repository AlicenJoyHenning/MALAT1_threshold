#' Define MALAT1 Expression Threshold
#'
#' This function calculates a minimum threshold value for MALAT1 expression. 
#' Given a vector of normalized MALAT1 expression values, this function identifies 
#' an optimal threshold to aid in filtering empty droplets in single-cell RNA sequencing (scRNA-seq) data.
#'
#' @name define_malat1_threshold
#' @param counts A numeric vector of normalized MALAT1 counts.
#' @param bw A numeric value specifying the bandwidth for density estimation. Default is 0.1.
#' @param lwd A numeric value specifying the line width for the red line in the histogram. Default is 2.
#' @param breaks An integer specifying the number of bins for the histogram. Default is 100.
#' @param chosen_min A numeric value specifying the minimum threshold for the MALAT1 peak. Default is 1.
#' @param smooth A numeric value specifying the smoothing parameter for the smooth spline. Default is 1.
#' @param abs_min A numeric value specifying the absolute minimum for the local minima. Default is 0.3.
#' @param rough_max A numeric value specifying the rough maximum for the local maxima. Default is 2.
#'
#' @import graphics 
#' @import stats
#' @importFrom utils globalVariables
#' 
#' @return A numeric value representing the calculated MALAT1 threshold.
#' @export
#'
#' @examples
#' counts <- rnorm(1000, mean = 1, sd = 0.5)
#' define_malat1_threshold(counts)

utils::globalVariables(c(
  "density_data", "fit", "first_derivative", "local_maxima", "local_minima",
  "biggest_y", "maxi", "mini", "delta", "df", "subset_df", "quad_model",
  "coefficients", "c", "b", "a", "discriminant", "x_intercept1", "x"
))

define_malat1_threshold <- function(counts, bw = 0.1, lwd = 2, breaks = 100,
                                    chosen_min = 1, smooth = 1, abs_min = 0.3,
                                    rough_max = 2) {
  tryCatch({
    # Calculate the density values
    density_data <- density(counts, bw = bw,
                            from = min(counts),
                            to = max(counts))
    # Fit a smooth line to the data
    fit <- smooth.spline(density_data$x, density_data$y, spar = smooth)
    # Predict the first derivative
    first_derivative <- predict(fit, density_data$x, deriv = 1)

    # Find the local maxima
    local_maxima <- density_data$x[which(diff(sign(first_derivative$y)) == -2)]
    if(length(local_maxima) == 0) {
      # Find x-val closest to rough_max
      local_maxima <- density_data$x[which(abs(density_data$x - rough_max) == min(abs(density_data$x - rough_max)))]
    }
    # Plot the density with the local maxima
    plot(density_data, main = "Density Plot with Local Maxima")
    lines(fit, col = "blue")
    points(local_maxima, predict(fit, local_maxima)$y, col = "red", pch = 19)

    # Find the local minima
    local_minima <- density_data$x[which(diff(sign(first_derivative$y)) == 2)]
    if(length(local_minima) == 0) {
      local_minima <- abs_min
    }
    # Plot the density with the local minima
    plot(density_data, main = "Density Plot with Local Minima")
    lines(fit, col = "blue")
    points(local_minima, predict(fit, local_minima)$y, col = "red", pch = 19)

    # Find biggest local maximum greater than desired value (default 1)
    biggest_y <- max(density_data$y[density_data$x %in% local_maxima[local_maxima > chosen_min]])
    maxi <- density_data$x[density_data$y == biggest_y]

    # Find local minimum closest to the left of that
    local_minima <- local_minima[local_minima < maxi]
    if(length(local_minima) < abs_min) {
      local_minima <- abs_min
    }
    mini <- local_minima[(maxi - local_minima) == min(maxi - local_minima)]

    # Calculate delta to get range of values to isolate for quadratic calculation
    delta <- maxi - mini

    # Subset dataframe
    df <- as.data.frame(cbind(density_data$x, density_data$y))
    colnames(df) <- c("x","y")
    subset_df <- df[df[,1] >= (maxi - delta) & df[,1] <= (maxi + delta), ]

    # Fit a quadratic model (y ~ x + I(x^2))
    quad_model <- lm(y ~ poly(x, 2, raw = TRUE), data = subset_df)

    # Plot quadratic
    plot(df$x, df$y,
         xlab = "Normalized MALAT1 expression",
         ylab = "Density value")
    # Add the extracted subset data
    points(subset_df$x, subset_df$y, pch = 16, col = "blue")
    # Add the fitted quadratic curve
    curve(predict(quad_model, newdata = data.frame(x = x)), add = TRUE, col = "red", lwd = 2)

    # Grab intercept
    # Extract coefficients from the summary
    coefficients <- summary(quad_model)$coefficients
    # Extract coefficients
    c <- coefficients[1, 1]
    b <- coefficients[2, 1]
    a <- coefficients[3, 1]
    # Calculate discriminant
    discriminant <- b^2 - 4 * a * c
    # Grab first intercept
    x_intercept1 <- (-b + sqrt(discriminant)) / (2 * a)
    if(x_intercept1 < 0) {
      x_intercept1 <- abs_min
    }

    # Plot histogram with threshold
    hist(counts, breaks = breaks,
         xlab = "Normalized MALAT1 expression",
         ylab = "Number of cells")
    abline(v = x_intercept1, col = "red", lwd = 2)

    return(x_intercept1)
  }, error = function(e) {
    # Code to execute if an error occurs
    message(" An error occurred: Please make sure you have use a vector of normalized counts as input. This may also indicate that you have no high MALAT1 peaks, meaning this particular sample may be poor quality. Please check your histogram of normalized MALAT1 counts to investigate: ", e$message)
    hist(counts, breaks = breaks,
         xlab = "Normalized MALAT1 expression",
         ylab = "Number of cells")
    return(2)
  })
}
