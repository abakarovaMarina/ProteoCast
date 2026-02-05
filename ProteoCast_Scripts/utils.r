segmentation <- function(f_proteocast, output, weight = 0.01, alpha = 2) {
  
  data <- read.csv(f_proteocast)
  prot <- sub("_[^_]*_[^_]*$", "", basename(f_proteocast))

  res <- tryCatch({
      # If pLDDT column is not available, we assume no pLDDT values
      if (is.null(data$pLDDT)) {
      # Segment GEMME without pLDDT values
      seg_gemme <- seg(y = data$GEMME_mean, y_plddt = rep(0, length(data$GEMME_mean)), w = weight, a = alpha)
    
      data.frame(
          start = seg_gemme$start,
          end = seg_gemme$end,
          mean = seg_gemme$mu,
          state = seg_gemme$state,
          type = rep("GEMME", length(seg_gemme$mu)),
          protein = prot
      )
      } else {
      # Segment GEMME and pLDDT
      seg_gemme <- seg(y = data$GEMME_mean, y_plddt = data$pLDDT, w = weight, a = alpha)
      seg_plddt <- seg(y = data$pLDDT, y_plddt = data$pLDDT, w = weight, a = alpha)

      data.frame(
          start = c(seg_gemme$start, seg_plddt$start),
          end = c(seg_gemme$end, seg_plddt$end),
          mean = c(seg_gemme$mu, seg_plddt$mu),
          state = c(seg_gemme$state, seg_plddt$state),
          type = rep(c("GEMME", "pLDDT"), c(length(seg_gemme$mu), length(seg_plddt$mu))),
          protein = prot
      )}
    }, error = function(e) {
      # If an error occurs, print the error message and the protein causing it
      message(paste("Error processing protein:", prot))
      message("Error message:", e$message)
      return(NULL)  # Return NULL to skip this protein in the final result
    })

  output <- file.path(output, paste0(prot, "_Segmentation.csv"))
  write.csv(res, output, row.names=FALSE)
}

################################################################################
# SEGMENTATION FUNCTION : WRAPPER AROUND FPOP
################################################################################

seg <- function(y, y_plddt, w, a) {
  
  # Check for NA/NaN/Inf values before proceeding
  if (any(is.na(y)) | any(is.na(y_plddt)) | any(is.infinite(y)) | any(is.infinite(y_plddt))) {
    stop("Input contains NA, NaN, or Inf values.")
  }
  
  # Estimate variance
  if (sum(y_plddt[y!=0] <= 0.70) >= 5) {
    sd_est <- fpopw::sdDiff(y[y_plddt <= 0.70], "HALL")
  } else {
    sd_est <- fpopw::sdDiff(y, "HALL")
  }

  #if (sum(y_plddt <= 0.70) < 5) {
  #  sd_est <- fpopw::sdDiff(y, "HALL")
  #} else {
  #  sd_est <- fpopw::sdDiff(y[y_plddt <= 0.70], "HALL")
  #}
  if (length(y) == 0) stop("Input 'y' has zero length.")
  
  # estimate penalty
  lambda <- a*log(length(y))
  
  # changepoint inference by optimizing likelihood
  cp     <- fpopw::Fpop_w(y/sd_est, ifelse(y_plddt>0.70, w, 1), lambda)$t.est
  
  # estimate the mean for each segment 
  s      <- c(0,cp[-length(cp)])+1
  mu     <- sapply(seq_along(cp), function(x) mean(y[s[[x]]:cp[[x]]]))
  
  # find state of each segment : 0 = valley; 2 = summit; 1 = intermediate steps
  if (length(mu)>1) {
    trends <- rle(diff(mu)>0)
  	states <- c(unlist(sapply(seq_along(trends$lengths), function(i) {
  		rep(
  		  c(2,0,1), 
  		  c(trends$values[[i]]==FALSE, trends$values[[i]]==TRUE, trends$lengths[[i]]-1)
  		)
  	})), ifelse(trends$values[length(trends$values)], 2, 0))
   } else {
    states <- 0
  }
  
  # return the start position, end position, mean and state of segments
  list(start=s, end=cp, mu=mu, state=states)
}
