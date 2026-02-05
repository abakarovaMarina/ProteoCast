source("/opt/ServersScripts/utils.r")

################################################################################
# EXPLORING MORE SEGMENTATIONS
################################################################################

# 1ST VERSION


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 1) {
  stop("Please provide ProteoCast.csv file")
}

# Assign arguments to variables
proteocast <- args[1]
proteocast_out <- paste(head(strsplit(args[1], "/")[[1]], -1), collapse = "/")



# contribution weight of large plddt values in loss = 0.01 and penalty hyperparameter alpha = 2
segmentation(
  f_proteocast   = proteocast,
  output    = proteocast_out,
  weight    = 0.1,
  alpha     = 1.4
)


