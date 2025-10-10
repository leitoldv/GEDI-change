# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# List of CRAN packages to be installed
cran_packages <- c(
  "s3", "foreach", "aws.s3", "stringr", "optmatch", "doParallel", "arrow"
)

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

# Load the packages
#library("s3")
