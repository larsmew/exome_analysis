args = commandArgs(trailingOnly=TRUE)

backupPackages <- function(filename = "R_packages"){
    tmp <- installed.packages()
    
    # Save list as R object for easy reinstall
    installedpkgs <- as.vector(tmp[,"Package"])
    save(installedpkgs, file=paste0(filename, ".rda"))
    
    # Save list of packages and version number as text file
    installedpkgs <- tmp[,c("Package", "Version")]
    write.table(installedpkgs, file=paste0(filename, ".txt"), row.names = F, col.names = F, sep = "=", quote = F)
}

# Update all packages from all sources
updateAllPackages <- function(){
    # Update from bioconductor    
    source("https://bioconductor.org/biocLite.R")
    biocLite(ask = FALSE)
    
    # Update from CRAN
    update.packages(ask = FALSE, repos = "https://cran.rstudio.com")
}

if (length(args) == 1){
    filename = args[1]
} else {
    filename = "R_packages"
}
backupPackages(filename)

updateAllPackages()
