library(argparse)
library(data.table)

#,Define command-line arguments
parser <- ArgumentParser(description='Computes Dxy from window-sfs file.')
parser$add_argument('--wsfs', required=TRUE, help='Multidimensional SFS file (MSFS / DSFS) with 2 populations and windows stored as replicate SFS (e.g. multiple lines) [required]')
parser$add_argument('--npop1', type="integer", required=TRUE, help='Number of diploid individuals in population 1')
parser$add_argument('--npop2', type="integer", required=TRUE, help='Number of diploid individuals in population 2')
parser$add_argument('--out', required=TRUE, help='name of outfile [required]')
args <- parser$parse_args()

# Read SFS file
wsfs <- fread(args$wsfs, fill = TRUE)
nsfs <- dim(wsfs)[1]
pop_sizes <- c(2*as.numeric(args$npop1), 2*as.numeric(args$npop2))

# Calculate weighting vector for SFS entries
idp <- expand.grid(0:pop_sizes[1], 0:pop_sizes[2])
wsfs_width <- dim(idp)[1]
wei <- apply(idp, 1, 
             function(x) 
               mean(abs(diff(
                 rbind(c(rep(0, x[1]), 
                         rep(1, pop_sizes[1] - x[1])), 
                       c(rep(0, x[2]), 
                         rep(1, pop_sizes[2] - x[2])))))))

head(wei)
# Loop through SFS file and calculate, dxy
out <- wsfs[,1:3]
for (i in 1:dim(wsfs)[1]) {
  sfs <- wsfs[i, 4:(wsfs_width +3)]
  out$nsites[i] <- sum(sfs, na.rm = T)
  out$nsnps[i] <- sum(sfs, na.rm = T) - sfs[[1]]
  
  if (sum(sfs, na.rm = T) > 0) {
    out$dxy[i] <- sum(wei * sfs, na.rm = T) / sum(sfs, na.rm = T)
  } else {
    out$dxy[i] <- NA
  }
}
fwrite(out, file = args$out, sep = "\t")
