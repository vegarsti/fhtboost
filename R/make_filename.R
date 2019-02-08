make_filename <- function(directory, descriptor, seed) {
  seed_string <- formatC(seed, width=4, flag="0")
  filename <- paste(seed_string, descriptor, sep="_")
  full_filename <- paste(directory, filename,  '.csv', sep='')
  return(full_filename)
}
