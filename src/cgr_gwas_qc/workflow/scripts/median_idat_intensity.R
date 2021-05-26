# Parse IDAT files and calculate median intensity.
require(illuminaio)

calc_median_intensity <- function(
                                  sample_id,
                                  red_idat,
                                  green_idat,
                                  output) {

  # Calculate median intenisty of red and green channels
  red <- illuminaio::readIDAT(red_idat)
  green <- illuminaio::readIDAT(green_idat)

  intensity <- red$Quants[, 1] + green$Quants[, 1]
  med_intensity <- median(intensity)

  # Save a csv "Sample_ID,Chip_ID,median_intensity"
  df <- data.frame(
    Sample_ID = sample_id,
    median_intensity = med_intensity
  )

  write.csv(df, file = output, quote = F, row.names = F)
}

if (exists("snakemake")) {
  sample_id <- snakemake@wildcards[["Sample_ID"]]

  calc_median_intensity(
    sample_id,
    snakemake@input[["red"]],
    snakemake@input[["green"]],
    snakemake@output[[1]]
  )

} else {

  args = commandArgs(trailingOnly=TRUE)
  calc_median_intensity(
    args[1],  # sample_id
    args[2],  # red pattern
    args[3],  # green pattern
    args[4]  # outfile
  )

}
