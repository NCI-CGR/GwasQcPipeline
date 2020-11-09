# Parse IDAT files and calculate median intensity.
require(illuminaio)

calc_median_intensity <- function(redIdat, greenIdat, output){

    red <- readIDAT(redIdat)
    green <- readIDAT(greenIdat)

    intensity <- red$Quants[,1] + green$Quants[,1]
    med_intensity <- median(intensity)

    write.table(med_intensity, file = output, quote = F, col.names = F, row.names = F)

}

if(exists("snakemake")){
    # The script was run using snakemake

    calc_median_intensity(snakemake@input[["red"]], snakemake@input[["green"]], snakemake@output[[1]])

} else {
    # The script was run from command line

    library("optparse")
    option_list = list(
        make_option("--red", type="character", default=NULL, help="Red intensities.", metavar="character"),
        make_option("--green", type="character", default=NULL, help="Green intensities.", metavar="character"),
        make_option("--out", type="character", default=NULL, help="Output file", metavar="character")
    )
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)

    calc_median_intensity(opt$red, opt$green, opt$out)

}
