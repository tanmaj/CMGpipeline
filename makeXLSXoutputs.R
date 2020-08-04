#!/usr/bin/env Rscript
library("optparse")
library("openxlsx")
options(stringsAsFactors = F)

# Parse input options
option_list = list(
  make_option(c("--sample_basename"), type="character", default=NULL, help="Exome id, for example PX5000", metavar="character"),
  make_option(c("--RARE_FUNCTIONAL"), type="character", default=NULL, help="RARE_FUNCTIONAL output", metavar="character"),
  make_option(c("--HET_DOMINANT"), type="character", default=NULL, help="HET_DOMINANT output", metavar="character"),
  make_option(c("--COMPHET_RECESSIVE"), type="character", default=NULL, help="COMPHET_RECESSIVE output", metavar="character"),
  make_option(c("--HOM_RECESSIVE"), type="character", default=NULL, help="HOM_RECESSIVE output", metavar="character"),
  make_option(c("--CLINVAR_PATHOGENIC"), type="character", default=NULL, help="CLINVAR_PATHOGENIC output", metavar="character"),
  make_option(c("--CLINVAR_FILTERED"), type="character", default=NULL, help="CLINVAR_FILTERED output", metavar="character"),
  make_option(c("--CLINVAR_ALL"), type="character", default=NULL, help="CLINVAR_ALL output", metavar="character"),
  make_option(c("--PANEL_FILTERED"), type="character", default=NULL, help="PANEL_FILTERED output", metavar="character"),
  make_option(c("--PANEL_ALL"), type="character", default=NULL, help="PANEL_ALL output", metavar="character"),
  make_option(c("--XLSX_OUTPUT"), type="character", default=NULL, help="Output XLSX", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$RARE_FUNCTIONAL)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Make a list of reports
reportList<-list()
reportList$RARE_FUNCTIONAL <- read.table(opt$RARE_FUNCTIONAL, sep="\t", header=T, quote="", dec = ".", fill=NA)
reportList$HET_DOMINANT <- read.table(opt$HET_DOMINANT, sep="\t", header=T, quote="", dec = ".", fill=NA)
reportList$COMPHET_RECESSIVE <- read.table(opt$COMPHET_RECESSIVE, sep="\t", header=T, quote="", dec = ".", fill=NA)
reportList$HOM_RECESSIVE <- read.table(opt$HOM_RECESSIVE, sep="\t", header=T, quote="", dec = ".", fill=NA)
reportList$CLINVAR_PATHOGENIC <- read.table(opt$CLINVAR_PATHOGENIC, sep="\t", header=T, quote="", dec = ".", fill=NA)
reportList$CLINVAR_FILTERED <- read.table(opt$CLINVAR_FILTERED, sep="\t", header=T, quote="", dec = ".", fill=NA)
reportList$CLINVAR_ALL <- read.table(opt$CLINVAR_ALL, sep="\t", header=T, quote="", dec = ".", fill=NA)
if(!is.null(opt$PANEL_FILTERED)) reportList$PANEL_FILTERED <- read.table(opt$PANEL_FILTERED, sep="\t", header=T, quote="", dec = ".", fill=NA)
if(!is.null(opt$PANEL_ALL)) reportList$PANEL_ALL <- read.table(opt$PANEL_ALL, sep="\t", header=T, quote="", dec = ".", fill=NA)

# Prepare outputs
hs <- createStyle(fontColour = "#ffffff", fgFill = "#4F80BD",
                  halign = "center", valign = "center", textDecoration = "Bold",
                  border = "TopBottomLeftRight", textRotation = 45)
bodyStyle <- createStyle(border="TopBottom", borderColour = "#4F81BD")

openxlsx::write.xlsx(reportList, file = opt$XLSX_OUTPUT, borders = "rows", headerStyle = hs, addFilter = T)


