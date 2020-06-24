#!/usr/bin/env Rscript
library("optparse")
library("openxlsx")

options(stringsAsFactors = F)

option_list = list(
  make_option(c("--sample_basename"), type="character", default=NULL, help="Exome id, for example PX5000", metavar="character"),
  make_option(c("--RARE_FUNCTIONAL"), type="character", default=NULL, help="RARE_FUNCTIONAL output", metavar="character"),
  make_option(c("--XLSX_OUTPUT"), type="character", default=NULL, help="Output XLSX", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$RARE_FUNCTIONAL)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

reportList<-list()
reportList$RARE_FUNCTIONAL <- read.table(opt$RARE_FUNCTIONAL, sep="\t", header=T, quote="", dec = ".", fill=NA)

hs <- createStyle(fontColour = "#ffffff", fgFill = "#4F80BD",
                  halign = "center", valign = "center", textDecoration = "Bold",
                  border = "TopBottomLeftRight", textRotation = 45)
bodyStyle <- createStyle(border="TopBottom", borderColour = "#4F81BD")

openxlsx::write.xlsx(reportList, file = opt$XLSX_OUTPUT, borders = "rows", headerStyle = hs, addFilter = T)


