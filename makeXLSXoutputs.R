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
  make_option(c("--MITOMAP"), type="character", default=NULL, help="MITOMAP output", metavar="character"),
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
if(!is.null(opt$MITOMAP)) reportList$MITOMAP <- read.table(opt$MITOMAP, sep="\t", header=T, quote="", dec = ".", fill=NA)

wb <- openxlsx::createWorkbook()
for(sheetName in names(reportList)){
  openxlsx::addWorksheet(wb, sheetName, zoom=90)
  openxlsx::writeData(wb, sheetName, reportList[[sheetName]])
  openxlsx::addFilter(wb, sheetName, row = 1, cols = 1:ncol(reportList[[sheetName]]))
  
  # Set the widths of columns and hide selected columns
  ColNames<-colnames(reportList[[sheetName]])
  
  ColWidths<-rep(8, ncol(reportList[[sheetName]]))
  ColWidths[grepl("CHR",ColNames)]<-6
  ColWidths[grepl("POS",ColNames)]<-12
  ColWidths[grepl("REF",ColNames)]<-5
  ColWidths[grepl("ALT",ColNames)]<-5
  ColWidths[grepl("QUAL",ColNames)]<-5
  ColWidths[grepl("..GT",ColNames)]<-3
  ColWidths[grepl("..AD",ColNames)]<-5
  ColWidths[grepl("..DP",ColNames)]<-4
  ColWidths[grepl("..GQ",ColNames)]<-4
  ColWidths[grepl("ANN....GENE",ColNames)]<-6
  ColWidths[grepl("Disease_name",ColNames)]<-24
  ColWidths[grepl("Categorization",ColNames)]<-14
  ColWidths[grepl("Inheritance",ColNames)]<-6
  ColWidths[grepl("HPO",ColNames)]<-20
  ColWidths[grepl("OMIM",ColNames)]<-20
  ColWidths[grepl("FEATUREID",ColNames)]<-12
  ColWidths[grepl("HGVS_C",ColNames)]<-12
  ColWidths[grepl("HGVS_P",ColNames)]<-12
  ColWidths[grepl("EFFECT",ColNames)]<-12
  ColWidths[grepl("SLOpopulation",ColNames)]<-4
  ColWidths[grepl("gnomAD",ColNames)]<-4
  ColWidths[grepl("CLNSIG",ColNames)]<-4
  ColWidths[grepl("CLNDN",ColNames)]<-18
  ColWidths[grepl("SpliceAI",ColNames)]<-18
  
  ColHide<-rep(FALSE, ncol(reportList[[sheetName]]))
  ColHide[grepl("RANK",ColNames)]<-TRUE
  
  openxlsx::setColWidths(wb, sheetName, cols=1:length(ColWidths), widths = ColWidths, hidden = ColHide, ignoreMergedCells = FALSE)
  
  ## Styling the sheet
  headerStyle <- createStyle(fontColour = "#636363", fgFill = "#9ecae1",
                             halign = "center", valign = "center", textDecoration = "Bold",
                             border = "TopBottom", textRotation = 25)
  
  addStyle(wb, sheet = sheetName, headerStyle, rows = 1, cols = 1:length(ColWidths), gridExpand = TRUE)
  bodyStyle <- createStyle(border = "TopBottom", borderColour = "#4F81BD")
  addStyle(wb, sheet = sheetName, bodyStyle, rows = 2:nrow(reportList[[sheetName]]), cols = 1:ncol(reportList[[sheetName]]), gridExpand = TRUE)
}

openxlsx::saveWorkbook(wb, file = opt$XLSX_OUTPUT, overwrite = TRUE)



