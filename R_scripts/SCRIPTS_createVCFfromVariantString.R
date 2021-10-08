# Load required libraries
library("optparse")

# Parse command line options
option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Sample identifier", metavar="character"),
  make_option(c("-v", "--variant"), type="character", default=NULL,
              help="External variant string in format chr7 117171161 A G HET", metavar="character")
)
opt_parser = optparse::OptionParser(option_list=option_list); opt = optparse::parse_args(opt_parser)

createVCFfromVariants <- function(VARIANT, EXOME, VCF_OUTDIR){
  # Test parameters
  # EXOME="PX0000"
  # VARIANT = "chr2 167055804 AACA TCC HET"
  # VARIANT = "chr2 167055804 AACA TCC HET"
  # VARIANT = "chr2 167055804 AACA TCC HET;chr2 167055804 AACA TCC HET"
  
  if(missing(VCF_OUTDIR)) {VCF_OUTDIR = ""}
  VCF<-paste0(VCF_OUTDIR, EXOME, ".IMPORTVARIANT.vcf")
  
  # Support multiple variant imports added
  VARIANT<-strsplit(VARIANT, split=";")[[1]]
  VARIANT<-as.list(VARIANT)
  VARIANT<-lapply(VARIANT, FUN = trimws )
  VARIANT<-lapply(VARIANT, FUN = function(x){ strsplit(x, split=" ")[[1]] } )
  VARIANT<-do.call(rbind, VARIANT)
  VARIANT<-data.frame(VARIANT)
  names(VARIANT)<-c("CHR", "POS", "REF", "ALT", "GT")
  
  for(i in 1:nrow(VARIANT)){
    if(VARIANT$GT[i]=="HET") VARIANT$GT[i] <- paste0(gsub("HET", "0/1", VARIANT$GT[i]), ":50,50:100:99:407,36,0")
    if(VARIANT$GT[i]=="HOM") VARIANT$GT[i] <- paste0(gsub("HOM", "1/1", VARIANT$GT[i]), ":0,100:100:99:407,36,0")
    
    if(VARIANT$REF[i]=="-") VARIANT$REF[i]<-"."
    if(VARIANT$ALT[i]=="-") VARIANT$ALT[i]<-"."
  }
  
  QUAL = "1000"
  INFO = "AC=1"
  GTFIELDS = "GT:AD:DP:GQ:PL"
  FILTER = "PASS"
  
  write("##fileformat=VCFv4.2", file=VCF, append = T)
  write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">', file=VCF, append = T)
  write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">', file=VCF, append = T)
  write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">', file=VCF, append = T)
  write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=VCF, append = T)
  write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">', file=VCF, append = T)
  write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">', file=VCF, append = T)
  write(paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", EXOME, ".IMPORTVARIANT"), file=VCF, append = T)
  write(paste(VARIANT$CHR, VARIANT$POS, ".", VARIANT$REF, VARIANT$ALT, QUAL, FILTER, INFO, GTFIELDS, VARIANT$GT, sep="\t"), file=VCF, append = T)
  
  return(VCF)
}

createVCFfromVariants(VARIANT = opt$variant, EXOME=opt$sample)