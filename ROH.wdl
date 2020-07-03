version 1.0
## Copyright CMG@KIGM, Ales Maver

task calculateBAF {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    String sample_basename

    File reference_fa

    File dbSNPcommon_bed
    File dbSNPcommon_bed_index

    # Runtime parameters
    String docker
  }
  
  command <<<
  set -e
  
  # The perl script converts mpileup to BAF wiggle file for IGV display
  echo '''#!/usr/bin/perl
  use strict;
  use warnings;
  use Getopt::Long qw(:config pass_through no_ignore_case);

  my ($min_reads) = (10);
  GetOptions (
          "min-reads:s" => \$min_reads,
  );

  while (<>) {
          my $line = $_;
          my @columns = split("\t",$line);
          my $chr = $columns[0];
          my $start = $columns[1];
          my $end = $start + 1;
          my $num_reads = $columns[3];
          my $calls = $columns[4];
          my $id = "mpileup_number_" . $.;
          if($num_reads < $min_reads){ # not enough coverage to have good confidence in the call
                  next;
          }
          my $num_ref = 0;
          while ($calls =~ /[,.]/g) { $num_ref++ }
          my $num_var = $num_reads - $num_ref;
          my $varAlleleFreq = ($num_var/$num_reads)*100;
          print("$chr\t$start\t$end\t$varAlleleFreq\n");
          }''' > mpileupToWig.pl

  samtools mpileup -q 15 -Q20 -f ~{reference_fa} -l ~{dbSNPcommon_bed} ~{input_bam} | perl mpileupToWig.pl --min-reads=20 > ~{sample_basename}.BAF.wig  
  >>>

  runtime {
    docker: docker
  }
  output {
    File output_BAF = "~{sample_basename}.BAF.wig"
  }
}

task CallROH {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    String sample_basename

    File reference_fa

    File dbSNPcommon_bed
    File dbSNPcommon_bed_index

    File gnomAD_vcf
    File gnomAD_vcf_index

    # Runtime parameters
    String docker
  }
  
  command <<<
  set -e
  bcftools mpileup -q 15 -Q20 -f ~{reference_fa} -T ~{dbSNPcommon_bed} ~{input_bam} | bcftools call -m | bcftools view -i 'DP>10 && QUAL>100' -V indels -Ob -o ~{sample_basename}.dbSNP.vcf.bgz
  tabix -p vcf ~{sample_basename}.dbSNP.vcf.bgz
  bcftools annotate -a ~{gnomAD_vcf} -c AF ~{sample_basename}.dbSNP.vcf.bgz -Ob -o ~{sample_basename}.dbSNP.AF.vcf.bgz
  bcftools --AF-tag AF -G30 -I ~{sample_basename}.dbSNP.AF.vcf.bgz | grep "^[^#]" | grep "^RG" | awk -F'\t' '{if($7>20 && $8>30 && $6>1000000)print $3,$4,$5,$8}' OFS='\t' > ~{sample_basename}.ROHcalls.wig
  >>>

  runtime {
    docker: docker
  }
  output {
    File ROH_calls = "~{sample_basename}.ROHcalls.wig"
    File BAF_vcf = "~{sample_basename}.dbSNP.AF.vcf.bgz"
  }
}

task CallPlink {
  input {
    # Command parameters
    File input_vcf
    String sample_basename

    # Runtime parameters
    String docker
  }
  
  command <<<
  set -e
  /usr/local/bin/plink1.9 --vcf ~{input_vcf} --make-bed --out ~{sample_basename} --no-sex --no-parents --no-fid --no-pheno --allow-extra-chr
  /usr/local/bin/plink1.9 -bfile ~{sample_basename} --homozyg --allow-extra-chr --out ~{sample_basename}
  >>>

  runtime {
    docker: docker
  }
  output {
    File ROHplink_calls = "~{sample_basename}.hom"
    File ROHplink_calls_windows = "~{sample_basename}.hom.summary"
    File ROHplink_calls_summary = "~{sample_basename}.hom.indiv"
  }
}