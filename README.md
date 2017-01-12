# DepthOfCoveragePerInterval

If you would like to calculate mean depth of coverage for custom intervals, GATK DepthOfCoverage -L has a serious limitation: it automatically merges overlapping intervals, making it impossible to calculate coverage for intervals such as known protein domains, which frequently overlap with each other.

This python tool calculates depth of coverage as GATK DepthOfCoverage -L does, but without merging intervals. You will need to run this tool AFTER running GATK DepthOfCoverage initially.

Note: The software does not do sample-by-sample comparisons; it only takes into account the total and mean coverage columns from the GATK DepthOfCoverage output file.

## Requirements

* Python. Python version 2.7.11 has been tested.
* Numpy. Version 1.11.0 has been tested.

### Generate initial DepthOfCoverage file
You will need .bam and .bai alignment files. As an example, you can download the mapped exome alignment file for a 1000 Genomes Project sample here: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/ 

Next, use GATK to generate a DepthOfCoverage file that contains the depth of coverage for each locus. E.g.:

    java -jar /path/to/gatk-current/target/GenomeAnalysisTK.jar \
      -T DepthOfCoverage \
      -R /path/to/gatk_resource_bundle/current/b37/human_g1k_v37_decoy.fasta \
      -o HG00096_DepthOfCoverage \
      -I /path/to/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam \
      -ct 1 -ct 4 -ct 10 \
      -L /path/to/interpro_annots_all_genetrack.interval_list
      
You will use the .sample_interval_summary output from the above command as input when running DepthOfCoveragePerInterval.

## Quickstart

    python depth_of_coverage_per_interval.py --help
    python depth_of_coverage_per_interval.py \
    -d HG00096_DepthOfCoverage \
    -o HG00096_DepthOfCoverage.notmerged.sample_interval_summary \
    -ct 1 -ct 4 -ct 10 \
    -L proteindomains.intervals

## Usage

```Options:
  -h, --help            show this help message and exit
  --depthofcoverage DEPTHOFCOVERAGE, -d DEPTHOFCOVERAGE
                        output of GATK -T DepthOfCoverage, containing coverage
                        at each locus
  --intervals INTERVALS, -L INTERVALS
                        file containing list of intervals. Sorted interval
                        file will be created here ending in .sorted
  --summaryCoverageThreshold CT, -ct CT
                        report the percentage of bases covered to an amount equal 
                        to or greater than this number (e.g. % bases >= CT)
  --no_sort             do not sort interval file (default:False)
  --no_header           depthofcoverage does not contain header line
                        (default:False)
  --output_file OUTPUT_FILE, -o OUTPUT_FILE
                        output file
```
## Example usage
    python depth_of_coverage_per_interval.py -d /path/to/dir/HG00096_DepthOfCoverage -L /path/to/interpro_annots_all_genetrack_nopatches.intervals -o /path/to/HG00096_DepthOfCoverage.notmerged.sample_interval_summary
