# DepthOfCoveragePerInterval
Calculate depth of coverage as GATK DepthOfCoverage does, but without merging intervals

## Requirements

* Python. Python version 2.7.11 has been tested.
* Numpy. Version 1.11.0 has been tested.

### Quickstart

python depth_of_coverage_per_interval.py -d HG00096_DepthOfCoverage -o HG00096_DepthOfCoverage.notmerged.sample_interval_summary -ct 1 -ct 4 -ct 10 -L proteindomains.intervals

### Usage

'''Options:
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
'''
