# Depth-of-Coverage-Per-Interval
GATK DepthOfCoverage, but without merging intervals

## Requirements

* Python. Python version 2.7.11 has been tested.
* 

### sorva.py
'''Options:
  -h, --help            show this help message and exit
  --depthofcoverage DEPTHOFCOVERAGE, -d DEPTHOFCOVERAGE
                        output of GATK -T DepthOfCoverage, containing coverage
                        at each locus
  --intervals INTERVALS, -L INTERVALS
                        file containing list of intervals. Sorted interval
                        file will be created here ending in .sorted
  --no_sort             do not sort interval file (default:False)
  --no_header           depthofcoverage does not contain header line
                        (default:False)
  --output_file OUTPUT_FILE, -o OUTPUT_FILE
                        output file
'''
