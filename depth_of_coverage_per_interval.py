#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
from collections import deque

class DepthOfCoveragePerIntervalError(Exception): pass

#-d /home/aliz/resources/HG00096_DepthOfCoverage -L /scratch1/tmp/aliz/resources/1kg/interpro_annots_all_genetrack_nopatches.intervals -o /home/aliz/resources/HG00096_DepthOfCoverage.notmerged.sample_interval_summary
#-d /home/aliz/resources/HG00096_DepthOfCoverage_1Mrows -L /home/aliz/resources/test.intervals -o /home/aliz/resources/HG00096_DepthOfCoverage.notmerged.test.sample_interval_summary

def file_exists(file_path):
    if isinstance(file_path, str) and os.path.isfile(file_path):
        try:
            open(file_path)
            return True
        except IOError as e:
            return False
    else:
        return False

def check_file_exists(file_path):
    if not file_exists(file_path):
        raise DepthOfCoveragePerIntervalError('file not found: ' + file_path)

chroms = {'': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, '10': 10, '11': 11, '12': 12,
          '13': 13, '14': 14, '15': 15, '16': 16,
          '17': 17, '18': 18, '19': 19, '20': 20, '21': 21, '22': 22, 'X': 23, 'Y': 24, 'MT': 25, 'GL000201.1': 26,
          'GL000242.1': 27, 'GL000237.1': 28,
          'GL000204.1': 29, 'GL000191.1': 30, 'GL000228.1': 31, 'GL000209.1': 32, 'GL000218.1': 33, 'GL000213.1': 34,
          'GL000205.1': 35, 'GL000219.1': 36,
          'GL000223.1': 37, 'GL000212.1': 38, 'GL000222.1': 39, 'GL000194.1': 40, 'GL000192.1': 41}

def compare_intervals(i1, i2):
    c1, pos1 = i1.split(':')
    c2, pos2 = i2.split(':')
    spos1, epos1 = pos1.split('-')
    spos2, epos2 = pos2.split('-')
    return cmp(chroms[c1], chroms[c2]) or cmp(int(spos1), int(spos2))

def output_final_values(i, l, ct_list, out_file):
    cvg_array = np.array([cvg for c, p, total, cvg in l])
    mean = np.mean(cvg_array)
    total = np.sum([total for c, p, total, cvg in l])
    ct_str = ''
    if ct_list:
        for ct in ct_list:
            ct_x = np.sum(cvg_array > int(ct)) * 100.0 / len(cvg_array)
            ct_str = ct_str + '\t' + format(ct_x, '.1f')
    out_file.write(i.strip() + '\t' + str(total) + '\t' + format(mean, '.2f') + ct_str + '\n')

def run(depthofcoverage_file, intervals_file, output_file, ct_list, no_sort=False, no_header=False):

    # Sort intervals file
    print 'Sorting intervals file...'
    in_file = open(intervals_file, 'r')
    intervals = in_file.readlines()
    in_file.close()

    if not no_sort:
        intervals.sort(compare_intervals)
        intervals_sorted_file = open(intervals_file+'.sorted', 'w')
        for i in intervals:
            intervals_sorted_file.write(i)
        intervals_sorted_file.close()

    # Open DepthOfCoverage (DOC) file and other preparations
    in_file = open(depthofcoverage_file, 'r')
    out_file = open(output_file, 'w')
    ct_str = ''
    if ct_list:
        for ct in ct_list:
            ct_str = ct_str + '\t' + '%_above_' + ct
    out_file.write('Target\ttotal_coverage\taverage_coverage' + ct_str + '\n')
    drow = 'header'
    if not no_header:
        drow = in_file.readline()
    dchrom = ''
    dpos = 0
    l_loci = deque()
    icount = 0

    # Read in intervals
    print 'Processing intervals...'
    for i in intervals:
        ichrom, ipos = i.split(':')
        ispos, iepos = ipos.split('-')
        icount = icount + 1

        # Read DOC until we reach starting point of interval
        while (dchrom != ichrom or ispos != dpos) and drow and chroms[dchrom] < chroms[ichrom] or (chroms[dchrom] == chroms[ichrom] and int(dpos) < int(ispos)):
            #print 'Comparing dchrom=' + dchrom + ', dpos=' + str(dpos) + ' with interval=' + i + '. ' + str(cmp(chroms[ichrom], chroms[dchrom])) + ', ' + str(cmp(int(ispos), int(dpos)))
            drow = in_file.readline()
            if drow:
                drow = drow.split('\t')
                dchrom, dpos = drow[0].split(':')

        if icount % 1000 == 0:
            print 'Processing interval ' + icount + ' (' + i.rstrip() + ')'

        #Are start positions aligned?

        #Case 1: Yes.
        if (dchrom == ichrom and dpos == ispos):

            l_loci.append([dchrom, dpos, int(drow[1]), float(drow[2])])

            # Traverse DOC file until we reach end of interval, storing positions we read along the way
            while (dchrom != ichrom or ispos != dpos) and drow and chroms[dchrom] < chroms[ichrom] or (chroms[dchrom] == chroms[ichrom] and int(dpos) < int(iepos)):
                # print 'Comparing dchrom=' + dchrom + ', dpos=' + str(dpos) + ' with interval=' + i + '. ' + str(cmp(chroms[ichrom], chroms[dchrom])) + ', ' + str(cmp(int(ispos), int(dpos)))
                drow = in_file.readline()
                if drow:
                    drow = drow.split('\t')
                    dchrom, dpos = drow[0].split(':')
                    l_loci.append([dchrom, dpos, int(drow[1]), float(drow[2])])

        #Case 2: interval overlaps with previously read interval, so start position is off
        if (dchrom == ichrom and dpos > ispos) or (len(l_loci) > 0 and ((chroms[l_loci[0][0]] < chroms[ichrom]) or (chroms[l_loci[0][0]] == chroms[ichrom] and l_loci[0][1] < ispos))):

            # Retrieve earlier read DOC data, pop everything we don't need
            while len(l_loci) > 0 and ((chroms[l_loci[0][0]] < chroms[ichrom]) or (chroms[l_loci[0][0]] == chroms[ichrom] and l_loci[0][1] < ispos)):
                l_loci.popleft()

        #Are end positions aligned?

        #Case 1: No, interval is fully within the boundaries of a previous interval
        if (dchrom == ichrom and int(dpos) > int(iepos)):
            it = iter(l_loci)
            l_current = []
            for locus in it:
                if int(locus[1]) <= int(iepos):
                    l_current.append(locus)
            if len(l_current) > 0:
                output_final_values(i, l_current, ct_list, out_file)

        # Case 2: No, interval overlaps with previous interval but extends further than it
        # Traverse DOC file until we reach end of interval, storing positions we read along the way
        while (dchrom != ichrom or ispos != dpos) and drow and chroms[dchrom] < chroms[ichrom] or (chroms[dchrom] == chroms[ichrom] and int(dpos) < int(iepos)):
            drow = in_file.readline()
            if drow:
                drow = drow.split('\t')
                dchrom, dpos = drow[0].split(':')
                l_loci.append([dchrom, dpos, int(drow[1]), float(drow[2])])

        # Case 3: Yes.
        if (dchrom == ichrom and int(dpos) == int(iepos) and len(l_loci) > 0):
            output_final_values(i, l_loci, ct_list, out_file)

    print 'Done.'

def main():
    # command line arguments
    parser = argparse.ArgumentParser(
        description='Calculate depth of coverage for each interval, similar to GATK DepthOfCoverage, but without merging intervals',
        epilog='depth_of_coverage_per_interval version 1.0?1 (c)2015-2016 Aliz R. Rao all rights reserved')
    parser.add_argument('--depthofcoverage', '-d', required=True,
                        help='output of GATK -T DepthOfCoverage, containing coverage at each locus')
    parser.add_argument('--intervals', '-L', required=True,
                        help='file containing list of intervals. Sorted interval file will be created here ending in .sorted')
    parser.add_argument('--summaryCoverageThreshold', '-ct', required=False, action='append',
                        help='report the percentage of bases covered to an amount equal to or greater than this number (e.g. % bases >= CT)')
    parser.add_argument('--no_sort', action='store_true', required=False,
                        help='do not sort interval file (default:False)')
    parser.add_argument('--no_header', action='store_true', required=False,
                        help='depthofcoverage does not contain header line (default:False)')
    parser.add_argument('--output_file', '-o', required=True,
                        help='output file')

    args = parser.parse_args()

    check_file_exists(args.depthofcoverage)
    check_file_exists(args.intervals)

    run(depthofcoverage_file=args.depthofcoverage, intervals_file=args.intervals, output_file=args.output_file, ct_list=args.summaryCoverageThreshold, no_sort=args.no_sort, no_header=args.no_header)

    return 0


if __name__ == "__main__": sys.exit(main())
