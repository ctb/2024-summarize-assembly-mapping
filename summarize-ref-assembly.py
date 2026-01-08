#! /usr/bin/env python
# OG location: https://github.com/ctb/2024-summarize-assembly-mapping/blob/main/summarize-ref-assembly.py
import os.path
import json
import csv
import sys
import argparse

import pandas
import sourmash
from sourmash import sourmash_args

class MetagenomeInfo:
    headers = ["accession", "assembly_f_unweighted", "assembly_f_weighted",
               "assembly_f_readmapped", "assembly_f_readmapped_w",
               "ref_f_unweighted", "ref_f_weighted", "f_reads_mapped"]
               
    def __init__(self, metag_acc, *, ksize=31, grist_dir, assembly_dir):
        self.metag_acc = metag_acc
        self.ksize = ksize
        self.grist_dir = grist_dir
        self.assembly_dir = assembly_dir

        self.metag_sig_path = os.path.join(grist_dir, "sigs",
                                           f"{self.metag_acc}.trim.sig.zip")
        self.metag_sig = sourmash_args.load_one_signature(self.metag_sig_path, ksize=self.ksize)
    

    def calc(self):
        "Calculate all the numbers."
        print(f'for accession {self.metag_acc}:')
        self.calc_assembly_stuff(self.assembly_dir)
        self.calc_ref_based_kmer_stuff(self.grist_dir)
        self.calc_mapping_stuff(self.grist_dir)

    def get_row(self):
        """Return a CSV row, starting with metagenome accession and
        containing the rest of the attributes in self.headers."""
        xx = [self.metag_acc]
        for x in self.headers[1:]:
            val = getattr(self, x)
            xx.append(val)
        return xx

    def calc_assembly_stuff(self, assembly_dir):
        sigfile = os.path.join(assembly_dir,
                               f"{self.metag_acc}.megahit.fa.gz.sig")
        assert os.path.exists(sigfile), sigfile

        self.assembly_sig = sourmash_args.load_one_signature(sigfile,
                                                             ksize=self.ksize)

        # percent of flat k-mers accounted for by assembly
        self.assembly_f_unweighted = self.metag_sig.contained_by(self.assembly_sig)
        print(f"assembly/unweighted: {self.assembly_f_unweighted*100:.1f}% (assembly_f_unweighted)")

        # abundance weighted version:
        assembly_mh = self.assembly_sig.minhash.flatten()
        metag_mh = self.metag_sig.minhash
        intersect = assembly_mh.intersection(metag_mh.flatten()).inflate(metag_mh)

        # now sum:
        total_weighted_sum = metag_mh.sum_abundances
        intersect_weighted_sum = intersect.sum_abundances
        print(f"assembly/weighted: {intersect_weighted_sum / total_weighted_sum * 100:.1f}% (assembly_f_weighted)")
        self.assembly_f_weighted = intersect_weighted_sum / total_weighted_sum

        sig_mapped = os.path.join(assembly_dir,
                                  f'{self.metag_acc}.x.ma.fq.gz.sig')
        ma_sig = sourmash_args.load_one_signature(sig_mapped, ksize=self.ksize)
        self.assembly_f_readmapped = self.metag_sig.contained_by(ma_sig)
        print(f"% k-mers in reads mapped to assembly (flat): {self.assembly_f_readmapped*100:.1f}% (assembly_f_readmapped)")

        # calculate weighted version.
        self.assembly_f_readmapped_w = self.metag_sig.contained_by_weighted(ma_sig)
        print(f"% k-mers in reads mapped to assembly (weighted): {self.assembly_f_readmapped_w*100:.1f}% (assembly_f_readmapped_w)")

    def calc_ref_based_kmer_stuff(self, grist_dir):
        gather_csv = os.path.join(grist_dir, "gather", f"{self.metag_acc}.gather.csv.gz")
        df = pandas.read_csv(gather_csv)

        self.ref_f_unweighted = df['f_unique_to_query'].sum()
        self.ref_f_weighted = df['f_unique_weighted'].sum()
        
        print(f"total ref k-mers found (abund): {self.ref_f_weighted * 100:.1f} (ref_f_weighted)")
        print(f"total ref k-mers found (flat): {self.ref_f_unweighted * 100:.1f} (ref_f_unweighted)")

    def calc_mapping_stuff(self, grist_dir):
        leftover_csv = os.path.join(grist_dir, 'leftover',
                                    f"{self.metag_acc}.summary.csv")
        df = pandas.read_csv(leftover_csv)
        total_mapped_reads = df['n_mapped_reads'].sum()
        print(f"total mapped reads: {total_mapped_reads}")

        read_stats_file = os.path.join(grist_dir, 'trim', f"{self.metag_acc}.trim.json")
        with open(read_stats_file, 'rb') as fp:
            read_stats = json.load(fp)
        total_reads = read_stats['summary']['after_filtering']['total_reads']
        f_mapped = total_mapped_reads / total_reads
        print(f"fraction of reads that mapped: {f_mapped*100:.1f}% (f_mapped)")
        self.f_reads_mapped = f_mapped

def main(argv):
    p = argparse.ArgumentParser(argv)
    p.add_argument('accs', nargs='+')
    p.add_argument('-o', '--output-csv', help='output CSV here',
                   required=True)
    p.add_argument('-t', '--top-level-directory', required=True,
                   help="directory containing assembly/, atta/, and grist/outputs/, e.g. '/group/ctbrowngrp4/2025-zyzhao-assemloss/85/'")
    
    args = p.parse_args()

    tld = args.top_level_directory.rstrip('/') + '/'
    grist_dir = os.path.join(tld, 'grist/outputs')
    assembly_dir = os.path.join(tld, 'assembly/')

    # check all the paths
    for path in (tld, grist_dir, assembly_dir):
        assert os.path.isdir(path), f'{path} is not a directory or does not exist?'

    # gather info for each accession
    results = []
    for metag in args.accs:
        info = MetagenomeInfo(metag,
                              grist_dir=grist_dir,
                              assembly_dir=assembly_dir)
        info.calc()
        results.append(info.get_row())

    # write out CSV
    with open(args.output_csv, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(MetagenomeInfo.headers)
        for rr in results:
            w.writerow(rr)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
