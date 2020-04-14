import argparse
import sys
from primer_processing import PrimerProcess
from context_manager import cd
from bam_processing import BamProcess
from core_calculation import Calculation
import os


primer_dir = "primer_align"
virus_dir = "virus_align"
cal_dir = "result"
offset = 10

def get_reference_name(ref_name):
    inf = open(ref_name, "r")
    line = inf.readline()
    buff = line.strip()
    ref = buff[1:]
    inf.close()
    return ref

def main():
    parser = argparse.ArgumentParser("reads_align2gene")
    parser.add_argument('-g', '--gff', help="GFF format gene annotation file for 2019-nCoV")
    parser.add_argument('-r', '--ref', help="Reference fasta for 2019-nCoV")
    parser.add_argument('-f', '--fastq', help="fastq format Nanopore sequenced reads")
    parser.add_argument('-p', '--primers', help="Primer tsv file, see example")
    parser.add_argument('-m', '--minimap2_path', help="Path for Minimap2, default is /usr/local/bin/minimap2")
    parser.add_argument('-b', '--blat_path', help="Path for Blat, default is /usr/local/bin/blat")
    parser.add_argument('-s', "--samtools_path", help ="Path for samtools, default is /usr/local/bin/samtools")
    parser.add_argument('-B', "--bcftools_path", help ="Path for bcftools, default is /usr/local/bin/bcftools")
    parser.add_argument('-o', '--out', help="output_prefix")
    options = parser.parse_args()

    primer_path = os.path.abspath(options.primers)
    ref_path = os.path.abspath(options.ref)
    fastq_path = os.path.abspath(options.fastq)
    gff_path = os.path.abspath(options.gff)
    ref_name = get_reference_name(ref_path)
    pcr_pools = set()
    primer_info = dict()
    if not os.path.exists(primer_dir):
        os.makedirs(primer_dir)
    with cd(primer_dir):
        p = PrimerProcess(options.blat_path, primer_path, ref_path, options.out)
        pcr_pools = set(p.pools)
        primer_info = dict(p.primer_info)
    if not os.path.exists(virus_dir):
        os.makedirs(virus_dir)
    with cd(virus_dir):
        BamProcess(options.minimap2_path, options.samtools_path, fastq_path, ref_path, options.out, pcr_pools, primer_info, offset)
    if not os.path.exists(cal_dir):
        os.makedirs(cal_dir)
    with cd(cal_dir):
        Calculation(gff_path, options.out, pcr_pools, offset, ref_name)
    return


if __name__ == "__main__":
    main()