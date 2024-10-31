#! /usr/bin/env python3
import argparse
from uuid import uuid4
import csv
import os
from virus_ngs import *
# import pathogenprofiler as pp
# import shutil
# import yaml
import logging
from rich.logging import RichHandler
from .kraken2 import run_kraken_reads
from .sourmash import sourmash_get_refrence
from .protein import get_protein_variants

def main():

    parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-1','--read1',type=str,help='Forward read',required = True)
    parser.add_argument('-2','--read2',type=str,help='Reverse read')
    parser.add_argument('-p','--prefix',type=str,help='Prefix for output files',required = True)
    parser.add_argument('--platform',type=str,choices=['illumina','nanopore'],help='Sequencing platform',required=True)
    parser.add_argument('-t','--threads',type=int,help='Number of threads',default=4)
    parser.add_argument('--min-dp',type=int,default=50,help='Minimum depth for consensus')
    parser.add_argument('--reference-assignment-method',type=str,choices=['genotype','sourmash'],default='sourmash',help='Minimum depth for consensus')
    parser.add_argument('--kraken-db',type=str,help='Kraken2 database directory',default='kraken2')
    parser.add_argument('--fix-ref',help='Force a reference instead of building one')
    parser.add_argument('--consensus-variant-frequency',default=0.01,type=float,help='Minimum frequency of variant to be considered in final reference')
    parser.add_argument('--assemble',action="store_true",help='Try assembly')
    parser.add_argument('--conf',required=True,help='JSON file with conf')
    parser.add_argument('--debug',action="store_true",help='Debug mode')
    parser.add_argument('--keep-intermediate-files',action="store_true",help='Keep intermediate fastq files')
    parser.set_defaults(func=main)

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(
            level='DEBUG', format="%(message)s", datefmt="[%X]", handlers=[RichHandler()]
        )

    conf = json.load(open(args.conf))
    otu_conf = conf['otus']

    tmp = str(uuid4())
    args.data_dir = os.path.expanduser('~')+"/.virus-ngs/"
    args.kraken_db = args.kraken_db if os.path.isdir(args.kraken_db) else args.data_dir+"/kraken2"

    get_fastq_stats(read1=args.read1,read2=args.read2)
    
    top_otu, args.filtered_read1,args.filtered_read2 = run_kraken_reads(
        read1=args.read1,
        read2=args.read2,
        prefix=args.prefix,
        threads=args.threads,
        kraken_db=args.kraken_db,
        otu_conf=otu_conf
    )

    print(report)

    if args.reference_assignment_method=="sourmash":
    
        ref = sourmash_get_refrence(
            prefix=args.prefix,
            read1=args.filtered_read1,
            read2=args.filtered_read2
        )
    elif args.reference_assignment_method=="genotype":
        ref_seqs = {d['taxid']:d['reference'] for d in otu_conf}
        ref = get_ref_file(ref_seqs[top_otu])
    else:
        raise Exception("Invalid reference assignment method")
    logging.debug(ref)

    

    pilon_correct(
        ref=ref,
        r1=args.filtered_read1,
        r2=args.filtered_read2,
        platform=args.platform,
        consensus_name=args.prefix+".temp.fasta",
        bam_file=f"{args.prefix}.ref.bam",
        threads=args.threads,
        min_depth=args.min_dp
    )
    
    fasta_depth_mask(
        input=f"{args.prefix}.temp.fasta",
        output=f"{args.prefix}.temp.consensus.fasta",
        bam_file=f"{args.prefix}.ref.bam",
        depth_cutoff=args.min_dp,
        newchrom=args.prefix
    )
    os.remove(f"{args.prefix}.temp.fasta")

     

    sys.stderr.write("Consensus sequence generated\n")
    if args.platform=="illumina":
        run_cmd("bwa index %(prefix)s.temp.consensus.fasta" % vars(args))
        if args.read2:
            run_cmd("bwa mem -t %(threads)s -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:%(platform)s' %(prefix)s.temp.consensus.fasta %(filtered_read1)s %(filtered_read2)s | samtools sort -@ %(threads)s -o %(prefix)s.consensus.bam" % vars(args))
        else:
            run_cmd("bwa mem -t %(threads)s -R '@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:%(platform)s' %(prefix)s.temp.consensus.fasta %(filtered_read1)s | samtools sort -@ %(threads)s -o %(prefix)s.consensus.bam" % vars(args))
        remove_bwa_index(f"{args.prefix}.temp.consensus.fasta")
    else:
        run_cmd("minimap2 -ax map-ont -t %(threads)s %(prefix)s.temp.consensus.fasta %(filtered_read1)s | samtools sort -@ %(threads)s -o %(prefix)s.consensus.bam" % vars(args))
    run_cmd("samtools index %(prefix)s.consensus.bam" % vars(args))
    
    report['average_depth'] = get_average_depth(
        bam=args.prefix+".consensus.bam",
        ref=f"{args.prefix}.temp.consensus.fasta"
    )

    freebayes_correct(
        ref=f"{args.prefix}.temp.consensus.fasta",
        bam=args.prefix+".consensus.bam",
        platform=args.platform,
        prefix=args.prefix,
        output=f"{args.prefix}.final.consensus.fasta",
        threads=args.threads,
        min_depth=args.min_dp,
        min_freq=args.consensus_variant_frequency
    )


    print(otu_conf)
    tmp = {d['taxid']:d['coords'] for d in otu_conf if 'coords' in d}
    
    get_protein_variants(
        ref=ref,
        consensus=f"{args.prefix}.final.consensus.fasta",
        coords=tmp[top_otu],
        ref_aa=conf['ref_aa'],
        outfile=f"{args.prefix}.protein_variants.json"
    )

    json.dump(report,open(f"{args.prefix}.report.json",'w'))

    if not args.keep_intermediate_files:
        os.remove(args.filtered_read1)
        if args.read2:
            os.remove(args.filtered_read2)
        os.remove(f"{args.prefix}.koutput.txt")
        os.remove(f"{args.prefix}.ref.bam")
        os.remove(f"{args.prefix}.ref.bam.bai")



