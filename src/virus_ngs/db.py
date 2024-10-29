#! /usr/bin/env python3
import os
import subprocess as sp
from tqdm import tqdm
import sys
from virus_ngs import run_cmd
import multiprocessing as mp
import argparse





patterns = {
    "Hepatitis B virus genotype A":"hepB A",
    "Hepatitis B virus genotype B":"hepB B",
    "Hepatitis B virus genotype C":"hepB C",
    "Hepatitis B virus genotype D":"hepB D",
    "Hepatitis B virus genotype E":"hepB E",
    "Hepatitis B virus genotype F":"hepB F",
    "Hepatitis B virus genotype G":"hepB G",
    "Hepatitis B virus genotype G":"hepB G",
    "Hepatitis B virus genotype I":"hepB I",
    "Hepatitis B virus genotype J":"hepB J",

}



taxid = {
    "hepB A": 489455,
    "hepB B": 489460,
    "hepB C": 2764122,
    "hepB D": 2847137,
    "hepB E": 2847138,
    "hepB F": 2847139,
    "hepB G": 2847140,
    "hepB H": 2847141,
    "hepB I": 2847142,
    "hepB J": 2847143,
}

def stream_fasta(f):
    seq = ""
    header = ""
    for l in open(f):
        if l.startswith(">"):
            if seq!="":
                yield(header,seq,genotype)
            header = l.strip().split()[0][1:]
            genotype = None
            for pattern in patterns:
                if pattern.lower() in l.lower():
                    genotype = patterns[pattern]
            seq = ""
        else:
            seq+=l.strip()
    yield(header,seq,genotype)




def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--accessions",required = True,type=str)
    parser.add_argument("--threads",default=mp.cpu_count()//4,type=int)
    parser.add_argument("--no-kraken",action="store_true")
    parser.add_argument("--add-human",action="store_true")
    parser.add_argument("--kmcp",action="store_true")

    args = parser.parse_args()

    args.accessions = os.path.abspath(args.accessions)

    data_dir = os.path.expanduser('~')+"/.virus-ngs/"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    os.chdir(data_dir)
    ref_dir = data_dir + "/ref/"
    kraken_ref_dir = data_dir + "/kraken_ref/"
    if not os.path.exists(ref_dir):
        os.makedirs(ref_dir)
    if not os.path.exists(kraken_ref_dir):
        os.makedirs(kraken_ref_dir)




    
    sys.stderr.write("Downloading reference files\n")
    run_cmd("datasets download virus genome accession --inputfile %s" % args.accessions)
    run_cmd("unzip -o ncbi_dataset.zip")

    sys.stderr.write("Processing reference files\n")
    id2tax = {}
    for name,seq,genotype in tqdm(stream_fasta('ncbi_dataset/data/genomic.fna')):
        print(name,genotype)
        if genotype is None:
            with open(ref_dir + name + ".fasta",'w') as O:
                O.write(f">{name}\n{seq}\n")
        else:
            id2tax[name] = taxid[genotype]
            with open(kraken_ref_dir + name + ".fasta",'w') as O:
                O.write(f">{name}|kraken:taxid|{taxid[genotype]}\n{seq}\n")
            with open(ref_dir + name + ".fasta",'w') as O:
                O.write(f">{name}\n{seq}\n")




    with open("taxid.map",'w') as O:
        for name in id2tax:
            O.write("%s\t%s\n" % (name,id2tax[name]))

    sys.stderr.write("Creating sourmash signature\n")
    run_cmd("sourmash sketch dna -p scaled=10 ~/.virus-ngs/ref/*.fasta --name-from-first -o virus.sig")
    # run_cmd("sourmash index hepB.sig ref/*.sig")

    if args.kmcp:
        sys.stderr.write("Creating kmcp database\n")
        run_cmd(f"kmcp compute -j {args.threads} --circular -k 31 -n 10 -l 150 -I ref -O refs.tmp --force")
        run_cmd(f"kmcp index -j {args.threads} -f 0.01 -I refs.tmp/ -O refs.kmcp --force")
        run_cmd("cp taxid.map refs.kmcp/")

    run_cmd("wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
    run_cmd("mkdir -p taxdump")
    run_cmd("tar -zxvf taxdump.tar.gz -C taxdump/")

    if not args.no_kraken:
        sys.stderr.write("Creating kraken2 database\n")
        run_cmd(f"kraken2-build --threads {args.threads}  --download-taxonomy --skip-maps --db kraken2 --use-ftp")
        if args.add_human:
            run_cmd(f"kraken2-build --threads {args.threads} --download-library human --db kraken2 --use-ftp")
        run_cmd("ls %s/ | parallel -j %s --bar kraken2-build --add-to-library %s/{} --db kraken2" % (kraken_ref_dir,args.threads,kraken_ref_dir))
        run_cmd(f"kraken2-build --build --db kraken2 --threads {args.threads}")
        run_cmd("kraken2-build --clean --db kraken2")

    print("\nDone!\n")
    quit(args.accessions)