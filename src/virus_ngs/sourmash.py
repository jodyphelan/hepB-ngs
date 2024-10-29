from os.path import expanduser
from virus_ngs import run_cmd, report
import os.path
import csv
import argparse


def sourmash_get_refrence(prefix: str, read1: str, read2: str):
    db = expanduser("~")+"/.virus-ngs/virus.sig"
    refdir = expanduser("~")+"/.virus-ngs/ref/"
    print(db)
    if read2:
        run_cmd(f"sourmash sketch dna -p abund,scaled=10 {read1}  {read2} --merge  {prefix}  -o {prefix}.sig" )
    else:
        run_cmd(f"sourmash sketch dna -p abund,scaled=10 {read1}  --merge  {prefix}  -o {prefix}.sig" )
    run_cmd(f"sourmash gather --threshold-bp 1000 {prefix}.sig {db} -o {prefix}.gather.csv" )
    if not os.path.isfile(f"{prefix}.gather.csv"):
        raise Exception("Can't find reference\n")
    rows = [row for row in csv.DictReader(open(f"{prefix}.gather.csv"))]
    if len(rows)==0:
        quit("Can't find reference\n")
    ref = rows[0]["name"].split(" ")[0]+".fasta"
    
    report["reference"] = ref
    
    ref = f"{refdir}/{ref}"
    return ref

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--read1', type=str, help='Forward read', required=True)
    parser.add_argument('-2', '--read2', type=str, help='Reverse read')
    parser.add_argument('-p', '--prefix', type=str, help='Prefix for output files', required=True)
    args = parser.parse_args()

    ref = sourmash_get_refrence(
        prefix=args.prefix,
        read1=args.read1,
        read2=args.read2
    )

    print(ref)