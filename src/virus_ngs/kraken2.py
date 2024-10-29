import os
import subprocess as sp
from tqdm import tqdm
import json
from virus_ngs import run_cmd
from os.path import expanduser
import pysam
from typing import List, Tuple
import argparse
import logging

class TaxonTree:
    """
    A data structure to store a taxonomic tree.
    
    Attributes
    ----------
    tree : dict
        A dictionary representing the tree. Keys are parent nodes and values are lists of child nodes.
    """
    def __init__(self,node_file='/Users/jody/.taxonkit/nodes.dmp'):
        """
        Parameters
        ----------
        node_file : str
            Path to the NCBI nodes.dmp file
        """
        self.tree = {}
        for l in open(node_file):
            row = [f.strip() for f in l.strip().split('|')]
            if row[0]==row[1]:
                continue
            child = int(row[0])
            parent = int(row[1])
            if parent not in self.tree:
                self.tree[parent] = []
            self.tree[parent].append(child)
    def find_descendants(self,node):
        """
        Recursively find all descendants of a given node.
        
        Parameters
        ----------
        node : int
            The node to find descendants of.
        
        Returns
        -------
        descendants : set
            A set of all descendants of the given node.
        """
        descendants = []
        children = self.tree.get(node, [])
        for child in children:
            descendants.append(child)
            # Recursively find descendants of the child
            descendants.extend(self.find_descendants(child))
        return set(descendants)
    


def kreport_extract_human(kreport_file: str):
    """
    Extract the human reads from a kraken report file.
    """
    human_reads = 0
    for line in open(kreport_file):
        if line.strip().split("\t")[4] == "9606":
            human_reads = float(line.strip().split("\t")[0])
    return {"Read percent human":human_reads}

def kreport_extract_otus(kreport_file:str, otus: list):
    """
    Extract the hepB reads from a kraken report file.
    """
    results = {int(otu):0 for otu in otus}
    for line in open(kreport_file):
        otu = int(line.strip().split("\t")[4])
        if otu in otus:
            results[otu] = float(line.strip().split("\t")[0])
    return results
    
def filter_fastq_by_taxon(kraken_output: str, otus: List[int], parent_otu: int, keep_otu: int,reads: str,output: str):
    """
    Filter a fastq file by taxon.
    """
    nodedump = expanduser("~")+"/.virus-ngs/taxdump/nodes.dmp"
    tree = TaxonTree(nodedump)

    include = [parent_otu]
    exclude = otus
    exclude.remove(keep_otu)

    taxon_to_keep = set()
    for node in include:
        taxon_to_keep.add(int(node))
        taxon_to_keep.update(tree.find_descendants(int(node)))

    if exclude:
        for node in exclude:
            if int(node) in taxon_to_keep:
                taxon_to_keep.remove(int(node))
                taxon_to_keep.difference_update(tree.find_descendants(int(node)))

    read_names = set()
    for i,l in tqdm(enumerate(open(kraken_output))):
        row = l.strip().split('\t')
        if int(row[2]) in taxon_to_keep:
            read_names.add(row[1])

    percentage_kept = round(len(read_names)/i*100)
    logging.info(f"Keeping {percentage_kept}% from {len(read_names)} reads")

    with pysam.FastxFile(reads) as fin, open(output, mode='w') as fout:
        for entry in tqdm(fin,total=i):
            if entry.name in read_names:
                fout.write(str(entry) + '\n')

def run_kraken_reads(read1: str, read2: str, prefix: str, threads: int, kraken_db: str, otu_conf: str) -> Tuple[int, str,str]:



    if read2:
        run_cmd(f"kraken2 --db {kraken_db} --report {prefix}.kreport.txt --output {prefix}.koutput.txt --paired --threads {threads} {read1} {read2}")
    else:
        run_cmd(f"kraken2 --db {kraken_db} --report {prefix}.kreport.txt --output {prefix}.koutput.txt --threads {threads} {read1}")

    
    otus_designation = {int(d['taxid']):d['name'] for d in otu_conf}
    parent_otu = {v:k for k,v in otus_designation.items()}['parent']
    otus_designation.pop(parent_otu)
    
    otus = list(otus_designation)

    kraken_results = kreport_extract_otus(f"{prefix}.kreport.txt", otus_designation)
    
    major_otu = sorted([(x,kraken_results[x]) for x in otus],key=lambda x:x[1],reverse=True)[0][0]
    logging.info(f"Major OTU: {major_otu}, {otus_designation[major_otu]}")
    filter_fastq_by_taxon(
        f"{prefix}.koutput.txt", 
        otus, 
        parent_otu, 
        major_otu, 
        read1, 
        f"{prefix}.kraken_filtered.1.fq"
    )
    if read2:
        filter_fastq_by_taxon(
            f"{prefix}.koutput.txt", 
            otus, 
            parent_otu, 
            major_otu, 
            read1, 
            f"{prefix}.kraken_filtered.1.fq"
        )
        return major_otu, f"{prefix}.kraken_filtered.1.fq", f"{prefix}.kraken_filtered.2.fq"

    return major_otu, f"{prefix}.kraken_filtered.1.fq", None

def cli():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-1', '--read1', type=str, help='Forward read', required=True)
    argparser.add_argument('-2', '--read2', type=str, help='Reverse read')
    argparser.add_argument('-p', '--prefix', type=str, help='Prefix for output files', required=True)
    argparser.add_argument('-t', '--threads', type=int, help='Number of threads', default=4)
    argparser.add_argument('-db', '--kraken_db', type=str, help='Kraken2 database directory',  default=expanduser("~")+"/.virus-ngs/kraken2")
    argparser.add_argument('--otus', type=str, help='JSON file with OTU structure',  required=True)

    args = argparser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(message)s", datefmt="[%X]"
    )

    run_kraken_reads(args.read1, args.read2, args.prefix, args.threads, args.kraken_db, args.otus)


    