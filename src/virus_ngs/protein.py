import sys
import argparse
import pysam
from Bio.Seq import Seq
from Bio import Align
from pydantic import BaseModel
import json
from typing import List
import logging



class Pos(BaseModel):
    start:int
    end:int

class Coords(BaseModel):
    P: List[Pos]
    S: List[Pos]
    X: List[Pos]
    C: List[Pos]




def get_protein_seqs(fa,chrom,c):
    res = {}

    res["P"] = str(Seq(fa.fetch(chrom,c.P[0].start-1,c.P[0].end) + fa.fetch(chrom,c.P[1].start-1,c.P[1].end)).translate())[:-1]
    res["S"] = str(Seq(fa.fetch(chrom,c.S[0].start-1,c.S[0].end) + fa.fetch(chrom,c.S[1].start-1,c.S[1].end)).translate())[:-1]
    res["X"] = str(Seq(fa.fetch(chrom,c.X[0].start-1,c.X[0].end)).translate())[:-1]
    res["C"] = str(Seq(fa.fetch(chrom,c.C[0].start-1,c.C[0].end)).translate())[:-1]

    return res

def load_fasta_file(filename: str) -> dict:
    """
    Load a fasta file and return a dictionary with the sequence
    """
    sequences = {}
    with open(filename) as f:
        for line in f:
            if line.startswith(">"):
                name = line.strip()[1:]
                sequences[name] = ""
            else:
                sequences[name] += line.strip()
    return sequences



def get_protein_changes(seq1, seq2):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    aligments = aligner.align(seq1, seq2)
    aln = sorted(aligments,key=lambda x: x.score)[0]
    aln1, aln2 = aln[0], aln[1]
    changes = []
    for i in range(len(aln1)):
        if aln1[i] != aln2[i]:
            changes.append((i+1, aln1[i], aln2[i]))
    return changes



def get_protein_variants(ref,consensus,coords,ref_aa,outfile):

    
    c = Coords(**coords)
    print(c)

    ref = pysam.FastaFile(ref)
    ref_chrom = ref.references[0]
    con = pysam.FastaFile(consensus)
    con_chrom = con.references[0]
    # ref_prot = get_protein_seqs(ref,ref_chrom)
    con_prot = get_protein_seqs(con,con_chrom,c)
    logging.debug(con_prot)
    ref_protein_seqs = ref_aa
    logging.debug(ref_protein_seqs)


    protein_changes = []
    for name in ref_protein_seqs:
        ref_seq = ref_protein_seqs[name]
        con_seq = con_prot[name]
        changes = get_protein_changes(ref_seq, con_seq)
        for change in changes:
            protein_changes.append({
                "gene": name,
                "pos": change[0],
                "ref": change[1],
                "alt": change[2]
            })

    json.dump(protein_changes, open(outfile, "w"))

def cli():
    parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref',type=str,help='Fasta file',required=True)
    parser.add_argument('--coords',type=str,help='Fasta file',required=True)
    parser.add_argument('--ref-aa',type=str,help='Fasta file',required=True)
    parser.add_argument('--consensus',type=str,help='Fasta file',required=True)
    parser.add_argument('--outfile',type=str,help='Fasta file',required=True)

    args = parser.parse_args()  
