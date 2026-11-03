import os
import subprocess
import warnings

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D

from motifcompass.sequenceConstruction.dna_regions import gff_to_fasta, read_fasta
from motifcompass.model.cross_validation import tune_sparse_rda
from motifcompass.featureConstruction.k_mer import generate_kmer_count_vector
from motifcompass.featureConstruction.pwm import pwm_logsumexp_matrix
from motifcompass.motif_de_novo.motif_construction import generate_PWMs_from_A, generate_meme_file


def generate_promoter(gff_file, fasta_file, output_fasta, upstream=1000, downstream=100):
    promoter_dict = gff_to_fasta(fasta_file, gff_file, output_fasta, upstream=upstream, downstream=downstream)
    return promoter_dict

def generate_X_pwm(promoter_dict, meme_file, pseudocount=0.01, background=(0.25, 0.25, 0.25, 0.25), both_strands=True):
    X_pwm = pwm_logsumexp_matrix(promoter_dict, meme_file, pseudocount=pseudocount, background=background, both_strands=both_strands)
    return X_pwm


if __name__ == "__main__":
    output_fasta = "data/promoter_features/promoters.fa"
    generate_promoter(
            gff_file="data/genomic_data/TAIR_12/annotation.gff",
            fasta_file="data/genomic_data/TAIR_12/col-0_keygene.fasta",
            output_fasta=output_fasta,
            upstream=500,
            downstream=200,
        )

    # genrate kmer matrix
    #X = generate_kmer_count_vector(
    #            promoters=read_fasta(output_fasta),
    #            k=6,
    #            both_strands=True,
    #        )
    #        # drop the coordinates (after ::)
    #X.index = [x.split("::")[0] for x in X.index]

    X = generate_X_pwm(
               promoter_dict=read_fasta(output_fasta),
               meme_file="data/binding/jaspar_plant_pfm.txt",
               pseudocount=0.01,
               background=(0.25, 0.25, 0.25, 0.25),
               both_strands=True,
           )