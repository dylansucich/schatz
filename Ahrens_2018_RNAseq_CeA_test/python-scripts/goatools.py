#!/usr/bin/env python
# Get http://geneontology.org/ontology/go-basic.obo
# Get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
# Data will be stored in this variable
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import goatools 
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GeneID2nt as GeneID2nt_mus
from goatools.go_enrichment import GOEnrichmentStudy


obo_fname = download_go_basic_obo()
gene2go = download_ncbi_associations()
obodag = GODag("go-basic.obo")
geneid2gos_mouse = read_ncbi_gene2go("gene2go", taxids=[10090])

geneid2symbol = {}

print("{N:,} annotated mouse genes".format(N=len(geneid2gos_mouse)))
print(GeneID2nt_mus.keys().head())

goeaobj = GOEnrichmentStudy(
        GeneID2nt_mus.keys(), # List of mouse protein-coding genes
        geneid2gos_mouse, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method