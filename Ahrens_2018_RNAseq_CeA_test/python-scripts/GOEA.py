#!/usr/bin/env python
# Get http://geneontology.org/ontology/go-basic.obo
# Get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
# Data will be stored in this variable
# GOEA.py 

import sys
sys.path.append("/Users/cmdb/miniconda3/lib/python3.6/site-packages/goatools-0.8.9-py_0/site-packages/goatools")
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
# from goatools.test_data.genes_NCBI_10090_ProteinCoding import GeneID2nt as GeneID2nt_mus
# import GeneID2nt as GeneID2nt_mus
from goatools.go_enrichment import GOEnrichmentStudy


# ! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)
# !
# !       1.  DB                       MGI
# !       2.  DB Object ID             MGI:xxxx
# !       3.  DB Object Symbol
# !       4.  Qualifier
# !       5.  GO ID                    GO:xxxx
# !       6.  DB:Reference(s)          MGI:MGI:xxxx|PMID:xxxx
# !       7.  Evidence Code            3-digit (not ECO:xxxx)
# !       8.  With (or)From            optional
# !       9.  Aspect (GO DAG Abbreviation (F, P, C))
# !       10. DB Object Name           optional
# !       11. DB Object Synonym(s)     optional
# !       12. DB Object Type
# !       13. Taxon                    taxon:10090
# !       14. Date                     YYYYMMDD
# !       15. Assigned By
# !       16. Annotation Extension     same as GPAD/col 11
# !       17. Gene Product Form ID     Isorform
# ['MGI', 'MGI:3612240', 'Lsmem2', '', 'GO:0008150', 'GO_REF:0000015|MGI:MGI:2156816', 'ND', '', 'P', 'leucine-rich single-pass membrane protein 2', 'BY080835', 'protein', 'taxon:10090', '20100202', 'MGI', '', '\n']
genes = []
gene_symbol = []
GeneOntologyID = []
association = {}



with open(sys.argv[1], "r") as filehandle:  
    for line in filehandle:
        currentPlace = line[:-1]
        genes.append(currentPlace)
        


for line in open( sys.argv[2] ):
    info = line.split("\t")
    gene_symbol.append(info[2])
    GeneOntologyID.append(info[4])
    association[info[2]] = info[4]
    continue

# for i in range(len(gene_symbol)):
#     print(gene_symbol[i]

for gene_symbol, GeneOntologyID in association.items():
    print(gene_symbol + "\t" + GeneOntologyID)

with open("population.out", 'w') as file:
    for gene_symbol, GeneOntologyID in association.items():
        file.write(str(gene_symbol) + "\t" + str(GeneOntologyID.upper()) + "\n")






obo_fname = download_go_basic_obo()

gene2go = download_ncbi_associations()

obodag = GODag("go-basic.obo")

geneid2gos_mouse = read_ncbi_gene2go("gene2go", taxids=[10090])


geneid2symbol = {}

goeaobj = GOEnrichmentStudy(
        genes, # List of mouse protein-coding genes
        geneid2gos_mouse, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method

# geneids_study = geneid2symbol.keys()
goea_results_all = goeaobj.run_study(genes)
goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]