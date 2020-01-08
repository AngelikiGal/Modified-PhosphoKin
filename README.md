# Modified-PhosphoKin
A modification of the PhosphoKin tool (https://github.com/AngelikiGal/PhosphoKin) to identify deleterious missense cancer mutations that could obstruct the predicted phosphorylation of the sequence of laminin γ1-chain (LAMC1).

Authors: Panagiota-Aggeliki Galliou (ag.gal.work@gmail.com), Kleio-Maria Verrou (kleioverrou@yahoo.com)  

# Description

## This code:

### A) Does the same as PhosphoKin:
It searches kinase binding sites and potential phosphorylated residues in LAMC1. Also, it categorizes the phosphorylated residues according to their location relatively to the known active sites of LAMC1 (inside of, close to within six residues proximity and outside of active sites). Also, it links kinases to active sites according to the phosphorylated residues they are responsible for and shows their activity inside of, close to and outside of active sites. Also, it finds the possibly responsible kinases for the experimentally observed phosphorylated residues in LAMC1 according to the predicted phosphorylated residues of LAMC1.

### B) Identifies the deleterious missense cancer mutations in LAMC1 that could obstruct its phosphorylation:
It searches on a given .xlxs file the LAMC1 deleterious missense mutations. Then, it examines each LAMC1 deleterious missense mutation on whether it could obstruct the protein's predicted phosphorylation, according to kinases' binding motifs. A mutation is considered to obstruct LAMC1 phosphorylation when it occurs directly on predicted phosphorylated residue as well as when it occurs to the required for kinase binding surrounding residues of the predicted phosphorylated residue. Then, it writes on a new .xlx file those identified mutations in LAMC1 that could obstruct the protein's predicted phosphorylation, along with the mutations' information in the following order:

> ID of sample, in which the mutation was found
>
> Cancer type, in which the mutation was found
> location of the mutation in LAMC1
>
> Mutated residue
>
> Location of mutation, relatively to active sites (inside, close with six residues proximity, outside)
>
> Corresponding active site
>
> Number of kinases' phopshorylation the mutation could obstruct
>
> Name of possibly obstructed kinase
>
> Reason the mutation obstructed the kinase's motif (phosphorylated or required residue)

# Requirements

PhosphoKin is a tool written in Python 3 and, therefore, it can be run in most operating systems.

# How to run

## Needed files

PhosphoKin uses as input five .txt files; one containing LAMC1 in fast format, one containing the LAMC1 active sites, one containing the experimentally assigned phosphorylated residues in LAMC1 and one containing the kinases along with their recognition motifs and one all LAMC1 mutations found in cancer as downloaded by CBioPortal. The input files used are uploaded here: https://www.dropbox.com/sh/5p0b4fkri1uisbl/AABm9MqYvhaxEinaPKVTgfMGa?dl=0.

### 1) The Active Sites file:
The file must have the following format: Start_of_active_site_1 - End_of_active_site_1, Start_of_active_site_2 - End_of_active_site_2, etc.

e.g.

>13-22, 150-162, 1147-1458
### 2) The Experimentally Observed Phosphorylations file:
The file must have the following format: phosphorylated_residue_1, phosphorylated_residue_2, phosphorylated_residue_3, etc.

e.g.

>S4, S156, T445
### 3) The Motifs file:
The file must have the following format: Name_of_kinase_1[space]protein type: motif_1-motif_2-motif_3 Name_of_kinase_2[space]protein type: motif_1-motif_2-etc.

e.g.

>H1K Kinase:[pS/pT]P[R/K]-[pS/pT]PX[R/K]-[R/K][pS/pT]P
>
>ATM Kinase:pSQ-[P/L/I/M]X[L/I/E/D]pSQ-LpSQE
   
### 4) The Protein Sequence file:
The file must have a faste format (https://en.wikipedia.org/wiki/FASTA_format).

e.g.

>`>`sp|P11047|LAMC1_HUMAN Laminin subunit gamma-1 OS=Homo sapiens OX=9606 GN=LAMC1 PE=1 SV=3 MRGSHRAAPALRPRGRLWPVLAVLAAAAAAGCAQAAMDECTDEGGRPQRCMPEFVNAAFN VTVVATNTCGTPPEEYCVQTGVTGVTKSCHLCDAGQPHLQHGAAFLTDYNNQADTTWWQS QTMLAGVQYPSSINLTLHLGKAFDITYVRLKFHTSRPESFAIYKRTREDGPWIPYQYYSG SCENTYSKANRGFIRTGGDEQQALCTDEFSDISPLTGGNVAFSTLEGRPSAYNFDNSPVL

### 5) The Mutations file:
The file should have a specific format as downloaded from CBioPortal. The input file used can be used as a template and it is uploaded here: https://www.dropbox.com/s/e2paiz1jvhc8t07/LAMC1-Mutations.xlsx?dl=0.

# In Terminal

`python3 Modified-PhosphoKin-LAMC1.py`

Remember to have the required input files in the same directory as PhosphoKin or insert the whole path of the required input files.

# The output

As output the tool produces five .txt files:
### 1) One for the prediction of phosphorylation sites and phosphorylated residues in LAMC1:
LAMC1-Kinases_Motifs.txt

### 2) One for the identification of possibly responsible kinases for the experimentally observed phosphorylated residues in LAMC1:
LAMC1-Possily_responsible_kinases.txt

### 3) One for the categorization of phosphorylated residues according to their location relatively to active sites:
LAMC1-Categorization_of_phosphorylations_comparatively_to_active_sites.txt

### 4) One for the association of kinases with active sites:
LAMC1-Link_ActiveSite_with_kinases.txt

>_For the analysis of the output, please read the reference of Phosphokin:_
>
>_Galliou, P.A., Verrou, K.M., 2019. An in silico method for studying the phosphorylation in association to active sites. Aristotle Biomedical Journal 1, 48–59._

### 5) One for the LAMC1 deleterious missense mutations that could obstruct the predicted phosphorylation of LAMC1:
LAMC1-Mutations_disturbing_Phosphorylations.xls

# Got a Question?

Please do not hesitate to ask any question via email: Panagiota-Angeliki Galliou (ag.gal.work@gmail.com).
