# MSA for phylogenetic inference
This is a workflow for reconstruction of phylogeny of ABCFs from _Firmicutes_.  For inferring trees I use `iqtree` and `RAxML`. To carry  bootstrap values over from  `RAxML` output to `iqtree` tree, OTU names must be identical as well as number of OTUs must be the same.  Requirements for sequence names and composition must be satisfied for both programs, so that they don't change OTU names or remove some sequences.  
_ OTU  - Operational Taxonomic Unit, it is basically the same sequence name_
`iqtree` - don't care about identical sequences but replaces some symbols  like `=` or `|` with underscore `_` form the flay. `RAxML` removes identical sequences as well as incompatible symbols to newik format. To keep output trees compatible these operations have to be done beforehand - that means checking must performed just before tree inferring, i.e. after trimming. That also means that in out case we can't ete3 workflow. Let  follow now step by step recipe.  

## Sequences
Increasing number of used sequences sets it's own requirements, if possible, add reference sequences. Reference sequences is the set of sequences with known class/sub-family. So, for final alignment we compose dataset from: a) study set; b) reference set
### Retrieve sequences (study set)
There are many ways retrieving desired data set. I keep my data in the table, including sequences. That table contains many additional columns what helps todo more operations than simply retrieving sequences. Sequences in this table   are found after HMM search using and refining results to identify most appropriate model for sequence. More than 4.600 _Firmicutes_ proteomes were scanned through (Exact procedure is described elsewhere). 

```bash
# get seq from table
table_to_fastA.py \
-i firmicutes_table_1e-70-80_with_seq_v3.csv \ # imput tabl 
-he 3 \ # Fasta heading  (0-4) bigger number, more detail header
-sl 420 \ # minimal sequence length 
> abcf_tbl_out.fa 
```
This script makes some cleanup in header name by removing or replacing symbols what are not accepted in newick format.

### Reduce redundancy
Removing identical and highly similar sequences can reduce dataset considerably leaving the whole variation in.
```bash
# make nr
usearch -cluster_fast abcf_tbl_out.fa -id 0.7  -centroids abcf_firmicutes_1e80_nr70_v3.fa
```

### Add reference set and clean names

```bash
# add reference
cat PoxtA_OpttA_victoriia.fa reference_ABCFs.fa   abcf_firmicutes_1e80_nr70_v3.fa > abcf_firmicutes_1e80_nr70_references_v3.fa
# remove funny characterhs from name
cat abcf_firmicutes_1e80_nr70_references_v3.fa  |  tr -d "=()'" | tr "[]:;,|.", "_"  > output.fa 
```
Sequence names in reference set might contain undesirable characters. Clean it again

### Multiple sequence alignment (MSA)
Is a part of ETE workflow 
`ete3 build -w mafft_linisi-trimal05-none-iqtree_bestmodel ...` 
but can be run as separate command. When running in server long time use `nohup`
```bash
# mafft-linsi
nohup mafft-linsi --reorder --quiet --thread 16 abcf_firmicutes_1e80_nr70_references_v3.fa  > abcf_firmicutes_1e80_nr70_refs_linsi_v3.fa  &
```

### Trimming alignment
Is a part of ETE workflow 
`ete3 build -w mafft_linisi-trimal05-none-iqtree_bestmodel ...` 
but can be run as separate command. 

```bash
# trimal 0.5
trimal -in abcf_firmicutes_1e80_nr70_refs_linsi_v3.fa -colnumbering  -gt 0.5 -fasta > abcf_firmicutes_1e80_nr70_refs_linsi_trimal05_v3.fa 
```
Output FastA file contains one line at the end what is not a part of FastA format. That needs to be deleted.

### Checking alignment
After trimming gapy regions are chopped and as result from some sequences identical part was left. To identify that use `raxml -check` option.
```bash
# raxml test
 raxml-ng --check --msa abcf_firmicutes_1e80_nr70_refs_linsi_trimal05_v3.fa --model LG+R10 --prefix T1
```
reveals 10 identical sequences and creates `T1.reduced.py` without these. Convert phylip to FastA
```bash
# get back to fastA
cat T1.raxml.reduced.phy | awk  'NR>2 {print ">"$0}' | sed 's/ /\n/' > abcf_firmicutes_1e80_nr70_refs_linsi_trimal05_reduced_v3.fa
# insenity check 
cat abcf_firmicutes_1e80_nr70_refs_linsi_trimal05_reduced_v3.fa |  tr -d "=()'" | tr "[]:;,|.", "_"  > output.fa  
```

### Check conserved  motifs of NBD
Common part for all ABCFs are two nucleotide binding domains - NBD1 and NBD2 (or ABC\_1 & ABC\_2) separated by linker region. NBD contains 5 conserved motifs what are forming nucleotide binding pocket: Wlaker-A (P-loop in 3D structure); Gln-loop; Signature; Walker-B; and His-loop. Walker-A, Gln-loop and Signature motifs are well conserved but Walker-B and His-loop might not so easy to recognise from alignment. The presence of first 3 domains in both NBDs was used as criteria to consider sequence as ABCF protein. After editing 2649 sequences left in the final alignment `abcf_firmi_2649_v4.fa`.

### Computing trees
We estimated the best fitting mode by   `iqtree` : __LG+R10__ 
Cipres portal was used for inferring phylogenetic history of ABCF family. Running time 72h:
1. iqtree - UFB (-bb 1100) 
2. iqtree - (-b 100)
3. raxml - (100 bootstraps)

	 
