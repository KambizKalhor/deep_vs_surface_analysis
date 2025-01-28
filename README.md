# deep_vs_surface_analysis

## Step 0: Download sample metagenomic data from the NCBI (paired-end Illumina reads). 

Download data from 2 samples:
```
wget
# OR
module load

```

Unzip the data
```
gunzip *.gz
```
do dump-spliiiiiit
```
ehtnsde
```


## Step 1: Run metaWRAP-Read_qc to trim the reads and remove human contamination
Note: you will need the bmtagger hg38 index to remove the human reads. You may also use another host genome to filter against with the `-x` option. Alternatively, use the `--skip-bmtagger` flag of of the ReadQC module to only do the read trimming.

```
mkdir READ_QC
metawrap read_qc -1 RAW_READS/forward_read_1.fastq -2 RAW_READS/reverse_read_2.fastq -t 24 -o READ_QC/output_directory

```

Lets have a glance at one of the output folders: `ls READ_QC/tsdhntdntdntd`

These are html reports of the read quality before and after QC:
```
post-QC_report
pre-QC_report
```

Original raw reads:
![Read quality before QC](https://i.imgur.com/x8aaFWs.png)
Final QC'ed reads:
![Read quality before QC](https://i.imgur.com/drJfxC9.png)


These are the final trimmed and de-contaminated reads:
```
final_pure_reads_1.fastq
final_pure_reads_2.fastq
```



## Step 2: Assembling the metagenomes with the metaWRAP-Assembly module


Assemble the reads with the metaSPAdes option flag (usually prefered over MegaHIT unless you have a very large data set):
```
metawrap assembly -1 final_pure_reads_1.fastq -2 final_pure_reads_1.fastq -m 200 -t 96 --use-metaspades -o ASSEMBLY
```


## Step 3: Run Kraken module on both reads and the assembly




## Step 4: Bin the co-assembly with three different algorithms with the Binning module

The initial binning process with CONCOCT, MaxBin, and metaBAT will be the more time intensive steps (especially CONCOCT and MaxBin), so I would advise you to run the Binning module with each of the algorithms seperately. However, metaWRAP supports running all three together. Our dataset is reasonably small, so I will run all three binning predictions at the same time.

If you are used to a different binning software(s), feel free to run them instead. The downstream refinement process (the Bin_refinement module) takes in up to 3 different bin sets, although you can get around this by running splitting your bin sets into groups and then recursively consolidating them.

Run the binning module with all three binners - notice how I put both the F and R read files at the end of the command.
```
metawrap binning -o INITIAL_BINNING -t 96 -a ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct CLEAN_READS/ERR*fastq
```

In the output folder, we see folders with the 3 final bin sets, and a file containing the average and standard deviation of the library insert sizes (this may or may not be usefull to you).  
```
insert_sizes.txt  concoct_bins	maxbin2_bins  metabat2_bins  work_files
```
Looking inside these folders reveals that we found 47, 29, and 20 bins with concoct, metabat2, and maxbin2, respectively. But we do not know how good these bins are yet (unless you used the `--run-checkm` flag). We will find out the quality of the bins in the next step!


## Step 5: Consolidate bin sets with the Bin_refinement module
Note: make sure you downloaded the CheckM database (see metaWRAP database instrucitons)

Now what you have metabat2, maxbin2, and concoct bins, lets consolidate them into a single, stronger bin set! If you used your own binning software, feel free to use any 3 bin sets. If you have more than 3, you can run them in groups. For example if you have 5 bin sets, try consolidating 1+2+3 and 4+5, and then consolidate again between the outputs.

When you do your refinement, make sure to put some thought into the minimum compleiton (-c) and maximum contamination (-x) parameters that you enter. During refinement, metaWRAP will have to chose the best version of each bin between 7 different versions of each bin. It will dynamically adjust to prioritize the bin quality that you desire. Consider this example: bin_123 comes in four versions in terms of completion/contamination: 95/15, 90/10, 80/5, 70/5. Which one is the best version? The high completion but high contamination, or the less complete but more pure bin? This is subjective and depends on what you value in a bin and on what your purposes for bin extraction are. 

By default, the minimum completion if 70%, and maximum contamination is 5%. However, because of the relatively poor depth of these demonstration samples, we will set minimum completion to 50% and maximum contamination to 10%, but feel free to be much more picky. Parameters like `-c 90 -x 5` are not unreasonable in some data (but you will get fewer bins, of course).

Run metaWRAP's Bin_refinement module:
```
metawrap bin_refinement -o BIN_REFINEMENT -t 96 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 50 -x 10
```

In the output directory, you will see the three original bin folders we fed in, as well as `metaWRAP` directory, which contains the final, consilidated bins. You will also see a `Binning_refiner` bin set - this is an internal benchmark for the bins produced by using another consolidation software (Binning_refiner). You can ignore this set - it will likely have low contamination and completion. You will also see .stats files for each one of the bin directories. 
```
concoct_bins.stats	maxbin2_bins.stats	metabat2_bins.stats	metawrap_50_10.stats	Binning_refiner.stats
concoct_bins		maxbin2_bins		metabat2_bins		metawrap_50_10_bins	Binning_refiner
```

Note that the '_50_10_' part of the naming refers to the `-c` and `-x` options you chose for the algorithm optimization. You can repeat the run with the same output directory using different options to get different results and see what works best on your sample (using '_90_5_' for example to get more near-complete genomes). The re-calculation will use the existing binning outputs, greatly reducing run time.

The .stat files contain usefull information about each bin, inscuding its completeness and contamination. For example, `cat BIN_REFINEMENT/metawrap_50_10.stats`:
```
bin	completeness	contamination	GC	lineage	N50	size	binner
bin.5	100.0	1.6	0.311	Euryarchaeota	12686	1705532	binsO.checkm
bin.4	99.32	1.342	0.408	Clostridiales	58825	2083650	binsO.checkm
bin.14	86.69	5.896	0.293	Bacteria	3754	2199676	binsO.checkm
bin.6	86.22	2.348	0.371	Clostridiales	4283	2055792	binsO.checkm
bin.8	83.16	2.516	0.446	Clostridiales	2723	1467846	binsO.checkm
bin.2	80.34	0.0	0.469	Bacteria	11936	3579466	binsO.checkm
bin.9	76.57	2.648	0.425	Selenomonadales	3155	1796524	binsO.checkm
bin.13	74.82	1.710	0.435	Bacteroidales	7456	3643185	binsO.checkm
bin.3	74.53	0.377	0.284	Clostridiales	10440	1241933	binsO.checkm
bin.10	65.78	0.0	0.263	Bacteria	3045	1159966	binsO.checkm
bin.11	64.85	3.776	0.417	Bacteroidales	2086	3103352	binsO.checkm
bin.1	57.36	0.0	0.430	Bacteria	4628	2673426	binsO.checkm
bin.7	52.94	1.724	0.501	Bacteria	3614	1465011	binsO.checkm
```

To evaluate how many "good bins" (based on out >50% comp., <10% cont. metric) metaWRAP produced, we can run
```
cat BIN_REFINEMENT/metawrap_50_10_bins.stats | awk '$2>50 && $3<10' | wc -l
13
```

By inspecting the other files, we find that metaBAT2, MaxBin2, CONCOCT, and metaWRAP produced 11, 7, 10, and 13 bins, respectively. So metaWRAP produced 2 more bins that the best single binner. Not bad! But this is just the number of bins. 

To closer compare the bin sets in terms of completion and contamination, we can look at the plots in `BIN_REFINEMENT/figures/`:
![Bin_refinement](https://i.imgur.com/m6RRJxi.jpg)

Note: This graph no longer has `Binning_refiner` in it, to reduce confusion. If you want to see Binning_refiner's performance, look at binsABC (or binsAB if you have two bin sets) in the other figure.

As you can see, the refinment process signifficantly produced the best bin set in terms of both compleiton and contamination. Keep in mind that these improvements are even more dramatic in more complex samples.


## Step 6: Visualize the community and the extracted bins with the Blobology module
Lets use the Blobology module to project the entire assembly onto a GC vs Abundance plane, and annote them with taxonomy and bin information. This will not only give us an idea of what these microbial communities are structured like, but will also show us our binning success in a more visual way. 

NOTE: In order to annotate the blobplot with bins with the `--bins` flag, you **MUST use the non-reassembled bins**! In other words, use the bins produced by the Bin_Refinment module, not the Reassemble_bins module. 

Note2: You will need R for this module.

```
metawrap blobology -a ASSEMBLY/final_assembly.fasta -t 96 -o BLOBOLOGY --bins BIN_REFINEMENT/metawrap_50_10_bins CLEAN_READS/ERR*fastq
```

You will find that the output has a number of blob plots of our communities, annotated with different levels of taxonomy or their bin membership. Note that to help with visualizing the bins, some of the plots only contain the contigs that were successfully binned (these files have `.binned.blobplot.` in their names). 

Phyla taxonomy of the entire assembled community:
![Phyla](https://i.imgur.com/VihLGWb.jpg)

Bin membership of all the contigs:
![bins](https://i.imgur.com/GDmIYe5.jpg)


## Step 7: Find the abundaces of the draft genomes (bins) across the samples
We would like to know how the extracted genomes are distributed across the samples, and in what abundances each bin is present in each sample. The Quant_bin module can give us this information. It used Salmon - a tool conventionally used for transcript quantitation - to estimate the abundance of each scaffold in each sample, and then computes the average bin abundances.

NOTE: In order to run this module, it is **highly** recomended to use the non-reassembled bins (the bins produced by the Bin_Refinment module, not the Reassemble_bins module) and provide the entire non-binned assembly with the `-a` option. This will give more accurate bin abundances that are in context of the entire community. 

Lets run the Quant_bins module:
```
metawrap quant_bins -b BIN_REFINEMENT/metawrap_50_10_bins -o QUANT_BINS -a ASSEMBLY/final_assembly.fasta CLEAN_READS/ERR*fastq
```

The output contains several usefull files. First, there is the `bin_abundance_heatmap.png` - a quick heatmap made to visualize the bin abundances across the samples. 
![heatmap](https://i.imgur.com/K1RaPUT.png)


The raw data for this plot (as you will most likely want to make your own heatmaps to analyze) is in `bin_abundance_table.tab`. Note that the abundances are expressed as "genome copies per million reads", and are calculated with Salmon in a simmilar way like TPM (transcripts per million) is calculated in RNAseq analysis. As such, they should be already standardized to the individual sample size. 

```
Genomic bins	ERR011349	ERR011348	ERR011347
bin.9	0.113912116828	35.851964987	39.8440491514
bin.10	0.273774684856	9.52869077293	39.988244574
bin.1	7.87827599808	31.3262582417	72.4475075589
bin.4	1.11852631889	100.052540293	111.213423224
bin.2	42.0242612674	69.0094806385	80.200001212
bin.5	2.16260151787	22.06396779	43.7720962538
bin.11	64.2884105466	25.3703846834	29.5444322752
bin.6	517.890689122	0.379711918465	0.834196723864
bin.7	0.499019812767	61.2121001057	82.5953338481
bin.14	5.49635966692	14.631433905	32.98399834
bin.13	0.230760165209	56.0018529273	91.6502833521
bin.8	98.4767064505	38.0691238971	22.8857472565
bin.3	349.730007621	0.0911113402849	0.196554603409
```

Finally, you can view the abundances of all the individual contigs in all the samples in the `quant_files` folder.

## Step 8: Re-assemble the consolidated bin set with the Reassemble_bins module
Now that we have our final, consilidated bin set in `BIN_REFINEMENT/metawrap_bins`, we can try to further improve it with reassembly. The Reassemble bins module will collect reads belonging to each bun, and then reassemble them sepperately with a "permissive" and a "strict" algorithm. Only the bins that improved through reassembly will be altered in the final set.

Let us run the Reassemble_bins module with all the reads we have:
```
metawrap reassemble_bins -o BIN_REASSEMBLY -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -t 96 -m 800 -c 50 -x 10 -b BIN_REFINEMENT/metawrap_50_10_bins
```

Looking at the output in `BIN_REASSEMBLY/reassembled_bins.stats`, we can see that 3 bins were improved though strict reassembly, 6 improved thorugh permissive reassembly, and 4 bins could not be improved (`.strict`, `.permissive`, and `.orig` bin extensions, respectively):
```
bin	completeness	contamination	GC	lineage	N50	size	binner
bin.10.orig	65.78	0.0	0.263	Bacteria	3045	1159966	NA
bin.7.strict	54.94	0.671	0.501	Clostridiales	3947	1474089	NA
bin.4.permissive	99.32	1.342	0.408	Clostridiales	72135	2088821	NA
bin.2.permissive	82.06	0.0	0.469	Bacteria	18989	3604843	NA
bin.14.strict	85.84	3.066	0.293	Bacteria	4576	2201824	NA
bin.9.permissive	76.74	2.554	0.425	Selenomonadales	3601	1802438	NA
bin.13.permissive	78.37	1.357	0.435	Bacteroidales	9887	3675176	NA
bin.11.orig	64.85	3.776	0.417	Bacteroidales	2086	3103352	NA
bin.6.permissive	88.04	1.006	0.371	Clostridiales	6288	2070146	NA
bin.5.orig	100.0	1.6	0.311	Euryarchaeota	12686	1705532	NA
bin.3.permissive	74.91	0.396	0.284	Clostridiales	16578	1243641	NA
bin.1.orig	57.36	0.0	0.430	Bacteria	4628	2673426	NA
bin.8.strict	83.89	1.342	0.446	Clostridiales	3870	1474833	NA
```

But how much did our bin set really improve? We can look at the `BIN_REASSEMBLY/reassembly_results.png` plot to compare the original and reassembled sets. We can see that the bin reassembly modestly improved the bin N50 and completion metrics, and signifficantly reduced contamination. Fantastic!
![heatmap](https://i.imgur.com/V8IosYQ.jpg)


We can also view the CheckM plot of the final bins in `BIN_REASSEMBLY/reassembled_bins.png`:
![heatmap](https://i.imgur.com/Yx00fuQ.png)



## Step 9: Determine the taxonomy of each bin with the Classify_bins module
Note: you will need the NCBI_nt and NCBI_tax databases for this module (see Database Installation section).

We already got an idea for the approximate taxonomy of each bin from CheckM in the `.stats` files in the Bin_refinment and Reassemble_bins modules. We can do better than that, however. The Classify_bins module uses Taxator-tk to accurately assign taxonomy to each contig, and then consolidates the results to estimate the taxonomy of the whole bin. Of course the success and accuracy of our predictions will rely heavily on the existing database.

Estimate the taxonomy of our final, reassembled bins with the Classify_bins module:
```
metawrap classify_bins -b BIN_REASSEMBLY/reassembled_bins -o BIN_CLASSIFICATION -t 48
```

We can view the final estimated taxonomy in `BIN_CLASSIFICATION/bin_taxonomy.tab`:
```
bin.1.orig.fa	Bacteria;Firmicutes;Clostridia;Clostridiales
bin.5.orig.fa	Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobrevibacter;Methanobrevibacter smithii
bin.11.orig.fa	Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides
bin.2.permissive.fa	uncultured organism
bin.10.orig.fa	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae
bin.14.strict.fa	Bacteria;Firmicutes
bin.8.strict.fa	Bacteria
bin.9.permissive.fa	Bacteria
bin.6.permissive.fa	Bacteria;Firmicutes;Clostridia;Clostridiales
bin.3.permissive.fa	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae
bin.4.permissive.fa	Bacteria
bin.13.permissive.fa	Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae
bin.7.strict.fa	Bacteria
```
As you can see some of the bins are annotated very deeply, while others can only be classified as "Bacteria". This method is relatively trustworthy, however it often fails to annotate organisms that are very distant from anything in the NCBI database. For these tricky bins, manually looking at marker genes (such as ribosomal proteins) can result in much more sensitive taxonomy assignment.


## Step 10: Functionally annotate bins with the Annotate_bins module
Now that we have our finalized reassembled bins, we are ready to functionally annotate them for downstream functional analysis. This module simply annotates the genes with PROKKA - it cannot do the actual functional analysis for you.

Run the functional annotation module on the final, reassembled bins: 
```
metaWRAP annotate_bins -o FUNCT_ANNOT -t 96 -b BIN_REASSEMBLY/reassembled_bins/
```

You will find the functional annotations of each bin in GFF format in the `FUNCT_ANNOT/bin_funct_annotations` folder
```
head FUNCT_ANNOT/bin_funct_annotations/bin.1.orig.gff
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	2866	3645	.	-	0	ID=HMOHEJHL_00001;inference=ab initio prediction:Prodigal:2.6;locus_tag=HMOHEJHL_00001;product=hypothetical protein
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	3642	4478	.	-	0	ID=HMOHEJHL_00002;inference=ab initio prediction:Prodigal:2.6;locus_tag=HMOHEJHL_00002;product=hypothetical protein
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	4606	5859	.	-	0	ID=HMOHEJHL_00003;inference=ab initio prediction:Prodigal:2.6;locus_tag=HMOHEJHL_00003;product=hypothetical protein
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	5856	6575	.	-	0	ID=HMOHEJHL_00004;Name=ypdB;gene=ypdB;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P0AE39;locus_tag=HMOHEJHL_00004;product=Transcriptional regulatory protein YpdB
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	6603	7658	.	-	0	ID=HMOHEJHL_00005;eC_number=1.1.1.261;Name=egsA;gene=egsA;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P94527;locus_tag=HMOHEJHL_00005;product=Glycerol-1-phosphate dehydrogenase [NAD(P)+]
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	7843	11946	.	-	0	ID=HMOHEJHL_00006;eC_number=3.2.1.78;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:A1A278;locus_tag=HMOHEJHL_00006;product=Mannan endo-1%2C4-beta-mannosidase
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	12597	13586	.	-	0	ID=HMOHEJHL_00007;eC_number=3.2.1.89;Name=ganB_1;gene=ganB_1;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:Q65CX5;locus_tag=HMOHEJHL_00007;product=Arabinogalactan endo-beta-1%2C4-galactanase
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	13719	14564	.	-	0	ID=HMOHEJHL_00008;Name=malG_1;gene=malG_1;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P68183;locus_tag=HMOHEJHL_00008;product=Maltose transport system permease protein MalG
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	14565	15476	.	-	0	ID=HMOHEJHL_00009;Name=malF_1;gene=malF_1;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P02916;locus_tag=HMOHEJHL_00009;product=Maltose transport system permease protein MalF
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	15484	16704	.	-	0	ID=HMOHEJHL_00010;Name=malE;gene=malE;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P0AEY0;locus_tag=HMOHEJHL_00010;product=Maltose-binding periplasmic protein
```

You will also find the translated and unstranslated predicted genes in fasta format in `FUNCT_ANNOT/bin_translated_genes` and `FUNCT_ANNOT/bin_untranslated_genes` folders. Finally, you can find the raw PROKKA output files in `FUNCT_ANNOT/prokka_out`.
