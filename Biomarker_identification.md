## Step 1
### Gene annotation and run blast
##### Requirements: ` Prodigal 2.6.3 `  ` blast 2.13.0 `

##### Aim: Identify CDS in the bacterial genomes.  

In **Bash**:
```
# go to the genome folder
for i in `ls *.fasta`;do;prodigal -i $i -a $i".faa";done

# combine cds into a single file
cat isolate*.faa > cds.faa

# make a blast database of cds from genomes of interest
makeblastdb -in cds.faa -dbtype prot -out cds

# perform blast search on the protein database
blastp -db cds -query nitrogen_fixation_biomarker.faa -outfmt 7 -out nf_blast.out
blastp -db cds -query acid_biomarker.faa -outfmt 7 -out acid_blast.out
```

## Step 2
### Identify candidiates for N fixation in an acidic field.  
In **R**:  

##### read in and arrange data
```R
setwd("~/Documents/technical_challenge")
nf <- read.table("nf_blast.out", sep = "\t")
ac <- read.table("acid_blast.out", sep = "\t")
cn <- c("query acc.ver", "subject acc.ver", "% identity", "alignment length", "mismatches", 
        "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")
colnames(nf) <- cn
colnames(ac) <- cn
```
##### Identify isolates with both biomarker

```R
fnf <- nf[which(nf$evalue <= 1.0e-10),]
fac <- ac[which(ac$evalue <= 1.0e-10),]
dt <- rbind(fnf, fac)
fnf <- fnf[,c(1,2)]
fac <- fac[,c(1,2)]
# split hits to get isolate # only
fnf$isolate <-data.frame(do.call('rbind', strsplit(as.character(fnf[,2]), '_', fixed=TRUE)))[,1]
fac$isolate <-data.frame(do.call('rbind', strsplit(as.character(fac[,2]), '_', fixed=TRUE)))[,1]  
m <- merge(fnf, fac, by = 3)
```

##### Write output files
```R
write.csv(m, "N fixation in an acidic field candidate.csv", row.names = F, quote = F)
write.csv(dt, "combined_filtered_biomarker.csv", row.names = F, quote = F)
nf_id <- nf$`subject acc.ver`
# output hits for nitrogen fixation marker for next step
write.table(nf_id, "nf_id.txt", row.names = F, quote = F, col.names = F) 
```


## Step 3
### Extect protein sequences that match the nitrogen fixation biomarker protein from each of the isolates that had matches.  
In **Jupyter Notebook** or **python**:
##### read in all cds records
```python
from Bio import SeqIO
records = SeqIO.parse("cds.faa", "fasta")
```
##### find records that are hits to nitrogen fixation marker
```python
with open('nf_id.txt') as file:
    hits = [line.rstrip() for line in file]
```
##### store the recrods in a list 
```python
hits_seq = []
for record in records:
    if record.id in hits:
        hits_seq.append(record)
```
##### output in fasta format for multiple sequence alignment
```python
with open("nf_hits.faa", "w") as output_handle:
    count = SeqIO.write(hits_seq, output_handle, "fasta")
```

## Step 4
### perform multiple sequence alignment of all nitrogen fixation biomarker homologs

You can use a list local tools such as clustal omega, MAFFT, etc. 
To save the effort of installing the tools, I used the online version of clustal omega  
###### Outputformat: STOCKHOLM  
###### Outputfile: nf_hits_msa.sto

## Step 5
### Build hmm profile and screen the unknown protein file for nitrogen fixation proteins
##### Requirements: ` HMMER 3.3.2 `
In **Bash**:
```
hmmbuild nf.hmm nf_hits_msa.sto  
hmmsearch nf.hmm unknown_proteins.faa >> nf_unknown_prot.out  
```

## Step 6
### Find the taxonomic identity of Isolate 3 which was identified as a good candidate for N fixation in an acidic field
You can use GTBD-tk to find the taxonomy of the isolate of interest. Previously, I used AMPHORA2 for taxonomy identification, too. This time, given that there is only one genome, I used a webserver for taxonomy identification: http://fbac.dmicrobe.cn/.


