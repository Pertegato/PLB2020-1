###### BSGENOME
## EXERCISE 1
# - Show which genomes are available in the BSgenome package
library(BSgenome)
head(available.genomes(),6)

# - Install and load the "BSgenome.Ecoli.NCBI.20080805" genome

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("BSgenome.Ecoli.NCBI.20080805")
BiocManager::install("BSgenome.Ecoli.NCBI.20080805")
library(BSgenome.Ecoli.NCBI.20080805)

# - Show informations about the genome
seqinfo(BSgenome.Ecoli.NCBI.20080805)

# - Show the names of the sequences of the genome
seqnames(BSgenome.Ecoli.NCBI.20080805)

# - Put the first and the second sequence in a DNAStringSet "set"
set <- DNAStringSet(list(BSgenome.Ecoli.NCBI.20080805[[1]], BSgenome.Ecoli.NCBI.20080805[[2]]))
# if in the constructor i pass the two elements as a vector their are concatenated.
set

# - Show a View with the first 100 elements of the first sequence of set
view1 <- set[[1]][1:100]

# - Show 100 random Views from the second sequence of set
#   - the width of each view must be 50

view2 <- replicate(100,sample(set[[2]],50,replace = TRUE))

# - Show the base frequencies of the previous views 
#letterFrequency(view1, 'A')
#letterFrequency(view1, 'C')
#letterFrequency(view1, 'G')
#letterFrequency(view1, 'T')
oligonucleotideFrequency(view1, 1)
#letterFrequency(view2[[1]], 'A')
#letterFrequency(view2[[1]], 'C')
#letterFrequency(view2[[1]], 'G')
#letterFrequency(view2[[1]], 'T')
oligonucleotideFrequency(view2[[1]], 1)


# - Show the letter frequency of "ACG" in the sequences of set
letterFrequency(set, letters = "ACG")
oligonucleotideFrequency(set, 3)

# - Show the letter frequency of "ACG" in Ecoli
param <- new("BSParams", X = Ecoli, FUN = letterFrequency)
head(bsapply(param, letters = "ACG"))

# - Search the sequence "ATGATAGATAG" in all sequences of Ecoli  
vmatchPattern("ATGATAGATAG", Ecoli)
# a function to find all the occurrences of a given pattern

## EXERCISE 2
# Exploit the package BSgenome to retrieve the mus musculus genome (latest release) from UCSC.
library('BSgenome.Mmusculus.UCSC.mm10')
Mmusculus # short name of the loaded genome
# (a)For the chromosome with maximum length retrieve the CG content percentage.
max_chrom <- names(seqlengths(Mmusculus)[1])
max_chrom
letterFrequency(Mmusculus$chr1, letters = "CG", as.prob = TRUE)

# (b) Do the same for all chromosomes using bsapply.
params <- new('BSParams', X = Mmusculus, FUN = letterFrequency)
bsapply(params, letters = 'CG', as.prob = TRUE)

# (c) Check if pattern d<-DNAString("AATC") occurs in a chromosome of your choice. 
d <- DNAString('AATC')
matchPattern(d, Mmusculus$chr1)

# (d) Count the occurences of the reverse of the pattern in all chromosomes and save it in n_oc
d <- reverse(d)
n_oc <- vcountPattern(d, Mmusculus)
# (e) Now consider the chromosome with the lowest number of occurrences (but greater than 0) and
# retrieve the sequence that starts at the start of the first occurence and which ends in the last position
# of the last occurrence.
n_ordered <- order(n_oc$count, decreasing = FALSE)
chrom_ordered <- n_oc[n_ordered,][[1,1]]
m <- matchPattern(d, Mmusculus$chr4_JH584295_random)
seq <- Mmusculus$chr4_JH584295_random[start(m[1]):end(m[4])]
# order the chromosomes with the respect to the count of occurrences

###### ANNOTATIONHUB 

## EXERCISE 1
# Using the AnnotationHub package retrieve a gtf file for Homo Sapiens
library(AnnotationHub)
ah <- AnnotationHub()
ah <- query(ah, c("Homo Sapiens", "GTF"))
#   - use the Ensembl database
#   - consider the GRCh37 genome (h19)
ah <- subset(ah, ah$dataprovider == "Ensembl")
ah <- subset(ah, ah$genome == "GRCh37")
ah <- subset(ah, ah$rdataclass == 'GRanges')
obj <- ah[['AH7558']]

#   - use the makeTxDbFromGRanges function to convert the GRanges 
txdb <- makeTxDbFromGRanges(obj)
#     obtained from AnnotationHub to a TxDb object
#   - show the following features from the TxDb object:
#       - genes
#       - exons
#       - transcripts
#   - Select a gene from the TxDb object and find its transcripts
#   - show the exons grouped by genes
genes(txdb)
exons(txdb)
transcripts(txdb)
select(txdb, keys = 'ENSG00000000003', column = c('TXID','TXNAME'), keytype = 'GENEID')
exonsBy(txdb, by = 'gene')
## EXERCISE 2
# Use the AnnotationHub to extract UCSC data that is from Homo sapiens and also specifically from the hg19 genome. What happens to the hub object as you filter data at each step?
ah1 <- AnnotationHub()
ah1 <- query(ah1, "UCSC")
length(ah1)
ah1 <- subset(ah1, ah1$genome == 'hg19')
length(ah1)
ah1 <- subset(ah1, ah1$species == 'Homo sapiens')
length(ah1)

## EXERCISE 3 
# Now that you have basically narrowed things down to the hg19 annotations from UCSC genome browser, lets get one of these annotations. Find the oreganno track and save it into a local variable.
#This track displays literature-curated regulatory regions, transcription factor binding sites, and regulatory polymorphisms from ORegAnno (Open Regulatory Annotation).
ah1 <- query(ah1, 'oreganno')
oreg <- ah1[['AH5087']]

## EXERCISE 4
# Load the org.Hs.eg.db package
library(org.Hs.eg.db)
# Look at the help page for the different columns and keytypes values with: help("SYMBOL").
# Now use this information to look up the entrez gene and chromosome for the gene symbol "MSX2" (from org.Hs.eg.db object).
keys <- 'MSX2'
columns <- c('ENTREZID','CHR')
select(org.Hs.eg.db, keys, columns, keytype = 'SYMBOL')
#selezionare l'entrzid e il cromosoma di MSX2 avendo come chiave il symbol (nome)

## EXERCISE 5
# Extract all of the gene symbols and their associated entrez gene ids from the org.Hs.eg.db package. 
# Then check the symbols for redundancy.
src <- org.Hs.eg.db
src_symbol <- keys(src, keytype = 'SYMBOL')
egr <- select(src, keys = src_symbol, 'ENTREZID', 'SYMBOL') # use this to get all the gene symbol that matched to all gene entrez id
# associare il nome dei geni al loro id, si ottengono due colonne di un dataframe
# Hint: use duplicated function to subset the redundant symbols
symbol <- egr$SYMBOL
red <- symbol[duplicated(symbol)]
select(src, keys = red, 'ENTREZID', 'SYMBOL') # only redundant symbols

## EXERCISE 6
# Use the accessors for the TxDb.Hsapiens.UCSC.hg19.knownGene package to retrieve the gene id,
# transcript name and transcript chromosome for all the transcripts. 
# Do this using both the select() method and also using the transcripts() method. 
# What is the difference in the output?
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
columns(TxDb.Hsapiens.UCSC.hg19.knownGene)
res1 <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
               keys(TxDb.Hsapiens.UCSC.hg19.knownGene, keytype="TXID"),
               columns=c("GENEID","TXNAME","TXCHROM"), keytype="TXID")

res2 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene,
                    columns = c("gene_id","tx_name"))
#Notice that in the 2nd case we don't have to ask for the chromosome, as transcripts() returns a GRanges object, so the chromosome will automatically be returned as part of the object.

## EXERCISE 7
# Use the src_organism object with the transcripts method to look up the entrez gene IDs for all gene symbols that 
# contain the letter 'X'.
# grange filter?
library(Organism.dplyr)
src <- src_organism(TxDb.Hsapiens.UCSC.hg19.knownGene)
xk = head(keys(src, keytype="entrez", pattern="X", column="symbol"))
select(src, xk, "symbol", "entrez") #check

## EXERCISE 8
# Pull down GO terms for entrez gene id "1" from human by using the ensembl "hsapiens_gene_ensembl" dataset.
library(biomaRt)
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
ids <- c('1')
select(org.Hs.eg.db, keys = ids, columns = 'GO', keytype = 'ENTREZID')