#######  IRANGES  ######
library(IRanges)
library(GenomicRanges)
library(Biostrings)
## EXERCISE 1
# (a) Create an IRanges ir with 2 ranges
#   - the first called "range1" with start 1 and end 10
#   - the second  called "range2" with start 5 and end 14

ir <- IRanges(start = c(1,5) , end = c(10,14))
names(ir) <- paste("range",1:2,sep = "")
ir

# Show the width of the 2 ranges 
width(ir)

# Set a width of 8 only to the first range
width(ir[1]) <- 8
ir

# Create an IRanges ir2 with 1 range called "range3" with start 10 and end 19
ir2 <- IRanges(start = c(10), end = c(19))
names(ir2) <- paste("range", 3 , sep = "")
ir2

# Create an IRange ir3 that contains ir and ir2

ir3 <- c(ir,ir2) #concatenate function
ir3

# Create a "normal" IRange ir4 from ir3

ir4 <- reduce(ir3)
ir4

# Shift the range ir4 to the right of 5 positions

shift(ir4,5)

# Explain the result of these lines of code

flank(ir, 5, start = TRUE) #return the flanking region of the input Irange if start is set TRUE this means that we take the flaking (a fianco) region to the left of the input range
flank(ir, 5, start = FALSE) #we take the right flakin region (9-13, 15-19)
flank(ir, -5, start = TRUE) # in this case we have a negative width we indicate the overlap with the input Irange and if start is set as TRUE the overlap starts from the start (left) 

# Explain the result of this line of code

ix <- IRanges(start = c(10), end = c(100), names = c("range"))
ix
pix <- promoters(ix, downstream = 5, upstream = 5) 
pix 

#promoters function gives the promoter site of the IRanges passed as argument, downstream represents 
# the number of nucleotides on the 3' direction, while upstream represents the number of nucleotides on the 5'
# the full range in the promoters is defined: start(x)-upstream ; start(x) + downstream - 1

# Create a promoters ranges for ir3 using 2 nucleotides both up/down stream
pix2 <- promoters(ir3, upstream = 2, downstream = 2)
pix2

# Restrict the ir4 range between 5 and 10
restrict(ir4, start = 5, end = 10)
#--------------------------------------------------------------------------------

# EXERCISE 2
# Define the following ranges:
ir1 <- IRanges(start = c(1,21,26), end = c(20,25,30))
ir2 <- IRanges(start = c(15,5,22), end = c(20,30,28))

# Show how many overlaps are there between ir1 and ir2
#   - the ranges must have a match on the end and can 
#     have a max number of gaps equal to 2
ov <- findOverlaps(ir1,ir2)
ov 

# Show just the number of overlaps
countOverlaps(ir1,ir2) # 2 2 2 each range of ir1 has 2 overlaps
  
# Show which ranges of ir1 have overlaps
unique(queryHits(ov))

#--------------------------------------------------------------------------------

#EXERCISE 3 
#Create an IRange ir with starts (1,3,4) and width 2.

ir <- IRanges(start = c(1,3,4) , width=2)
ir

#(a) Reduce ir to its "normal" representation using reduce

reduce(ir) # reduce first orders ranges in x from left to right then merges the overlapping or adjacent ones


#(b) Use disjoin on ir. What is the difference with respect to reduce?

disjoin(ir) #it seems to be the opposite of the reduce function

#--------------------------------------------------------------------------------

#EXERCISE 4
#Create an IRange ir with starts (10,15,20) and width (2,5,7)
ir <- IRanges(start = c(10,15,20), width = c(2,5,7))
ir
#(a) Name your ranges in the following way: range1,range2,range3.
names(ir) <- paste("range", 1:3, sep = "")
ir
#(b) Resize ir by setting 3 as the new width and by keeping the end fixed.
resize(ir, width = 3, fix="end", use.names = TRUE)

#(c) Create an IRange ir2 with starts (1,3,10) and width 8.  
ir2 <- IRanges(start = c(1,3,10), width = 8)
ir2
#(d) Find the overlaps between ir1 and ir2 and save it in ov.
ov <- findOverlaps(ir1, ir2)
ov
#Extract the range corresponding to the second row of ov.
intersect(ir2[subjectHits(ov)[2]], ir1[queryHits(ov)[1]]) # 2? row of ov is between 1 ir1 & 2 ir2

#(e) Count the overlaps in two different ways.
countOverlaps(ir1, ir2)
unique(queryHits(ov))
#unique(subjectHits(ov))

#(f) Use the function gaps() on ir and explain what it does.
gaps(ir) #return the number of position between the ranges where there is a hole.

#--------------------------------------------------------------------------------


#######  GRANGES  ######

# EXERCISE 5 
# Q: what is the difference between IRanges and GRanges?

# GRanges are like IRanges with chromosomes and strands. Chromosomes in GRanges are called seqnames.

# Create a GRanges object with the following informations:
#   - 3 IRanges sequences with width equal to 10
#   - the sequences refers to the genome "h19"
#   - the sequences refers to 3 different non-circular chromosomes("chrx"),
#     each with length 100
#   - all the sequences are in the positive strand
#   - each sequence has a metadata information (ex: SYMBOL)

gr<-GRanges(seqnames=c("chrX","chrY","chrZ"),ranges=IRanges(start=c(9,22,10), width=10),strand=c("+","+","+"))
seqlengths(gr) <- c("chrX"=100, "chrY"=100, "chrZ"=100)
genome(gr) <- "h19"
values(gr) <- DataFrame(symbol = c("gene1","gene2","gene3"))
seqinfo(gr)
gr


#Change the chromosome names of the GRanges object according to the "dbSNP" database
seqlevelsStyle(gr) <- "dbSNP"
seqlevelsStyle(gr)
#--------------------------------------------------------------------------------

#EXERCISE 6
#Create a GRanges gr starting from ir and setting the chromosome to chr1. 
gr <- GRanges(seqnames = "chr1", ir)
gr
#(a) Set the strand to +. 
strand(gr) <- "+"
gr
#(b) Change the chromosome of the last sequence to chr8.
seqlevels(gr) <- c("chr1","chr8")
seqnames(gr) <- c("chr1", "chr1", "chr8")

#let expand GRanges with another chromosome
seqlevels(gr) <- c("chr1", "chr8", "chrX")


#(c) Now consider
gr2<-GRanges(seqnames=c("chr3","chr8","chr1"),ranges=IRanges(start=c(9,22,10),end=c(11,25,20)),strand=c("+","+","-"))
gr2
#Find the gaps of gr2 
gaps(gr2)
#(d) Find the overlaps between gr and gr2 ignoring the strand information. Explain the result.
ov <- findOverlaps(gr,gr2)
ov
#(e) Create gr2, which contains only the ranges of gr in chr8. (Use either dropSeqlevels or keepSeqLevels)
gr2 <- keepSeqlevels(gr, "chr1", pruning.mode = "coarse") #per? con chrom8
#(f) Map the chromosome names of gr to the ones used by Ensembl. Do it in a  single line of code.
gr <- renameSeqlevels(gr, mapSeqlevels(seqlevels(gr), "Ensembl"))
gr
#Is there any difference?

# Yes, the seqnames are not anymore chrom1 but they present the name 1!

#(g) Retrieve the coverage of gr and explain the result.
rl <- coverage(gr)
rl
#--------------------------------------------------------------------------------

#######  DNA STRING  ######

# EXERCISE 7
# Using the vector DNA_ALPHABET generate 2 random DNAString named "dna1" and "dna2"
#   - the length of the sequences must be 50

seq1 <- paste(sample(DNA_ALPHABET, size = 50, replace = TRUE), collapse = "")
dna1 <- DNAString(seq1)
dna1 <- paste(dna1, sep = " ")


seq2 <- paste(sample(DNA_ALPHABET, size = 50, replace = TRUE), collapse = "")
dna2 <- DNAString(seq2)
dna2 <- paste(dna2, sep = " ") 

# Create a DNAStringSet with the 2 previous sequences
sequence <- c(dna1, dna2)
setSequence <- DNAStringSet(sequence)
setSequence

# Do the following operations with the sequences:
#   - show the length of each sequence of the set
width(setSequence[1])
width(setSequence[2])

#   - sort the sequences in ascending order
sort(dna1, decreasing = FALSE)
sort(dna2, decreasing = FALSE)

#   - compute the reverse of the sequences
reverse(setSequence[1])
reverse(setSequence[2])

#   - compute the reverse complement of the sequences
reverseComplement(setSequence[1])
reverseComplement(setSequence[2])

#   - how is the frequency of "CGA" in the second sequence?
letterFrequency(setSequence[[2]], letters = 'CGA')
oligonucleotideFrequency(setSequence[[2]], width = 3)
#--------------------------------------------------------------------------------

# EXERCISE 8
#Create a DNAString d using sample and the IUPAC_CODE_MAP.
seq <- paste(sample(IUPAC_CODE_MAP, size=50, replace = TRUE), collapse = "")
mergeIUPACLetters(seq)
d <- DNAString(seq)
d <- paste(d, sep=" ")
set <- DNAStringSet(d)

#(a) What is its reverse complement?
reverse(set)
reverseComplement(set)

#(b) Concatenate d to a set of 3 strings (not created using DNAStringSet ). Is it possible?
s1 <- 'ACCT'
s2 <- 'GCC'
s3 <- 'CCAT'

s <- c(s1,s2,s3)
paste(d,s,collapse = "") # compatta tutto in un unica stringa

#Solve the "issue" in order to have a DNAStringSet called d_set
s <- append(d,s) # mantengo i miei oggetti separati 
d_set <- DNAStringSet(s)

#(c) Retrieve the number of sequences and their lengths.
length(d_set)
width(d_set[1])
width(d_set[2])
width(d_set[3])
width(d_set[4])
#(d) Retrieve the 3rd letter of the 3rd string
letter(d_set[[3]],3)
#(e) Show the number of occurrences of each nucleotide in d_set. 
letterFrequency(d_set, "A")
letterFrequency(d_set, "C")
letterFrequency(d_set, "G")
letterFrequency(d_set, "T")


#(f) Show the frequency of "GAT" in d_set
letterFrequency(d_set, "GAT")
#--------------------------------------------------------------------------------

########  RLE  ######

# EXERCISE 9
#Define r as the following vector rep(c(1:5),c(2:6)) and extract the Rle representation.
r <- rep(c(1:5), c(2:6)) # create an array from 1 to 5 with repetitions from 2 to 6
rl <- Rle(r) 
#(a) Explain what you have obtained.
#shows the values of the array and the number of times each value has a read

#(b) Extract: the run lengths, the run values and the number of runs in r.
runLength(rl)
runValue(rl)
sum(runLength(rl))
#(c) Use width(r). There's a difference with resect to runLength? 
width(rl) # there is no difference
#--------------------------------------------------------------------------------

# EXERCISE 10
# Create a GRanges object with sequence names (in order) chr1,chr2, chr1, strand all positive 
# and ranges with start 1,3,5 and width 10.
gr <- GRanges(seqnames = c("chr1", "chr2", "chr1"), strand = c("+", "+", "+"), ranges = IRanges(start = c(1,3,5), width = 10))
# Counts, for each integer, how many ranges overlap the integer (consecutive intervals with the same value) 
# by creating an Rle object.
rl <- coverage(gr)
rl
# Show the length of each runs and then the values.
runLength(rl)
runValue(rl)
# Create an IRanges object with start 2 and 3 and width 3.
ir <- IRanges(start = c(2,3) , width = 3)
ir
# Calculate the median only of chromosome 1 of the GRanges through a floating window consisting of the created IRanges. 
gr1 <- gr[which(seqnames(gr) == 'chr1')]
gr1 <- dropSeqlevels(gr1, 'chr2')
rl1 <- coverage(gr1)
aggregate(rl1, ir, FUN = median) ##??
