#######  IRANGES  ######

## EXERCISE 1
# (a) Create an IRanges ir with 2 ranges
#   - the first called "range1" with start 1 and end 10
#   - the second  called "range2" with start 5 and end 14

# Show the width of the 2 ranges 

# Set a width of 8 only to the first range

# Create an IRanges ir2 with 1 range called "range3" with start 10 and end 19

# Create an IRange ir3 that contains ir and ir2

# Create a "normal" IRange ir4 from ir3

# Shift the range ir4 to the right of 5 positions

# Explain the result of these lines of code
flank(ir, 5, start = TRUE)
flank(ir, 5, start = FALSE)
flank(ir, -5, start = TRUE)

# Explain the result of this line of code
ix <- IRanges(start = c(10), end = c(100), names = c("range"))
pix <- promoters(ix, downstream = 5, upstream = 5)
pix

# Create a promoters ranges for ir3 using 2 nucleotides both up/down stream

# Restrict the ir4 range between 5 and 10

#--------------------------------------------------------------------------------

# EXERCISE 2
# Define the following ranges:
ir1 <- IRanges(start = c(1,21,26), end = c(20,25,30))
ir2 <- IRanges(start = c(15,5,22), end = c(20,30,28))

# Show how many overlaps are there between ir1 and ir2
#   - the ranges must have a match on the end and can 
#     have a max number of gaps equal to 2

# Show just the number of overlaps

# Show which ranges of ir1 have overlaps

#--------------------------------------------------------------------------------

#EXERCISE 3 
#Create an IRange ir with starts (1,3,4) and width 2.
#(a) Reduce ir to its "normal" representation using reduce
#(b) Use disjoin on ir. What is the difference with respect to reduce?

#--------------------------------------------------------------------------------

#EXERCISE 4
#Create an IRange ir with starts (10,15,20) and width (2,5,7)
#(a) Name your ranges in the following way: range1,range2,range3.
#(b) Resize ir by setting 3 as the new width and by keeping the end fixed. 
#(c) Create an IRange ir2 with starts (1,3,10) and width 8.  
#(d) Find the overlaps between ir1 and ir2 and save it in ov.
#Extract the range corresponding to the second row of ov.
#(e) Count the overlaps in two different ways.
#(f) Use the function gaps() on ir and explain what it does.

#--------------------------------------------------------------------------------


#######  GRANGES  ######

# EXERCISE 5 
# Q: what is the difference between IRanges and GRanges?
# Create a GRanges object with the following informations:
#   - 3 IRanges sequences with width equal to 10
#   - the sequences refers to the genome "h19"
#   - the sequences refers to 3 different non-circular chromosomes("chrx"),
#     each with length 100
#   - all the sequences are in the positive strand
#   - each sequence has a metadata information (ex: SYMBOL)

# Change the chromosome names of the GRanges object according to the "dbSNP" database

#--------------------------------------------------------------------------------

#EXERCISE 6
#Create a GRanges gr starting from ir and setting the chromosome to chr1. 
#(a) Set the strand to +. 
#(b) Change the chromosome of the last sequence to chr8.
#(c) Now consider
gr2<-GRanges(seqnames=c("chr3","chr8","chr1"),ranges=IRanges(start=c(9,22,10),end=c(11,25,20)),strand=c("+","+","-"))
#Find the gaps of gr2 
#(d) Find the overlaps between gr and gr2 ignoring the strand information. Explain the result.
#(e) Create gr2, which contains only the ranges of gr in chr8. (Use either dropSeqlevels or keepSeqLevels)
#(f) Map the chromosome names of gr to the ones used by Ensembl. Do it in a  single line of code. 
#Is there any difference?
#(g) Retrieve the coverage of gr and explain the result.

#--------------------------------------------------------------------------------

#######  DNA STRING  ######

# EXERCISE 7
# Using the vector DNA_ALPHABET generate 2 random DNAString named "dna1" and "dna2"
#   - the length of the sequences must be 50

# Create a DNAStringSet with the 2 previous sequences

# Do the following operations with the sequences:
#   - show the length of each sequence of the set
#   - sort the sequences in ascending order
#   - compute the reverse of the sequences
#   - compute the reverse complement of the sequences
#   - how is the frequency of "CGA" in the second sequence?

#--------------------------------------------------------------------------------

# EXERCISE 8
#Create a DNAString d using sample and the IUPAC_CODE_MAP.
#(a) What is its reverse complement?
#(b) Concatenate d to a set of 3 strings (not created using DNAStringSet ). Is it possible?
#Solve the "issue" in order to have a DNAStringSet called d_set.
#(c) Retrieve the number of sequences and their lengths. 
#(d) Retrieve the 3rd letter of the 3rd string
#(e) Show the number of occurrences of each nucleotide in d_set. 
#(f) Show the frequency of "GAT" in d_set

#--------------------------------------------------------------------------------

########  RLE  ######

# EXERCISE 9
#Define r as the following vector rep(c(1:5),c(2:6)) and extract the Rle representation.
#(a) Explain what you have obtained.
#(b) Extract: the run lengths, the run values and the number of runs in r.
#(c) Use width(r). There's a difference with resect to runLength? 

#--------------------------------------------------------------------------------

# EXERCISE 10
# Create a GRanges object with sequence names (in order) chr1,chr2, chr1, strand all positive 
# and ranges with start 1,3,5 and width 10.
# Counts, for each integer, how many ranges overlap the integer (consecutive intervals with the same value) 
# by creating an Rle object.
# Show the length of each runs and then the values.
# Create an IRanges object with start 2 and 3 and width 3.
# Calculate the median only of chromosome 1 of the GRanges through a floating window consisting of the created IRanges. 

