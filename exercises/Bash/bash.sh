#Es 1
!/bin/bash
filename="Mycoplasma_genitalium_g37.ASM2732v1.37.gff3"
if [ -f "${filename}.gz" ]; then
rm ${filename}.gz
fi
wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_13_collection/mycoplasma_genitalium_g37/${filename}.gz
gunzip ${filename}.gz
mv $filename myco.gff3

head -n 20 myco.gff3

#Es2
grep -v -P "^#" myco.gff3 | cut -f1 | sort | uniq | wc -l

grep -v -P "^#" myco.gff3 | grep -P "\tgene\t" | wc -l

grep -v -P "^#" myco.gff3 | cut -f2 | sort | uniq

grep -v -P "^#" myco.gff3 | cut -f3 | sort | uniq -c

fields=`grep -v -P "^#" myco.gff3 | grep -P "\tgene\t" | cut -f9`
echo "$fields" | sed s/\;/\\\n/g | grep -P "^Name=" | sed s/Name=//g | sort |uniq | wc -l

grep -v -P "^#" myco.gff3 | grep -P "\tgene\t" | grep -v "Name=" | wc -l

#Es 3
!/bin/bash
rm Homo_sapiens.GRCh38.98.gff3.gz
wget ftp://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/Homo_sapiens.GRCh38.98.gff3.gz
gunzip Homo_sapiens.GRCh38.98.gff3.gz
mv Homo_sapiens.GRCh38.98.gff3 grch38.gff3

head -n 10 grch38.gff3

grep -v "#" grch38.gff3 | head -n 25

function stats {
ifile="$1"
echo -n "number of chomosomes: "
grep -v -P "^#" $ifile | grep -P "\tchromosome\t" | wc -l
grep -v -P "^#" $ifile | grep -P "\tchromosome\t"
echo ""
echo -n "number of genes: "
grep -v -P "^#" $ifile | grep -P "\tgene\t" | wc -l
echo ""
echo "sources:"
grep -v -P "^#" $ifile | cut -f2 | sort | uniq
echo ""
echo "feature types:"
grep -v -P "^#" $ifile | cut -f3 | sort | uniq -c
echo ""
echo -n "number of gene names: "
fields=`grep -v -P "^#" $ifile | grep -P "\tgene\t" | cut -f9`
echo "$fields" | sed s/\;/\\\n/g | grep -P "^Name=" | sed s/Name=//g | sort␣
, → | uniq | wc -l
echo ""
echo -n "number of genes without names: "
grep -v -P "^#" $ifile | grep -P "\tgene\t" | grep -v "Name=" | wc -l
}

stats "grch38.gff3"

tmp=`mktemp`
echo $tmp
grep -P "^1\t" grch38.gff3 > $tmp
head -n 10 $tmp
stats $tmp
rm $tmp
echo "done"

#Es 4
mkdir GRCH38
#grep -v -P "^#" grch38.gff3 | cut -f1 | sort | uniq
chrs=`grep -v -P "^#" grch38.gff3 | cut -f1 | sort | uniq | grep -P "^([1-9]+|X|Y|MT)$"`
for chr in $chrs
do
echo $chr
ofile="GRCH38/grch38.${chr}.gff3"
head -n 1 grch38.gff3 > $ofile
grep "##sequence-region
${chr} " grch38.gff3 >> $ofile
#cat $ofile
grep -P "^${chr}\t" grch38.gff3 >>$ofile
#head -n 10 $ofile
done

head -n 10 GRCH38/grch38.1.gff3

#Es 5
tmp=`mktemp`
chrs=`grep -v -P "^#" grch38.gff3 | cut -f1 | sort | uniq | grep -P "^([1-9]+|X|Y|MT)$"`
echo "# chr file nof_chrs nof_genes nof_sources not_ftypes nof_namesnof_nonames"
for chr in $chrs
do
ifile="GRCH38/grch38.${chr}.gff3"
#echo $ifile
stats $ifile > $tmp
#cat $tmp
c1=`grep "number of chomosomes: " $tmp | sed s/number\ of\ chomosomes:\ //g`
c2=`grep "number of genes: " $tmp | sed s/number\ of\ genes:\ //g`
start=`grep -n "sources:" $tmp | cut -d":" -f1`
end=`grep -n "feature types:" $tmp | cut -d":" -f1`
(( c3 = $end - $start - 2 ))
start=`grep -n "feature types:" $tmp | cut -d":" -f1`
end=`grep -n "number of gene names:" $tmp | cut -d":" -f1`
(( c4 = $end - $start - 2 ))
c5=`grep "number of gene names: " $tmp | sed s/number\ of\ gene\ names:\ //g`
c6=`grep "number of genes without names: " $tmp | sed s/number\ of\ genes\without\ names:\ //g`
echo "# $chr $ifile $c1 $c2 $c3 $c4 $c5 $c6"
done
rm $tmp

#Es 6
grep ">" data/families.fa | cut -f1 | sort | uniq | wc -l

grep ">" data/families.fa | cut -f2 | sort | uniq | wc -l

grep ">" data/families.fa | sed s/\>//g | cut -f1 | sort | uniq -c | head -n 20

grep ">" data/families.fa | cut -f3 | sort | uniq -c | head -n 20

tmp=`mktemp`
head -n 1000 data/families.fa | sed s/\>.*$/\>/g | tr -d '\n' | sed s/\>/\\n/g >$tmp
#sed s/\>.*$/\>/g data/families.fa | tr -d '\n' | sed s/\>/\\n/g >$tmp
function getlengths {
while read line
do
#echo "@" $line
echo $line | wc -c
done < $1
}
lengths=`getlengths $tmp`
n=0
avg=0
avgi=0
for l in $lengths
do
(( avgi = (($avgi * $n) + $l) / ($n + 1) ))
avg=`echo "(($avg * $n ) + $l ) / ($n + 1.0)" | bc -l`
#echo $l $avgi $avg
(( n = n + 1 ))
done
echo $n $avgi $avg
rm $tmp

tmp=`mktemp`
head -n 1000 data/families.fa | sed s/\>.*$/\>/g | tr -d '\n' | sed s/\>/\\n/g >$tmp
#sed s/\>.*$/\>/g data/families.fa | tr -d '\n' | sed s/\>/\\n/g >$tmp
function getlengths {
while read line
do
#echo "@" $line
echo $line | wc -c
done < $1
}
lengths=`getlengths $tmp | sort | uniq -c`
echo "$lengths"
rm $tmp

