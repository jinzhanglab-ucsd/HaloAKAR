## plasmid from Seq Primer Rd1 to Sew Primer Rd2
## Left Forward: GT'CGCGCCACCCTGGTGGACGGCGGCACCGGCGGCAGCGAGCTCAGCNNKNNKNNKGAGACCTTCCAGGCCTTCCGCACC'ACCGACGTCGGCCGCAAGCTGATCATCGATCAGAACGTTTTTATCGAGGGTACGCTGCCGATGG; 
## Left Reverse: CC'ATCGGCAGCGTACCCTCGATAAAAACGTTCTGATCGATGATCAGCTTGCGGCCGACGTCGGTGGTGCGGAAGGCCTGGAAGGTCTCMNNMNNMNNGCTGAGCTCGCTGCCGCCGGTGCCGCCGTCCACCAGGGTGGCGCGAC

## Right Forward: CG'CGTCAAAGGTATTGCATTTATGGAGTTCATCCGCCCTATCCCGACCTGGGACGAATGGNNKNNKNNKAAGTTTTCTCA'AGAACAGATCGGCGAAAACATTGTGTGCAGGGTCATTTGTACCACGGGTCAAATTCC 
## Right Reverse: GG'AATTTGACCCGTGGTACAAATGACCCTGCACACAATGTTTTCGCCGATCTGTTCTTGAGAAAACTTMNNMNNMNNCC'ATTCGTCCCAGGTCGGGATAGGGCGGATGAACTCCATAAATGCAATACCTTTGACGCG

## Generate the comparison document with all the 21*21*21 options for codons (9261)
for i in Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Stop Thr Trp Tyr Val; do for j in Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Stop Thr Trp Tyr Val; do 
    for k in Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Stop Thr Trp Tyr Val; do
        echo -e $i'|'$j'|'$k'\t0' >> tmpCombIndex; done done done
mv tmpCombIndex CombIndex

##Right: 
##Change the filename and change the extension at the end fo line 17 to fastq.gz and double check if line 28 says million
for FL in PKA-right-treat-Med_S106_L002; do
  cutadapt -g CGTCAAAGGTATTGCATTTATGGAGTTCATCCGCCCTATCCCGACCTGGGACGAATGG -o ${FL}.s1.R1.fq.gz --discard-untrimmed ${FL}*R1*.fastq.gz
  cutadapt -a AAGTTTTCTCA -m 9 -M 9 -o ${FL}.s2.R1.fq.gz --discard-untrimmed ${FL}.s1.R1.fq.gz
  
  zcat ${FL}.s2.R1.fq.gz | awk 'NR%4==2' > ${FL}.singleDNA
  ##zcat ${FL}.s2.R2.fq.gz | awk 'NR%4==2' | tr ATCG TAGC | rev > ${FL}.singleDNA
  paste <(cut -c 1-3 ${FL}.singleDNA) <(cut -c 4-6 ${FL}.singleDNA) <(cut -c 7-9 ${FL}.singleDNA) | gzip -c > ${FL}.divideDNA.gz
  rm ${FL}.s[12].R[12].fq.gz && rm ${FL}.singleDNA
  zcat ${FL}.divideDNA.gz | sed -f DNA2Amino - | awk '{OFS="\t"; print $1 "|" $2 "|" $3}' | gzip -c > $FL.Amino.gz

  zcat $FL.Amino.gz | sort | uniq -c |  sed 's/^[ \t]*//g' | sed s/' '/'\t'/g | awk '{OFS="\t"; print $2, $1}' > $FL.CT
  awk -v g=$FL.CT 'BEGIN{OFS="\t"; while((getline < g)>0) {G[$1]=$2}; close(g)} {a=G[$1]>0?G[$1]:0; print $1,a}' CombIndex > $FL.COUNT && rm $FL.CT
  scale=`echo "scale=5;1000000/$(awk '{sum += $2} END {print sum}' $FL.COUNT)" | bc`; echo $scale
  awk '{OFS="\t"; print $1,$2*"'$scale'"}' $FL.COUNT > ${FL}_R1.CPM && rm $FL.COUNT
done