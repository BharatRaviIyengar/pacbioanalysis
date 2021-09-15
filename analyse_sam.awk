#!/usr/bin/gawk -f

function mul(x,y){
	if(x%y==0)
		return y
	else
		return x%y
}

function ceil(x){
	if(x==int(x))
		return x
	else
		return int(x)+1
}

BEGIN{
	PROCINFO["sorted_in"]="@ind_num_asc"
	for(i=33;i<=126;i++)
		phr33[sprintf("%c",i)]=i
	OFS="\t"
	rseq="ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGCGTATGGTCTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA"

	rbs="ATTTAAGAAGGAGATATACAT"
	rfull=rbs rseq

	for(i=1;i<=3;i++)
		rinit[i]=i+length(rbs);


	## GENETIC CODE ##
	while((getline < "/home/bharat/Documents/Lab_resources/Data/Selection/PacSeq/genetic_code.txt") > 0)
		prot[$1]=$3

}
FNR==1{
	fname=FILENAME
	sub(/_.*/,"",fname)
}

FNR>3{
	n=q=r=start=len=0
	n=split($6,aln,"[A-Z]|=",typ)
#	split($6,typ,"[0-9]+")

	for(i=1;i<=n-1;i++){
		if(typ[i]!="D")
			q+=aln[i]
		qpos[i]=q
		if(typ[i]=="=" || typ[i]=="X" || typ[i]=="D")
			r+=aln[i]
		rpos[i]=r+$4-1
	}
	if(r+$4-1!=length(rfull) || $4>length(rbs))
		next
	rec++

	split($1,rname,"/")
	
	# MAP through all "events" #
	for(i=1;i<=n-1;i++){
		prevcodon=-1
		# SUBSTITUTIONS #	
		if(typ[i]=="X" && rpos[i-1]>=length(rbs)){
			rabspos=rpos[i-1]-length(rbs)
			for(j=1;j<=aln[i];j++){
				allsnv++							# nreads SNPs in sample
				nuclvar=nuclvar rabspos+j substr($10,qpos[i-1]+j,1)		# Genotype
				snv[rabspos+j substr($10,qpos[i-1]+j,1)]++			# SNV
				rrf=mul(rabspos+j,3);						# Position within codon
				codonidx=ceil((rabspos+j)/3)					# Codon number
				qcodonx[codonidx][rrf]=substr($10,qpos[i-1]+j,1)	
			}
		}
	}

	if(nuclvar!=""){
		nnvar++
		allnuclvar[nuclvar]++
	}
	for(c in qcodonx){
		for(r=1;r<=3;r++){
			if(!(r in qcodonx[c])){
				qcodonx[c][r]=substr(rseq,r+(c-1)*3,1)
			}
		}
		qcodon=qcodonx[c][1] qcodonx[c][2] qcodonx[c][3]				# query codon
		rcodon=substr(rseq,1+(c-1)*3,3)							# reference codon
		if(prot[qcodon]!=prot[rcodon]){							# non-synonymous substitutions
			sav[c prot[qcodon]]++
			codv[c prot[qcodon]]=rcodon OFS qcodon 
			protvar=protvar c prot[qcodon]
		}
		if(prot[qcodon]==prot[rcodon]){							# synonymous substitutions
			syn[c qcodon]++
			codv[c qcodon]=rcodon OFS qcodon OFS prot[qcodon] OFS prot[rcodon]
		}
	}
	
	print protvar > fname_"haptype.txt"							# geno(haplo)types

	if(protvar!=""){									# protein variants
		npvar++
		allprotvar[protvar]++
	}
	
	nuclvar=protvar=""
	delete qcodonx

}

END{

	for(i in syn)
		print i, syn[i], 1000*syn[i]/rec, codv[i] > fname"_SYN.txt"
	for(i in sav)
		print i, sav[i], 1000*sav[i]/rec, codv[i] > fname"_SAV.txt"
	for(i in allnuclvar){
		print i, allnuclvar[i], 1000*allnuclvar[i]/rec > fname"_nuclvar.txt"
		nsum+=allnuclvar[i]
	}
	print "WT", rec-nsum, 1000*(rec-nsum)/rec > fname"_nuclvar.txt"
	for(i in allprotvar)
		print i, allprotvar[i], 1000*allprotvar[i]/rec > fname"_protvar.txt"
		psum+=allprotvar[i]
	}
	print "WT", rec-psum, 1000*(rec-psum)/rec > fname"_protvar.txt"
} 
