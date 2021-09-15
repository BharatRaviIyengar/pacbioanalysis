#!/usr/bin/mawk -f

function lena(arr, j, x){
	for(j in arr)
		x++
	return x
}

function dist(hap1,hap2, u, i, j, h1, h2){
	gsub(/[A-Z]/,"& ",hap1)
	gsub(/[A-Z]/,"& ",hap2)
	n1=split(hap1,h1," ")
	n2=split(hap2,h2," ")
	for(j in h1)
		u[h1[j]]
	for(j in h2){
		if(h2[j] in u)
			i[h2[j]]
		u[h2[j]]
	}
	return(lena(u)-lena(i))
}

BEGIN{
	FS="\t"
	OFS=","
	PROCINFO["sorted_in"]="@ind_num_asc"
	d=0
}

function wt(h){
	if(h=="")
		return "WT"
	else
		return h
}

FNR==1{
	s[$0]
	next
}

{
	for(j in s){
		d+=dist($0,j)
	}
	s[$0]
}

END{
	avd = 2*d/(FNR*(FNR-1))
	printf("%s\t%0.3f",substr(FILENAME,1,3),avd)
	printf("\n")
}
