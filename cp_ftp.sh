#for i in $(less prokaryotes.csv|cut -d "," -f 15|cut -d "/" -f 14);do sed -n "/$i/p" prokaryotes.csv|cut -d "," -f 15|sed 's!$!/'$i'_genomic.fna.gz!g' >> ftp.txt;done

for i in $(less ftp.txt);do wget $i ;done
for i in $(ls *.gz);do gzip -d $i;done

