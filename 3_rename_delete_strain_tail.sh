for i in $(ls *.gff);do j=$(ls *.gff|sed -n "/$i/p"|cut -d "." -f 1|cut -d "_" -f 1,2);mv $i $j.gff;done
