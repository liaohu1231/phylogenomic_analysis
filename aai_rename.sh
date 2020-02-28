for i in $(less name.txt|cut -f 1);do j=$(less name.txt|sed -n "/$i/p"|cut -f 2);new=$i_$j;echo $i;sed -i 's/"$i"/"$new"/g' aai_colum3.txt;done
