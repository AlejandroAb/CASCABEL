 awk -F'\t' '{if($0 ~ "#" || NR == 1){print $0}else{for(i=2;i<=NF;i++){t[i-1]=t[i-1]+$i;}}}END{tt=t[1];for(i=2;i<=length(t);i++){tt=tt"\t"t[i]}print tt; }' $1
