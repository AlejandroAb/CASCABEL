#!/bin/bash

abundance="rel"  #print relative abundance: "rel" 
                 #something different print absolute abundances
miss="T"         #print missing lines from first OTU table
                 #Set to "F" to avoid this


 awk -F'\t' -v rel_abun=${abundance} -v print_missing=${miss} 'function abs(v){return v < 0 ? -v : v} 
         FNR==NR{
                  if(NR==2){
                    samples="";
                    for(i=2;i<=NF;i++){
                      samp[i-1]=$i;  
                      samples=samples"\t"$i; 
                    }
                    print samples;
                  }else if(NR>2){
                    for(i=2;i<=NF;i++){
                       asv[$1,samp[i-1]]=$i;  
                       t[i-1]=t[i-1]+$i; 
                    }
                  }
                  next;
                }
                {
                  if(FNR>2){
                    mergedLine=$1;
                    split($0,abun,"\t");
                    taxa[$1]="OK";
                    for(i=2;i<=NF;i++){
                      if(asv[$1,samp[(i-1)]]>=0){
                        if(rel_abun == "rel"){
                          mergedLine=mergedLine "\t" (abs(asv[$1,samp[i-1]] - $i) / t[i-1]) * 100;
                        }else{
                          mergedLine=mergedLine "\t" abs(asv[$1,samp[i-1]] - $i); 
                        }
                      }else{
                        if(rel_abun == "rel"){ 
                          if(i==2){
                            mergedLine=$1"(NL)";
                          }
                          mergedLine=mergedLine "\t" (abun[i-2] / t[i-1]) * 100;
                        }else{
                          mergedLine=$0;
                        }
                      }
                    }
                    print mergedLine;
                  }
                } 
             END{
                 if(print_missing == "T"){
                   for(s in asv){
                     split(s,sep,SUBSEP);  
                     if(taxa[sep[1]]!="OK"){
                       taxa[sep[1]]="OK";
                       mergedLine=sep[1]"(NR)";
                       for(i=1;i<=length(samp);i++){
                         if(rel_abun == "rel"){
                           mergedLine=mergedLine"\t"(asv[sep[1],samp[i]]/t[i])*100; 
                         }else{
                           mergedLine=mergedLine"\t"(asv[sep[1],samp[i]]);
                         }
                       }
                       print mergedLine;
                     }                     
                   }
                  }
                 }'  $1 $2
