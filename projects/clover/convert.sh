#!/bin/bash

date
obs=threept
obs=twopt


m_list=(0.00072   0.0148  0.0185  0.0222  0.1745  0.1962  0.2181  0.2399  0.2617  0.2884  0.3316  0.3814  0.4386  0.5044  0.5800  0.6670   )
th_list=(0.00)  #L16
#-0.5557 0.5557   

file_head_twist=1
file_head_nf=2
file_head_nsrc=0
file_head_l0=128
file_head_l1=64
file_head_l2=64
file_head_l3=64
file_head_nk=$((${#m_list[@]} ))
#?? pure negativi per fare r=-1
file_head_nmoms=${#th_list[@]}

file_head_beta=1.778
file_head_ksea=0.1394267
file_head_musea=0.00072
file_head_csw=1.69

sanfo_out="../../out_run3"
contractions=(V0P5 V1P5 V2P5 V3P5 P5P5) #follow the order
contractions=(P5P5) #follow the order
contractions=(P5S0 P5V1 P5V2 P5V3 P5V0 P5P5 P5A1 P5A2 P5A3 P5A0 P5T1 P5T2 P5T3 P5B1 P5B2 P5B3 S0P5 V1P5 V2P5 V3P5 V0P5 P5P5 A1P5 A2P5 A3P5 A0P5 T1P5 T2P5 T3P5 B1P5 B2P5 B3P5 S0S0 V0V0 A0A0 V1V1 V2V2 V3V3 A1A1 A2A2 A3A3 T1T1 T2T2 T3T3 V1T1 V2T2 V3T3 T1V1 T2V2 T3V3 B1B1 B2B2 B3B3 A1B1 A2B2 A3B3 B1A1 B2A2 B3A3 V0S0 S0V0)
contractions=(P5P5 V0P5) #follow the order

if [ $obs = "twopt" ]
then
size=$((  4*( file_head_nk*(file_head_nk+1)/2)*file_head_nmoms*file_head_nmoms*${#contractions[@]}*file_head_l0*2))
correlators=$((size/(file_head_l0*2)))
fi
if [ $obs = "threept" ]
then
size=$((file_head_nk*file_head_nk*file_head_nk*file_head_nmoms*file_head_nmoms*${#contractions[@]}*file_head_l0*2))
fi


write_header_ASCI ()
{
echo $file_head_twist 
echo $file_head_nf    
echo $file_head_nsrc  
echo $file_head_l0    
echo $file_head_l1    
echo $file_head_l2    
echo $file_head_l3    
echo $file_head_nk    
echo $file_head_nmoms 

echo $file_head_beta  
echo $file_head_ksea  
echo $file_head_musea 
echo $file_head_csw   



for((ik=0;ik<$((${#m_list[@]}));ik++))
do
    k=${m_list[${ik}]}  
    echo $file_head_ksea 
#    echo $file_head_ksea >> to_read
done
for((ik=0;ik<$((${#m_list[@]}));ik++))
do
    k=${m_list[${ik}]}  
    echo $k 
 #   echo -$k >> to_read
done


for((ith=0;ith<$((${file_head_nmoms}));ith++))
do
    th=${th_list[${ith}]}  
    mom=`echo $th/2. | bc -l | awk '{printf "%.5f", $1}'`
    echo  0.5  $mom  $mom  $mom  
done
echo  $size  
echo $Nconfs
}





#confs=`ls $sanfo_out | grep -v 0900 | grep -v 02900 `
#confs=`ls $sanfo_out`
#confs="01800 02200 02300 02500 02600 02700 02800 02900 03000 03100 03300 03400 03500 03600 03900 04000 04200 04300 04400 04600 05400 05600 05700 05800 07000 07200 07300"
#confs="00900  01100  01300  01500  01700  01900  02100  02300  02500  02700  02900  03100  03300  03500  03700  03900  04100  04300  04500  04700  04900  05100 01000  01200  01400  01600  02000  02200  02400  02600  02800  03000  03200  03400  03600  03800  04000  04200  04400  04600  04800  05000  05200"

#confs=`du -h -s $sanfo_out/* | grep -E  '58M' | awk '{print $2}'  `
confs=`ls $sanfo_out  | grep 0  `
confs=${confs//"$sanfo_out/"}
Nconfs=`ls $sanfo_out | grep 0 | wc -l `

echo  the configuration found are $Nconfs :
echo $confs
#echo total number of configuration  `ls $sanfo_out | grep -v 0900 | grep -v 02900 | wc -l `
#echo total number of configuration  `ls $sanfo_out |  wc -l `
echo ""


make
list=${contractions[0]}
for ((ic=1;ic<${#contractions[@]};ic++))
do
	list="$list|${contractions[$ic]}"
done

echo $list
#if [ $obs = "threept" ]
#then
write_header_ASCI > to_read_ll
write_header_ASCI > to_read_sl
write_header_ASCI > to_read_ls
write_header_ASCI > to_read_ss

for conf in $confs
do 
    printf '%s\t'   $conf  
    echo $conf  >> to_read_ll
    echo $conf  >> to_read_sl
    echo $conf  >> to_read_ls
    echo $conf  >> to_read_ss
#    for (( j=0;j<${correlators};j++ ))
#    do
	    cat  $sanfo_out/$conf/mes_contr_2pts_ll |   awk '/'${list}'/{ if( $2 != "P5P5") f='$((file_head_l0+1))'; else {i++; if (i%2==0) f='$((file_head_l0+1))';} } { if (f<'$((file_head_l0+1))' && f>0)  print ; f--; }' >> to_read_ll
	    cat  $sanfo_out/$conf/mes_contr_2pts_sl |   awk '/'${list}'/{ if( $2 != "P5P5") f='$((file_head_l0+1))'; else {i++; if (i%2==0) f='$((file_head_l0+1))';} } { if (f<'$((file_head_l0+1))' && f>0)  print ; f--; }' >> to_read_sl
	    cat  $sanfo_out/$conf/mes_contr_2pts_ls |   awk '/'${list}'/{ if( $2 != "P5P5") f='$((file_head_l0+1))'; else {i++; if (i%2==0) f='$((file_head_l0+1))';} } { if (f<'$((file_head_l0+1))' && f>0)  print ; f--; }' >> to_read_ls
	    cat  $sanfo_out/$conf/mes_contr_2pts_ss |   awk '/'${list}'/{ if( $2 != "P5P5") f='$((file_head_l0+1))'; else {i++; if (i%2==0) f='$((file_head_l0+1))';} } { if (f<'$((file_head_l0+1))' && f>0)  print ; f--; }' >> to_read_ss
#	    for ((ic=1;ic<=${#contractions[@]};ic++))
#    	    do
#	    if [ "${contraction[$ic]}" = "P5P5" ] ; then
		    
#		    cat  $sanfo_out/$conf/mes_contr_2pts_ll |  awk '/'${contractions[$ic]}'/{i++; if (i=='$((j*2+1))') f='$((file_head_l0 +1))';} f{print;f--;}' | tail -n +2>> to_read
#		    else
#		    cat  $sanfo_out/$conf/mes_contr_2pts_ll |  awk '/'${contractions[$ic]}'/{i++; if (i=='$((j+1))') f='$((file_head_l0 +1))';} f{print;f--;}' | tail -n +2 >> to_read
#		    fi
#	    done
#    done
done
date
echo converting in to binary
./convert  to_read_ll
./convert  to_read_sl
./convert  to_read_ls
./convert  to_read_ss
mkdir data
mv to_read*bin*  data/
rm to_read*
#mv to_read_bin.dat   meas_3pt.dat
#rm tmp to_read
#sbatch scr.sh

#fi
   
#if [ $obs = "twopt" ]
#then
#date
#for((ic=0;ic<${#contractions[@]};ic++))
#do
#   echo ${contractions[ic]} 
#   write_header_ASCI
#   for conf in $confs
#   do  
#       printf '%s\t'   $conf  
#       echo $conf  >> to_read
#       sed '/^$/d' $sanfo_out/$conf/mes_contr_2pts_SM  > tmp   #remove empty lines
#       cat tmp | grep -A$file_head_l0 ${contractions[ic]}  | sed '/--/d' | grep -v  ${contractions[ic]} > tmp1
#       cat tmp1 >> to_read
       
#   done
#   date
#   echo `wc -l to_read`
#   ./convert
#   mv to_read_bin.dat   meas_2pt_${contractions[ic]}.dat
#done
#rm tmp tmp1 to_read

#fi

date  
