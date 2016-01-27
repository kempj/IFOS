#lamboot -v lamhosts  
lamboot
mpirun -np 8 nice -19 ../bin/IFOS2D in_and_out/IFOS_INV.json | tee in_and_out/IFOS_INV.out

#cd su

#for ((i=1; i < ($1+1); i++)) ; do

#	cat IFOS_y.su.shot$i.* > IFOS_y.su.shot$i
#        cat IFOS_x.su.shot$i.* > IFOS_x.su.shot$i

#done

#cd ..			
