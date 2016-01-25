rm -f su/IFOS* su/measured_data/IFOS*

lamboot
#lamboot -v lamhosts
mpirun -np 8 nice -19 ../bin/IFOS in_and_out/IFOS_FW.json | tee in_and_out/IFOS_FW.out
cd model
cp model_Test_rho_it_0.bin model_Test_vp_it_0.bin model_Test_vs_it_0.bin ../model_true/.
cd ..

cd su

for ((i=1; i < ($1+1); i++)) ; do

	mv IFOS_y.su.shot$i.it1 measured_data/IFOS_y.su.shot$i
	mv IFOS_x.su.shot$i.it1 measured_data/IFOS_x.su.shot$i

done

cd ..
