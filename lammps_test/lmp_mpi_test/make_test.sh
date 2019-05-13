#!/bin/sh
test_case (){
	sed "s/TEST_CASE/$1/g" lmp_template/eta.in > lmp.in;
	echo "AWPMD:" 
	../lmp -in lmp.in | grep "Step TotEng" -A1
	mpirun -n 4 ../lmp -in lmp.in | grep "Step TotEng" -A1
	mpirun -n 8 ../lmp -in lmp.in | grep "Step TotEng" -A1
}

test_case_eff (){
	sed "s/TEST_CASE/$1/g" lmp_template/eta_hartree.in > hartree.in;
	echo "HARTREE:" 
	../lmp -in hartree.in | grep "Step TotEng" -A1
	mpirun -n 4 ../lmp -in hartree.in | grep "Step TotEng" -A1
	mpirun -n 8 ../lmp -in hartree.in | grep "Step TotEng" -A1
}

test_case_eta (){
	echo "ETA:" 
	../energy_optimization inconfig/$1.ini | grep "Total (Min)"
}


for ((i=1;i<$1;i++))
do
echo "*****************************************$i****************************************"
python ./create_configs.py $i
test_case $i
test_case_eff $i
test_case_eta $i
echo "***********************************************************************************"
done
