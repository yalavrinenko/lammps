#!/bin/sh
test_case (){
	sed "s/TEST_CASE/$2/g" lmp_template/eta.in > lmp.in;
	echo "AWPMD:" 
	$1 -in lmp.in | grep "Step TotEng" -A1
	mpirun -n 4 $1 -in lmp.in | grep "Step TotEng" -A1
	mpirun -n 8 $1 -in lmp.in | grep "Step TotEng" -A1
}

test_case_eff (){
	sed "s/TEST_CASE/$2/g" lmp_template/eta_hartree.in > hartree.in;
	echo "HARTREE:" 
	$1 -in hartree.in | grep "Step TotEng" -A1
	mpirun -n 4 $1 -in hartree.in | grep "Step TotEng" -A1
	mpirun -n 8 $1 -in hartree.in | grep "Step TotEng" -A1
}

test_case_eta (){
	echo "ETA:" 
	$1 inconfig/$2.ini | grep "Total (Min)"
}


for ((i=$2;i<$3;i++))
do
echo "*****************************************$i****************************************"
python ./create_configs.py $i
test_case $1 $i
test_case_eff $1 $i
test_case_eta ../energy_optimization $i
echo "***********************************************************************************"
done
