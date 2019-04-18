#!/bin/sh
test_case (){
	sed "s/TEST_CASE/$1/g" lmp_template/eta.in > lmp.in;
	../lmp -in lmp.in | grep "Step TotEng" -A1
	mpirun -n 4 ../lmp -in lmp.in | grep "Step TotEng" -A1
	mpirun -n 8 ../lmp -in lmp.in | grep "Step TotEng" -A1
	../energy_optimization inconfig/$1.ini | grep "Total (Min)"
}

for ((i=1;i<$1;i++))
do
echo "*****************************************$i****************************************"
python ./create_configs.py $i
test_case $i
echo "***********************************************************************************"
done