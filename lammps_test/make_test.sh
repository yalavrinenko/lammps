#!/bin/sh
test_case (){
	sed "s/TEST_CASE/$1/g" lmp_template/eta.in > lmp.in;
	./l_opt -in lmp.in | grep "Step TotEng" -A1
	mpirun -n 4 ./l_opt -in lmp.in | grep "Step TotEng" -A1
	optirun ./e_opt $1.ini | grep "Total (Min)"
}

for ((i=1;i<1000;i++))
do
echo "*****************************************$i****************************************"
python ./create_configs.py $i
test_case $i
echo "***********************************************************************************"
done
