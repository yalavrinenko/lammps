#!/bin/sh

No_inter=$($1 -in none.in | grep "Begin interaction-End interaction" | awk '{print $4}')
II_inter=$($1 -in ii.in | grep "Begin interaction-End interaction" | awk '{print $4}')
EI_inter=$($1 -in ei.in | grep "Begin interaction-End interaction" | awk '{print $4}')
EE_inter=$($1 -in ee.in | grep "Begin interaction-End interaction" | awk '{print $4}')
All_inter=$($1 -in all.in | grep "Begin interaction-End interaction" | awk '{print $4}')
echo $No_inter "|" $II_inter "|" $EI_inter "|" $EE_inter "|" $All_inter