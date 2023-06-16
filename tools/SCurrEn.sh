for i in $(seq 1 9) 
do 
  echo $i
    echo "Final energy from 1st step"
    energy=$(grep "e  e" pos/POSCAR_00$i/OUTCAR | tail -1)
    energyp=$(echo $energy | sed -e "s/.* = //g")
    echo $energyp >> currens.txt
    #echo "Current energy from 2nd step"
    #grep "e  e" pos/POSCAR_00$i2/RELAX/OUTCAR | tail -1
    vectora=$(sed -n 3"p" <pos/POSCAR_00$i/CONTCAR)
    vectorap=$(echo $vectora | sed 's/\s.*$//')
    vectorb=$(sed -n 4"p" <pos/POSCAR_00$i/CONTCAR)
    vectorbp=$(echo $vectorb | sed -e "s/0.0000000000000000 //g" |  sed 's/\s.*$//')
    vectorc=$(sed -n 5"p" <pos/POSCAR_00$i/CONTCAR)
    vectorcp=$(echo $vectorc | sed -e 's/.* //g')
    echo $vectorap $vectorbp $vectorcp
    volume=$(echo $vectorap*$vectorbp*$vectorcp | bc -l)
    echo $volume $energyp >> volen.txt





done 

for i in $(seq 10 99) 
do 

  echo $i
    echo "Final energy from 1st step"
    energy=$(grep "e  e" pos/POSCAR_0$i/OUTCAR | tail -1)
    energyp=$(echo $energy | sed -e "s/.* = //g")
    echo $energyp >> currens.txt
    #echo "Current energy from 2nd step"
    #grep "e  e" pos/POSCAR_00$i2/RELAX/OUTCAR | tail -1
    vectora=$(sed -n 3"p" <pos/POSCAR_0$i/CONTCAR)
    vectorap=$(echo $vectora | sed 's/\s.*$//')
    vectorb=$(sed -n 4"p" <pos/POSCAR_0$i/CONTCAR)
    vectorbp=$(echo $vectorb | sed -e "s/0.0000000000000000 //g" |  sed 's/\s.*$//')
    vectorc=$(sed -n 5"p" <pos/POSCAR_0$i/CONTCAR)
    vectorcp=$(echo $vectorc | sed -e 's/.* //g')
    echo $vectorap $vectorbp $vectorcp
    volume=$(echo $vectorap*$vectorbp*$vectorcp | bc -l)
    echo $volume $energyp >> volen.txt


done
 
for i in $(seq 100 700) 
do 
  echo $i
    echo "Final energy from 1st step"
    energy=$(grep "e  e" pos/POSCAR_$i/OUTCAR | tail -1)
    energyp=$(echo $energy | sed -e "s/.* = //g")
    echo $energyp >> currens.txt
    #echo "Current energy from 2nd step"
    #grep "e  e" pos/POSCAR_00$i2/RELAX/OUTCAR | tail -1
    vectora=$(sed -n 3"p" <pos/POSCAR_$i/CONTCAR)
    vectorap=$(echo $vectora | sed 's/\s.*$//')
    vectorb=$(sed -n 4"p" <pos/POSCAR_$i/CONTCAR)
    vectorbp=$(echo $vectorb | sed -e "s/0.0000000000000000 //g" |  sed 's/\s.*$//')
    vectorc=$(sed -n 5"p" <pos/POSCAR_$i/CONTCAR)
    vectorcp=$(echo $vectorc | sed -e 's/.* //g')
    echo $vectorap $vectorbp $vectorcp
    volume=$(echo $vectorap*$vectorbp*$vectorcp | bc -l)
    echo $volume $energyp >> volen.txt


done 
