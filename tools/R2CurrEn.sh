for i in $(seq 1 9) 
do 
    echo $i
    echo "Final energy from 1st step"
    grep "e  e" pos/POSCAR_00$i/RELAX/RELAX2/OUTCAR | tail -1
    #echo "Current energy from 2nd step"
    #grep "e  e" pos/POSCAR_00$i2/RELAX/OUTCAR | tail -1
done 

for i in $(seq 10 99) 
do 
    echo $i
    grep "e  e" pos/POSCAR_0$i/RELAX/RELAX2/OUTCAR | tail -1
done
#echo "hello" 
#for i in $(seq 100 300) 
#do 
#    echo $i
#    grep "e  e" pos/POSCAR_$i/RELAX/OUTCAR | tail -1
#done 
