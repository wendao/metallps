for i in $(seq 1 101)
do
echo $i
mkdir ZN_${i}
cd ZN_${i}
grep ^ATOM ../input.pdb > target.pdb
head -n ${i} ../zn.pdb | tail -n 1 >> target.pdb
python ../../../scripts/near-Zn.py target.pdb > pos.txt
cd ..
done
