#put it in dir of LigandMPNN

work_dir=/home/wendao/work/LLPS/metallps/fill_zinc/${1}

for i in $(seq 1 101)
do
python run.py --model_type "ligand_mpnn" --pdb_path ${work_dir}/ZN_${i}/target.pdb --out_folder ${work_dir}/ZN_${i}/ --redesigned_residues "$(tr -d '\n' < ${work_dir}/ZN_${i}/pos.txt)" --checkpoint_ligand_mpnn "./model_params/ligandmpnn_v_32_020_25.pt" --batch_size 10 --number_of_batches 10
done

