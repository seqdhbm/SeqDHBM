#!/bin/bash
clear
sequence_home="~/Public/mauricio/IL-1_homology/"
fasta_file="IL-1_formatted.fasta"
yasara_exe="~/Downloads/yasara/yasara -txt"
macro_Homology="/home/ajay/Data/IL36/Docking/PC1/mcr/dock_run.mcr"
#Declare an array variable
#declare -a arr_1=("WT")
declare -a arr_1=("WT" "Y108S" "H109A" "Y108SH109A")
declare -a arr_2=("C136SP137A" "H109AC136SP137A" "Y108SC136SP137A" "Y108SH109AC136SP137A")

dock_no=($(seq 1 1 10))

cd $sequence_home

for i in "${arr_1[@]}"
do
	cd $sequence_home/"$i"_Docking
	for count in "${dock_no[@]}"
	do
		#echo $docking_home/"$i"_Docking/"$i"_"$count"_BlindDock
		cd $sequence_home/"$i"_Docking/"$i"_"$count"_BlindDock
		macro_target="MacroTarget='"${sequence_home}"/"${i}"_Docking/"${i}"_"${count}"_BlindDock/"${count}"'"
		echo ${macro_target}
		sed -i "s|^MacroTarget=.*|${macro_target}|" "${macro_BlindDock}"
		cd ~/YASARA/yasara
		./yasara -txt "${macro_Homology}" 
	done

	cd ..
done