#!/bin/bash

file="occ.input"
file1="occ_spin1.input"
file2="occ_spin2.input"
 #for spin=1  #
if [ -f "$file" ]; then         
    kpoint=$(awk 'NR==2 {print $1}' "$file")
    nd=$(awk 'NR==2 {print $2}' "$file")
    di=$(awk 'NR==2 {print $3}' "$file")
    q=$(awk 'NR==2 {print $4}' "$file")

    for ((j=1; j<=$nd; j++)); do
        for ((i=1; i<=$kpoint; i++)); do
            numdefect=$((di+j-1))
       	    echo -e "1\n1\n$i\nOUT.WG\n$numdefect\n" | convert_wg2rho.x >/dev/null
            mv OUT.WG2RHO IN.RHO$i-$j
	done
    done
    array=($(awk 'NR>2 {for (i=1; i<=NF; i++) print $i}' "$file"))
    index=0
    parameters_array=()
    for ((i=1; i<=$kpoint; i++)); do       
        for ((j=1; j<=$nd; j++)); do
            if [ $index -eq $((kpoint*nd-1)) ]; then
               parameters_array+=("IN.RHO$i-$j ${array[index]}\\nN\\n")
            else
               parameters_array+=("IN.RHO$i-$j ${array[index]}\\nY\\n")
            fi
            index=$((index+1))
        done
    done
   # for parameters in "${parameters_array[@]}"
   # do
    #   echo "$parameters"
   # done
    echo -e "${parameters_array[@]}"|convert_rho_diff.x >/dev/null
    mv OUT.RHO_diff GEN.RHO
 ###for spin =2 ###
else
  if [ -f "$file1" ]; then
    kpoint=$(awk 'NR==2 {print $1}' "$file1")
    nd1=$(awk 'NR==2 {print $2}' "$file1")
    di1=$(awk 'NR==2 {print $3}' "$file1")
    q=$(awk 'NR==2 {print $4}' "$file1")
    echo "kpoint:$kpoint"
    for ((j=1; j<=$nd1; j++)); do
        for ((i=1; i<=$kpoint; i++)); do
            numdefect=$((di1+j-1))
            echo -e "1\n1\n$i\nOUT.WG\n$numdefect\n" | convert_wg2rho.x >/dev/null
            mv OUT.WG2RHO IN.RHO_spin1_$i-$j
        done
    done
    array=($(awk 'NR>2 {for (i=1; i<=NF; i++) print $i}' "$file1"))
    index=0
    parameters_array=()
    for ((i=1; i<=$kpoint; i++)); do
        for ((j=1; j<=$nd1; j++)); do
            if [ $index -eq $((kpoint*nd1-1)) ]; then
               parameters_array+=("IN.RHO_spin1_$i-$j ${array[index]}\\nN\\n")
            else
               parameters_array+=("IN.RHO_spin1_$i-$j ${array[index]}\\nY\\n")
            fi
            index=$((index+1))
        done
    done
    echo -e "${parameters_array[@]}"|convert_rho_diff.x >/dev/null
    mv OUT.RHO_diff GEN_spin1.RHO
   fi
   if [ -f "$file2" ]; then
    kpoint=$(awk 'NR==2 {print $1}' "$file2")
    nd2=$(awk 'NR==2 {print $2}' "$file2")
    di2=$(awk 'NR==2 {print $3}' "$file2")
    q=$(awk 'NR==2 {print $4}' "$file2")
    echo "kpoint:$kpoint"
    for ((j=1; j<=$nd2; j++)); do
        for ((i=1; i<=$kpoint; i++)); do
            numdefect=$((di2+j-1))
            echo -e "1\n1\n$i\nOUT.WG_2\n$numdefect\n" | convert_wg2rho.x >/dev/null
            mv OUT.WG2RHO IN.RHO_spin2_$i-$j
        done
    done
    array=($(awk 'NR>2 {for (i=1; i<=NF; i++) print $i}' "$file2"))
    index=0
    parameters_array=()
    for ((i=1; i<=$kpoint; i++)); do
        for ((j=1; j<=$nd2; j++)); do
            if [ $index -eq $((kpoint*nd2-1)) ]; then
               parameters_array+=("IN.RHO_spin2_$i-$j ${array[index]}\\nN\\n")
            else
               parameters_array+=("IN.RHO_spin2_$i-$j ${array[index]}\\nY\\n")
            fi
            index=$((index+1))
        done
    done
    echo -e "${parameters_array[@]}"|convert_rho_diff.x >/dev/null
    mv OUT.RHO_diff GEN_spin2.RHO
    fi
    if [ ! -f "$file1" ] && [ ! -f "$file2" ]; then
    echo "error: no occ.input exists"
    exit 1
    fi
  if [ ! -f "$file1" ]; then
    mv GEN_spin2.RHO GEN.RHO
  elif [ ! -f "$file2" ]; then
    mv GEN_spin1.RHO GEN.RHO
  else
    echo -e "GEN_spin1.RHO 1\nY\nGEN_spin2.RHO 1\nN\n" | convert_rho_diff.x >/dev/null
  mv OUT.RHO_diff GEN.RHO
  fi
fi

