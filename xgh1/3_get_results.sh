#!/bin/bash
echo " "
echo "----------------------- Coulomb ingetrals -------------------------"
path0=$(pwd)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
vatom_script="$path0/vatom.py"
if [ ! -f "$vatom_script" ]; then
    vatom_script="$script_dir/vatom.py"
fi

# check defect.input
if [ ! -f "defect.input" ]; then
    echo "There is no 'defect.input'."
    exit 1
fi

q1=$(grep "charged " $path0/defect.input | awk '{print $3}')
sign=${q1:0:1}
q=${q1:1}
sq=$q
path_corr=$path0/image-corr_$sign$sq

E1=$(grep "E_Coul(eV)" $path_corr/REPORT | tail -n 1 | awk '{ printf "%8.5f \n", $2 }')
E0=$(grep "E_Coul(eV)" $path_corr/REPORT.0 |tail -n 1 | awk '{printf "%8.5f \n", $2}')
echo "E_\infty = $E1   E_P = $E0"
echo " "
echo " "
echo "------------------- Potential alignment -------------------"
pos=$(grep "defect" defect.input)
pos=${pos:6}
echo "The defect coordinates: $pos"
path_b=$(grep "bulk" defect.input | awk '{print $2}')
path_0=$(grep "neutral" defect.input | awk '{print $2}')
#echo $path_b $path_0

cd $path_0
# check OUT.VATOM
if [ ! -f "OUT.VATOM" ]; then
    echo "OUT.VATOM is not in $path_0"
    echo "Exit!!!"
    exit 1
fi

echo -e "neutral: \c"
python $vatom_script << EOF
$pos
EOF

cd $path_b

# check OUT.VATOM
if [ ! -f "OUT.VATOM" ]; then    
    echo "OUT.VATOM is not in $path_b"
    echo "Exit!!!"
    exit 1
fi

echo -e "bulk: \c"
python $vatom_script << EOF
$pos
EOF
echo " "
echo "--------------------Finished ( * ~ * )_Y-------------------"
echo " "
