#!/bin/bash

############### Attention ##################
# make sure you know the path of the file: #
#   defect.input			   #
#   job.pbs (the pbs script)	   	   #
#   pseudopotential			   #
# and the scriptes:			   #
#   generate_occ.py			   #
#   generate_defect_rho.x		   #
#   vatom.x				   #
############################################

################# steps ###################
#   1 get defect rho                      #
#   2 Coulomb integral			  #
#   3 PA 	                          #
###########################################

#1   Step 1
######## prepare defect information #######
# the format is as follows:               #
# dielectric constant 5.20                #
# defect 0.50000 0.50000 0.50000  	  #
# bound 0.00000 0.00000 0.00000           #
# charge state  1 (or                     #
# charge state -1 or                      #
# charge state  2 or                      #
# charge state -2 remain the same format!)#
# bulk /.....				  #
# neutral /.......			  #
###########################################

echo "                                 "
echo "                                 "
echo "-------- Step 1 Generate the defect charge density --------"
echo "Check the paths for generate_occ.py, generate_defect_rho.x, vatom.x"
echo "Check the files: defect.input, OUT.OCC0, OUT.OCC1, job.pbs."

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# check defect.input
if [ ! -f "defect.input" ]; then
    echo "defect.input is not here."
    echo "Exit!!!"
    exit 1
fi

# check OUT.OCC0 OUT.OCC1
if [ ! -f "OUT.OCC0" ] || [ ! -f "OUT.OCC1" ]; then
	echo "OUT.OCC0 or OUT.OCC1 is not here."
	echo "Exit!!!"
        exit 1
fi

# check job.pbs OUT.WG REPORT atom.config OUT.GKK OUT.RHO OUT.KPT OUT.SYMM
for fn in job.sh OUT.WG REPORT  OUT.GKK OUT.RHO OUT.KPT OUT.SYMM
do
    if [ ! -f "$fn" ]; then
	echo "$fn is not here!!!"
	echo "Exit!!!"
	exit 1
    fi
done
path0=`pwd`
echo "The current directory:  "$path0

# check generate_occ.py
path_generate_occ="${path0}/generate_occ.py"	        # check the path
if [ -f "$path_generate_occ" ]; then
        echo "generate_occ.py:  $path_generate_occ"
elif [ -f "$script_dir/generate_occ.py" ]; then
        path_generate_occ="$script_dir/generate_occ.py"
        echo "generate_occ.py:  $path_generate_occ"
else
	which generate_occ.py > /dev/null 2>&1
	if (($?));then
                echo "Can not find generate_occ.py"
                exit 1
        else
		path_generate_occ=`which generate_occ.py`
		echo "generate_occ.py:  $path_generate_occ"
#    if [ -f "$path_generate_occ" ]; then
#        echo "generate_occ.py:        $path_generate_occ"
#    else
#	echo "Please check the path of generate_occ.py."
#	echo "Exit!!!"
#	exit 1
    	fi
fi
#echo "generate_occ.py:        $path_generate_occ"

# check generate_defect_rho.sh
#path_generate_rho="${path0}/generate_defect_rho.sh"
#	echo "Please check the path of generate_defect_rho.x"
#	echo "Exit!!!"
#	exit 1

###################################################
###	check the path				###
###################################################
q1=`grep "charged " $path0/defect.input | awk '{print $3}'`

sign=${q1:0:1}
q=${q1:1}
sq=$q
path_corr=$path0/image-corr_$sign$sq  # path for image charge correction
if [ ! -r "$path_corr" ]; then
    mkdir $path0/image-corr_$sign$sq
fi
echo "charge state: $sign$sq        folder: image-corr_$sign$sq"


###################################################
###	generate IN.RHO for job=potential	###
###################################################

cp OUT.OCC0 OUT.OCC1 OUT.KPT OUT.SYMM $path_corr

cd $path_corr
echo "                                             "
echo "                                             "
echo "-------------------- generate_occ.py ----------------------"
python $path_generate_occ   # check the path
cp $path_corr/occ* $path0
cd $path0
echo "                                             "
echo "                                             "
echo "----------------- generate_defect_rho.sh -------------------"
if [ -f "$path0/generate_defect_rho.sh" ]; then
    path_generate_rho="$path0/generate_defect_rho.sh"
elif [ -f "$script_dir/generate_defect_rho.sh" ]; then
    path_generate_rho="$script_dir/generate_defect_rho.sh"
else
    path_generate_rho=""
fi

if [ -n "$path_generate_rho" ]; then
    echo "sh $path_generate_rho"
    sh "$path_generate_rho"
   # source "$path_generate_rho"
else
    echo "Cannot find generate_defect_rho.sh"
    exit 1
fi
mv $path0/GEN.RHO $path_corr/IN.RHO
rm IN.RHO*-*
echo " "
echo "------------ 1_get_rho.sh Finished ( * _ * )_Y ------------"
echo " "
