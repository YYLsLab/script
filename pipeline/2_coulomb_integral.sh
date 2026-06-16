#!/bin/bash
echo "  "
echo "  "
echo "----------------- Step 2 - Coulomb integral----------------"

###########################################################
### 		Step 2 - Coulomb integral		###
###########################################################


############################################################
###		confirm the path			 ###
###		in the folder for q=0			 ###
############################################################

path0=`pwd`

# check defect.input
if [ ! -f "defect.input" ]; then
    echo "There is no "defect.input"."
fi

# check REPORT job.sh  atom.config
for fn in REPORT job.sh atom.config
do
    if [ ! -f "$fn" ];then
	echo "$fn is not here!"
	echo "Exit!!!"
	exit 1
    fi
done

q1=`grep "charged " $path0/defect.input | awk '{print $3}'`
sign=${q1:0:1}
q=${q1:1}
sq=$q

path_corr=$path0/image-corr_$sign$sq
if [ ! -d "$path_corr" ]; then
    mkdir $path_corr      #folder for job=potential
fi
echo "Charge state: $q1     folder:  $path_corr"

############################################################
### generate etot.input.0 etot.input.1 for job=potential ###
############################################################

# still in $path0
PWGK=`sed -n "1p" etot.input`
PREC=`grep "PRECISION" REPORT`
ACCU=`grep "ACCURACY" REPORT`
CONV=`grep "CONVERGENCE" REPORT`
P1=`grep "IN.PSP1" REPORT`
P2=`grep "IN.PSP2" REPORT`
P3=`grep "IN.PSP3" REPORT`
P4=`grep "IN.PSP4" REPORT`
P5=`grep "IN.PSP5" REPORT`
P6=`grep "IN.PSP6" REPORT`
P7=`grep "IN.PSP7" REPORT`
P8=`grep "IN.PSP8" REPORT`
Ecut=`grep "Ecut " REPORT`
Ecut2=`grep "Ecut2 " REPORT`
Ecut2L=`grep "Ecut2L" REPORT | head -1`
N123=`grep " N123 " REPORT`
NS123=`grep "NS123" REPORT`
N123L=`grep "N123L" REPORT`
MP=`grep "MP_N123 " REPORT`
XCFUNC=`grep "XCFUNCTIONAL" REPORT`
cat > etot.input.0 << EOF
$PWGK
 JOB 	 =	POTENTIAL
$PREC
$ACCU
$CONV
 IN.ATOM = 	defect.config
$P1
$P2
$P3
$P4
$P5
$P6
$P7
$P8
$Ecut
$Ecut2
$Ecut2L
$N123
$NS123
$N123L
$MP
$XCFUNC
 NUM_ELECTRON = $q
 IN.RHO	 =	T
 OUT.RHO =	T
 OUT.VR	 =	T
EOF
cat etot.input.0 |tr -s '\n' > etot.input.0.$sq
rm etot.input.0
cp etot.input.0.$sq etot.input.1.$sq
# etot.input.0.$sq for periodic Coulomb integral
# etot.input.1.$sq for double box integral

echo " COULOMB = 	0" >> etot.input.0.$sq
bound=`grep "bound" $path0/defect.input`
bound=" COULOMB = 1"${bound:5}
echo $bound >> etot.input.1.$sq


###########################################################
### 		generate configuration file 		###
###########################################################
echo " 0" > defect.config
sed "1d" atom.config >> defect.config

mv etot.input.0.$sq etot.input.1.$sq defect.config $path_corr
cp OUT.KPT OUT.SYMM job.sh $path_corr
#cp slurm.sh $path_corr # for sbatch


###########################################################
###  		prepare for the potential files		###
########################################################### 
if [ ! -d "$path_corr/0" ]; then
    mkdir $path_corr/0 
fi


# check pseudopotentials
count=1
ps=`grep "IN.PSP1" $path0/REPORT | awk '{print $3}'`
echo " "
echo -e "pseudopotential: \c"
while [ -n "$ps" ]; do
    if [ -f "$ps" ]; then 
	cp $ps $path_corr
	cp $ps $path_corr/0
	echo -e "$ps \c"
    else
	if [ -f "$path0/$ps" ]; then
	    cp $path0/$ps $path_corr
	    cp $path0/$ps $path_corr/0
	    echo -e "$ps \c"
	else
	    echo " "
	    echo "$ps does not exsit."
	    echo "Exit!!!"
	    exit 1
	fi
    fi
    count=$((count+1))
    ps=`grep "IN.PSP"$count $path0/REPORT | awk '{print $3}'`
done
echo " "
echo " "

###########################################################
###		job = potential				###
###	etot.input, IN.RHO, OUT.KPT, OUT.SYMM, job.sh	###
###	potential file					###
###########################################################
echo " "
echo "--------------------- job = potential ---------------------"
echo -e "job=potential, periodic integral, \c"
cd $path_corr

# check IN.RHO OUT.KPT OUT.SYMM job.sh
for fn in IN.RHO OUT.KPT OUT.SYMM job.sh
do
    if [ ! -f "$fn" ]; then
	echo "$fn is not here!"
	echo "Exit!!!"
	exit 1
    fi
done

cd $path_corr/0
ln -sf ../etot.input.0.$sq etot.input
ln -sf ../IN.RHO .
ln -sf ../defect.config .
ln -sf ../OUT.KPT .
ln -sf ../OUT.SYMM .
cp ../job.sh .
#cp ../slurm.sh . # for sbatch
#sbatch slurm.sh  # for sbatch 
if grep -q "#SBATCH" job.sh; then
    echo "Detected Slurm script"
    sbatch job.sh
elif grep -q "#PBS" job.sh; then
    echo "Detected PBS script"
    qsub job.sh
else
    echo "Unrecognized script format"
fi

echo -e "job=potential, double box integral, \c"
cd $path_corr
cp etot.input.1.$sq etot.input
#sbatch slurm.sh # for sbatch
sbatch job.sh  #pbs script
ln -sf 0/OUT.VR_hion OUT.VR_hion.0
ln -sf 0/REPORT REPORT.0
echo " "
echo " "
echo " "
echo "------- 2_Coulomb_integral.sh Finished! ( * _ * )_Y -------"
echo " "
