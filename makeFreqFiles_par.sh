#!/bin/bash

## ------------------------------------------
# Ported to Bash from csh, 5/2015 - aps
## ------------------------------------------


# ASSIGN MEANINGFUL VARIABLE NAMES
# --------------------------------
# if ( $#argv != 9 ) then
#     echo number of input args: $#argv
#     echo  Wrong number of input arguments!
#     exit 1
# endif

centerLat=$1
rayMaker=$2
readOne=$3
targDir=$4
bot=$5
top=$6
makeRays=$7
runProg=$8
#L_TARG=($9)
prefix=$9
echo $targDir

# ASSIGN VARIABLE VALUES
# ----------------------
numRays=40    # Number of rays to span latitude range
deltLow=0     # deltLow, deltHigh set launch angle (0 = straight up)
deltHigh=0  
dur=5         # Total run time; seconds
latSpread=20  # degrees in latitude to spread numRays across
thisDir=`pwd` # Where we're at yo
divLatNum=50  #
freqStep=1

L_TARG=( 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 \
         2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 \
         3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 \
         4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 )

#L_TARG=( 2.2 )
#set freqs = ( 1000 )

# freqs=( 200 \
#             210 \
#             220 \
#             230 \
#             240 \
#             250 \
#             260 \
#             270 \
#             280 \
#             290 \
#             300 \
#             310 \
#             320 \
#             330 \
#             340 \
#             350 \
#             360 \
#             370 \
#             380 \
#             390 \
#             400 \
#             410 \
#             420 \
#             430 \
#             440 \
#             460 \
#             480 \
#             500 \
#             520 \
#             540 \
#             560 \
#             580 \
#             600 \
#             620 \
#             640 \
#             660 \
#             690 \
#             720 \
#             750 \
#             780 \
#             810 \
#             840 \
#             880 \
#             920 \
#             960 \
#             1000 \
#             1040 \
#             1080 \
#             1120 \
#             1160 \
#             1200 \
#             1250 \
#             1300 \
#             1360 \
#             1420 \
#             1480 \
#             1540 \
#             1600 \
#             1670 \
#             1740 \
#             1800 \
#             1880 \
#             1960 \
#             2050 \
#             2140 \
#             2230 \
#             2320 \
#             2410 \
#             2500 \
#             2600 \
#             2700 \
#             2800 \
#             2900 \
#             3000 \
#             3100 \
#             3200 \
#             3300 \
#             3450 \
#             3600 \
#             3750 \
#             3900 \
#             4050 \
#             4200 \
#             4400 \
#             4600 \
#             4800 \
#             5000 \
#             5200 \
#             5400 \
#             5600 \
#             5800 \
#             6000 \
#             6200 \
#             6400 \
#             6600 \
#             6900 \
#             7200 \
#             7500 \
#             7800 \
#             8100 \
#             8400 \
#             8800 \
#             9200 \
#             9600 \
#             10000 \
#             10500 \
#             11000 \
#             11500 \
#             12000 \
#             12500 \
#             13000 \
#             13500 \
#             14000 \
#             14500 \
#             15000 \
#             16000 \
#             17000 \
#             18000 \
#             19000 \
#             20000 \
#             22000 \
#             24000 \
#             26000 \
#             28000 \
#             30000 \
#             35000 \
#             40000 \
#             45000 \
#             50000 \
#             60000      )

# # Freqs: 32, log-spaced between 200 and 60000
freqs=( 200 \
            240 \
            289 \
            347 \
            418 \
            502 \
            603 \
            725 \
            872 \
            1048 \
            1259 \
            1514 \
            1819 \
            2187 \
            2629 \
            3160 \
            3798 \
            4565 \
            5487 \
            6596 \
            7928 \
            9530 \
            11455 \
            13769 \
            16550 \
            19893 \
            23912 \
            28742 \
            34549 \
            41528 \
            49916 \
            60000 ) 




# Now: 
# create and CD into the correct directory
# copy the correct raymaker file into this directory
# foreach frequency
#   make the correct newray.in file
#   run the raytracer
#   move the newray.dat file to newray[FREQUENCY].dat
#   remove .out, .in., etc. files
#   run damping on this frequency file
# end foreach



# IF DIRECTORY DOESN'T EXIST, CREATE IT
if [ ! -d "$targDir" ]; then
    echo "making output directory"
    mkdir $targDir
fi

if [ ! -d "logs" ]; then
    mkdir logs
fi

cd $targDir



# Copy rayMaker and damping into the right directory
# unalias cp
cp -f $thisDir/runfiles/$rayMaker .
cp -f $thisDir/runfiles/damping .
cp -f $thisDir/runfiles/$readOne .    
cp -f $thisDir/runfiles/newray .
cp -f $thisDir/readOne_par.pbs .
cp -f $thisDir/rayMaker_par.pbs .
cp -f $thisDir/jh.sh .

len=`echo ${freqs[@]} | wc -w`  
lenL=`echo ${L_TARG[@]} | wc -w`

echo "will do " $len "frequences"
echo "will do " $lenL "L-shells"

if [ $top -eq 0 ]; then
    top=$len    
fi


# ---- RayMaker ------
if [ $makeRays -eq 1 ]; then
    i=$bot

    while [ $i -le $top ]; do
        freq=${freqs[$i]}
        echo "i = " $i " freq = " $freq

    # --- Just run serially on single node ---

    # ./$rayMaker $freq $numRays $centerLat $deltLow $deltHigh $dur $latSpread
    # ./newray
    # mv newray.dat newray${freq}.dat
    # rm -f Lprofile.dat latprofile.dat newray.in newray.out 
    # ./damping $freq newray${freq}.dat 1    


    # --- Run in parallel on Nansen ---
        # Create temp directory (each instance of raymaker works with 'newray.in')
        if [ ! -d "$freq" ]; then
            mkdir $freq
        fi

        #cd $freq 
        #cp -f $thisDir/runfiles/newray . # Copy compiled raytracer into working dir
        #cp -f $thisDir/runfiles/damping .
        cp $thisDir/rayMaker_par.pbs ${freq}/rayMaker_par.pbs

        workingDir=`pwd`/${freq}
        codeDir=$thisDir/runfiles
        jobname=${prefix}_ray_${freq}

        # No clue why, but the raytracer segfaults half the time, but only on batchnew.
        # Works fine on batch!
	    qq=batch

#         # nQueued=$(./jh.sh ${prefix} | grep "Queued: " | awk -F ' ' '{print $2}')
#         # echo $nQueued

        qsub -N $jobname -j oe -o log_ray_${freq}.txt -l nodes=1:ppn=1 -l walltime=24:00:00 -q $qq\
            -v codeDir=$codeDir,workingDir=$workingDir,targDir=$targDir,A0=$rayMaker,A1=$freq,A2=$numRays,A3=$centerLat,A4=$deltLow,A5=$deltHigh,A6=$dur,A7=$latSpread \
            rayMaker_par.pbs
        #cd ..    
        ((i++))
    done 

    # Wait for all RayMaker jobs to complete:
    # echo "PWD is " `pwd`
    ./jh.sh ${prefix}    
    while [ $? -eq 1 ]; do
        sleep 10
        ./jh.sh ${prefix}
    done

fi


# Check for any segfaults!!! WHY ARE THERE SEGFAULTS
echo Raymaker Segfaults:
head $targDir/log_ray* | grep "fault" | wc -l


# echo "targDir is: " $targDir
# pwd

if [ $runProg -eq 1 ]; then
    i=$bot
    while [ $i -lt $top ]; do
        ii=$i+1
        for L  in ${L_TARG[@]}; do
            # --- Run in parallel on Nansen ---

            jobname=${prefix}_scatter_${freqs[i]}_$L
            qq=batchnew

            # ----------------- Ugh, be polite. ----------------------
            nQueued=$(./jh.sh ${prefix} | grep "Queued: " | awk -F ' ' '{print $2}')
            while [ $nQueued -gt 50 ]; do
                echo "idling... centerLat = " $centerLat " L = " $L " i= " $i
                sleep 60
                nQueued=$(./jh.sh ${prefix} | grep "Queued: " | awk -F ' ' '{print $2}')
            done



            qsub -N $jobname -j oe -o log_scatter_${freqs[i]}_$L.txt -l nodes=1:ppn=1 -l walltime=24:00:00 -q $qq \
            -v targDir=$targDir,A0=$readOne,A1=${freqs[$i]},A2=${freqs[$ii]},A3=$centerLat,A4=$divLatNum,A5=$freqStep,A6=$L \
            readOne_par.pbs
        
            # --- Just run serially on single node ---
            # ./$readOne ${freqs[$i]} ${freqs[$ii]} \
            #     $centerLat $divLatNum $freqStep $L
        done
        ((i++))
    done
fi

cd $thisDir

exit 0
