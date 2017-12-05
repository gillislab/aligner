#$ -e ./Logs/err.$JOB_ID.$TASK_ID
#$ -o ./Logs/out.$JOB_ID.$TASK_ID

set -vx

##
# make sure you use your local version of STAR
STAR=STAR
##
genome=/home/genome_dir/

inDir=$1
inRead1=( $2 )

# directory to check for completed runs
checkDir=$3

# all fixed star parameters
STARoptions=$4

# name of the series parameter
parName=$5

# values of the series parameter
parValues=( $6 )

# number of threads
nThreads=$7

# read1 name
read1name=$8

# read2 name
read2name=$9

nFiles=${#inRead1[*]}
nValues=${#parValues[*]}

echo Number of files $nFiles
echo Number of the parameter values nValues

outPrefix=`pwd`


inI=$((SGE_TASK_ID-1))

fileI=$((inI/nValues))
valueI=$((inI-fileI*nValues))

read1=${inRead1[$fileI]}
read2=${read1/$read1name/$read2name}

value=${parValues[valueI]}

echo value  =$valueI $value
echo fileI  =$fileI
echo Read1: $read1
echo Read2: $read2
ls -l $inDir/$read1
ls -l $inDir/$read2

readSuffix=${read1##*\.}

outDir=$read1
outDir=${outDir/$read1name/""}
outDir=${outDir/.*/""}
echo $read1 $readSuffix $outDir

outDir=$outDir/$value

mkdir -p $outDir
cd $outDir

ls -lh $checkDir/$outDir/Log.final.out

if [ -s "$checkDir/$outDir/Log.final.out" ]
then
   echo "Run was completed"
   exit 0
fi

if [ "$readSuffix" == "gz" ]
then
   readComm=zcat
   read1=$inDir$read1
   read2=$inDir$read2
fi
if [ "$readSuffix" == "tgz" ]
then
   echo `date` : decompressing tgz files

   mkdir $TMPDIR/R1
   mkdir $TMPDIR/R2
   tar xfz $inDir$read1 --directory=$TMPDIR/R1
   tar xfz $inDir$read2 --directory=$TMPDIR/R2
   read1=`find $TMPDIR/R1 -name "*.*" | sort | tr '\n' ','`
   read2=`find $TMPDIR/R2 -name "*.*" | sort | tr '\n' ','`
   readComm="cat"

   echo `date` : done
fi


$STAR --genomeDir $genome --readFilesIn $read1 $read2 --readFilesCommand $readComm --runThreadN $nThreads $STARoptions $parName $value



