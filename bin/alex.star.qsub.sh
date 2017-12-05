set -x

dirFASTQ="./"
filesIn=`cd $dirFASTQ; ls ERR188*_1.fastq.gz`
filesN=`echo $filesIn | wc -w`
parValues="0.55 0.66 0.75 0.85 0.9 0.92 0.94 0.96 0.98 0.99"
parValuesV=( $parValues )
valuesN=${#parValuesV[*]}

echo $filesN $valuesN $((filesN*valuesN))

optEnc="--outFilterType BySJout   --outFilterMultimapNmax 20   --outFilterMismatchNmax 999   --outFilterMismatchNoverReadLmax 0.04   \
--alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJoverhangMin 8   --alignSJDBoverhangMin 1   --sjdbScore 1"

qsub -j y -v PATH -pe threads 10 -l tmp_free=200G,m_mem_free=4G -cwd -t 1-$((filesN*valuesN)) \
  submSTAR.sh \
 $dirFASTQ \"$filesIn\" `pwd` \
 \" $optEnc --outSAMtype None --quantMode GeneCounts --outFilterMatchNminOverLread 0 \" \
 --outFilterScoreMinOverLread \"$parValues\"   8    _1    _2


