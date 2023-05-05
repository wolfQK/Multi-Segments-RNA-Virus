fasterq_dump=/ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/luoyunzhe/software/sratoolkit/sratoolkit.2.9.2-centos_linux64/bin/fasterq-dump
poutdir=/jdfssz1/ST_HEALTH/P20Z10200N0206/fengqikai/wuhan_Virology/segmented_virome/biye/Pub_data
sra_list=$1

cat $sra_list|while read sra_file
do
    sampledir=`dirname $sra_file`
    sra=`basename $sampledir`
    prefix=`basename $sra_file`
    shelldir=$sampledir/shell
    runshell=$shelldir/step00.${sra}.sra2fq.sh
    datadir=$sampledir/00.Raw_reads
    mkdir -p $shelldir $datadir
    echo -e "echo \"start at:\`date\`\"" > $runshell
    echo "$fasterq_dump $sra_file --split-3 --force --threads 3 --outdir $datadir --progress" >> $runshell
    raw_fq1=$datadir/${prefix}_1.fastq
    new_fq1=$datadir/${sra}_1.fastq.gz
    raw_fq2=$datadir/${prefix}_2.fastq
    new_fq2=$datadir/${sra}_2.fastq.gz
    echo "gzip -c $raw_fq1 > $new_fq1 && rm -f $raw_fq1" >> $runshell
    echo "gzip -c $raw_fq2 > $new_fq2 && rm -f $raw_fq2" >> $runshell
    echo -e "echo \"end at:\`date\`\"" >> $runshell
    cd $shelldir
    qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:3 -l vf=6G,num_proc=3 $runshell
done

