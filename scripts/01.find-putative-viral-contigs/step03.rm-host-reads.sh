seqkit=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/seqkit/2.0.0/bin/seqkit
samtools=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/samtools/1.10/bin/samtools
bowtie2=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/bowtie2/2.4.2/bin/bowtie2
index=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/Mosquitos/genome/Genus/AAC/Aedes-Anopheles-Culex.bt2
star=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/star/2.7.10a/bin/STAR
star_index=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/Mosquitos/genome/Genus/AAC/STAR

# remove host reads using bowtie2 and AAC(Aedes Anopheles Culex) genomes
helpFunction()
{
    echo ""
    echo "Usage: $0 -l sample_list -o outdir of sample -r qsub or not"
    echo -e "\t-l which 'sample_list' is used to run rm-host-reads"
    echo -e "\t-o which parent directory to write the results"
    echo -e "\t-r qsub or not (yes|no)"
    exit 1 # Exit script after printing help
}

while getopts "l:o:r:" opt
do
   case "$opt" in
        l ) sample_list="$OPTARG" ;;
        o ) poutdir="$OPTARG" ;;
        r ) qsub_status="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done
if [ -z $sample_list ] || [ -z $poutdir ] || [ -z $qsub_status ];then
    echo "Some or all of the parameters are empty";
    helpFunction
else
    cat $sample_list|while read sample
    do
        sampledir=$poutdir/$sample
        shelldir=$sampledir/shell
        stat_dir=$sampledir/03.stats
        #runshell=$shelldir/step06.${sample}.rm-host-reads.sh
        runshell=$shelldir/step06.${sample}.STAR.rm-host-reads.sh
        mkdir -p $shelldir $stat_dir
        echo -e "echo \"start at:\`date\`\"" > $runshell
        ## raw clean reads
        sortmerna_dir=$sampledir/01.QC/03.sortmerna
        other_prefix=$sortmerna_dir/out/${sample}_rmrRNA
        raw_clean_fq1=${other_prefix}_fwd.fq.gz
        raw_clean_fq2=${other_prefix}_rev.fq.gz
		new_clean_fq1=${other_prefix}_rm-host.1.fastq.gz
		new_clean_fq2=${other_prefix}_rm-host.2.fastq.gz
        ## bowtie2 remove host reads
		echo "# $bowtie2 -p 2 -x $index -1 $raw_clean_fq1 -2 $raw_clean_fq2|$samtools fastq -f 13 -1 $new_clean_fq1 -2 $new_clean_fq2 --threads 2 -" >> $runshell
		new_clean_fq1_stats=$stat_dir/${sample}.rm-host-reads.1.stats
		new_clean_fq2_stats=$stat_dir/${sample}.rm-host-reads.2.stats
		echo "# $seqkit stats $new_clean_fq1 > $new_clean_fq1_stats" >> $runshell
		echo "# $seqkit stats $new_clean_fq2 > $new_clean_fq2_stats" >> $runshell
		##STAR remove host reads
		star_prefix=${other_prefix}_rm-host.STAR
		tmp_clean_fq1=${star_prefix}Unmapped.out.mate1
		tmp_clean_fq2=${star_prefix}Unmapped.out.mate2
		final_clean_fq1=${other_prefix}_rm-host.STAR.1.fastq.gz
		final_clean_fq2=${other_prefix}_rm-host.STAR.2.fastq.gz
		echo "# $star --runThreadN 3 --runMode alignReads --readFilesCommand zcat --genomeDir $star_index --readFilesIn $new_clean_fq1 $new_clean_fq2 --outFileNamePrefix $star_prefix --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx " >> $runshell
		echo "gzip $tmp_clean_fq1" >> $runshell
		echo "gzip $tmp_clean_fq2" >> $runshell
		echo "mv ${tmp_clean_fq1}.gz $final_clean_fq1" >> $runshell
		echo "mv ${tmp_clean_fq2}.gz $final_clean_fq2" >> $runshell
		final_clean_fq1_stats=$stat_dir/${sample}.rm-host.STAR.1.stats
		final_clean_fq2_stats=$stat_dir/${sample}.rm-host.STAR.2.stats
		echo "$seqkit stats $final_clean_fq1 > $final_clean_fq1_stats" >> $runshell
		echo "$seqkit stats $final_clean_fq2 > $final_clean_fq2_stats" >> $runshell
        echo -e "echo \"end at:\`date\`\"" >> $runshell

        if [ $qsub_status == 'yes' ];then
            cd $shelldir
			rm -f ${runshell}.* core.*
            qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:2 -l vf=2G,num_proc=2 $runshell
            #qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:2 -l vf=30G,num_proc=2 $runshell
			#qsub -clear -cwd -q st_supermem.q -P P20Z10200N0206 -binding linear:3 -l vf=310G,num_proc=3 $runshell
        else
            ls $runshell
			rm -f ${runshell}.* core.*
        fi
    done
fi

