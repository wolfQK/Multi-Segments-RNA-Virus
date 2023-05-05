fastp=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/fastp/0.23.1/fastp
prinseq=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/prinseq-plus-plus/1.2.4/bin/prinseq++
sortmerna=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/sortmerna/4.3.4/bin/sortmerna
metaspades=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/spades/3.15.4/bin/metaspades.py
seqkit=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/seqkit/2.0.0/bin/seqkit
samtools=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/samtools/1.10/bin/samtools

rRNAdir=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Multi-segmentsVirusFinder/data/rRNA
indexdir=$rRNAdir/SortMeRNA_index
rfam_5p8s=$rRNAdir/fasta/rfam-5.8s-database-id98.dna.fasta
rfam_5s=$rRNAdir/fasta/rfam-5s-database-id98.dna.fasta
silva_arc_16s=$rRNAdir/fasta/silva-arc-16s-id95.fasta
silva_arc_23s=$rRNAdir/fasta/silva-arc-23s-id98.fasta
silva_bac_16s=$rRNAdir/fasta/silva-bac-16s-id90.fasta
silva_bac_23s=$rRNAdir/fasta/silva-bac-23s-id98.fasta
silva_euk_18s=$rRNAdir/fasta/silva-euk-18s-id95.fasta
silva_euk_28s=$rRNAdir/fasta/silva-euk-28s-id98.fasta

helpFunction()
{
    echo ""
    echo "Usage: $0 -l sample_list -o outdir of sample -r qsub or not"
    echo -e "\t-l which 'sample_list' is used to run QC2Assembly"
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
        sampledir=$poutdir/samples/$sample
        shelldir=$sampledir/shell
        stat_dir=$sampledir/03.stats
        runshell=$shelldir/step02.${sample}.QC2Assembly.sh
        mkdir -p $shelldir $stat_dir
        echo -e "echo \"start at:\`date\`\"" > $runshell
        ##raw data
        rawdata_dir=$sampledir/00.Raw_reads
		tmp_fq1=$sampledir/${sample}_1.fastq.gz
		tmp_fq2=$sampledir/${sample}_2.fastq.gz
        raw_fq1=$rawdata_dir/${sample}_1.fastq.gz
        raw_fq2=$rawdata_dir/${sample}_2.fastq.gz
		mkdir -p $rawdata_dir
		if [ -f $tmp_fq1 -a -f $tmp_fq2 ];then
			mv $tmp_fq1 $raw_fq1
			mv $tmp_fq2 $raw_fq2
		fi
        raw_fq1_stats=$stat_dir/${sample}.raw.1.stats
        raw_fq2_stats=$stat_dir/${sample}.raw.2.stats
        echo "#$seqkit stats $raw_fq1 > $raw_fq1_stats" >> $runshell
        echo "#$seqkit stats $raw_fq2 > $raw_fq2_stats" >> $runshell
        ##去除低质量和重复
        fastp_dir=$sampledir/01.QC/01.fastp
        mkdir -p $fastp_dir
        fastp_fq1=$fastp_dir/${sample}_1.fastp.fastq.gz
        fastp_fq2=$fastp_dir/${sample}_2.fastp.fastq.gz
        json_out=$fastp_dir/${sample}.fastp.json
        html_out=$fastp_dir/${sample}.fastp.html
        echo "#$fastp --thread 6 -i $raw_fq1 -o $fastp_fq1 -I $raw_fq2 -O $fastp_fq2 --detect_adapter_for_pe --dedup --dup_calc_accuracy 3 --qualified_quality_phred 20 --n_base_limit 5 --average_qual 20 --length_required 50 --low_complexity_filter --correction --json $json_out --html $html_out --report_title \"$sample\"" >> $runshell
        fastp_fq1_stat=$stat_dir/${sample}.fastp.1.stats
        fastp_fq2_stat=$stat_dir/${sample}.fastp.2.stats
        echo "#$seqkit stats $fastp_fq1 > $fastp_fq1_stat" >> $runshell
        echo "#$seqkit stats $fastp_fq2 > $fastp_fq2_stat" >> $runshell
        ##去除低复杂度
        prinseq_dir=$sampledir/01.QC/02.prinseq
        mkdir -p $prinseq_dir
        prinseq_fq1=$prinseq_dir/${sample}_1.prinseq.fastq.gz
        prinseq_fq2=$prinseq_dir/${sample}_2.prinseq.fastq.gz
        echo "#$prinseq -threads 6 -fastq $fastp_fq1 -fastq2 $fastp_fq2 -lc_entropy=0.5 -lc_dust=0.5 -out_gz -out_good $prinseq_fq1 -out_good2 $prinseq_fq2 -out_single /dev/null -out_single2 /dev/null -out_bad /dev/null -out_bad2 /dev/null" >> $runshell
        prinseq_fq1_stat=$stat_dir/${sample}.prinseq.1.stats
        prinseq_fq2_stat=$stat_dir/${sample}.prinseq.2.stats
        echo "#$seqkit stats $prinseq_fq1 > $prinseq_fq1_stat" >> $runshell
        echo "#$seqkit stats $prinseq_fq2 > $prinseq_fq2_stat" >> $runshell
        echo "#if [ -s $prinseq_fq1_stat -a -s $prinseq_fq2_stat ];then" >> $runshell
        echo "# rm -f $fastp_fq1 $fastp_fq2" >> $runshell
        echo "#fi" >> $runshell
        ##去除rRNA
        sortmerna_dir=$sampledir/01.QC/03.sortmerna
        mkdir -p $sortmerna_dir
        aligned_prefix=$sortmerna_dir/out/${sample}_rRNA
        other_prefix=$sortmerna_dir/out/${sample}_rmrRNA
        kvdb_dir=$sortmerna_dir/kvdb
        readb=$sortmerna_dir/readb
        echo "#rm -fr $sortmerna_dir/*" >> $runshell
        echo "#$sortmerna --threads 6 --reads $prinseq_fq1 --reads $prinseq_fq2 --idx-dir $indexdir --ref $rfam_5p8s --ref $rfam_5s --ref $silva_arc_16s --ref $silva_arc_23s --ref $silva_bac_16s --ref $silva_bac_23s --ref $silva_euk_18s --ref $silva_euk_28s --no-best --num_alignments 1 --workdir $sortmerna_dir --fastx --aligned $aligned_prefix --other $other_prefix --paired_out --out2 " >> $runshell
        echo "#rm -fr $kvdb_dir" >> $runshell
        echo "#rm -fr $readb" >> $runshell
        clean_fq1=${other_prefix}_fwd.fq.gz
        clean_fq2=${other_prefix}_rev.fq.gz
        clean_fq1_stats=$stat_dir/${sample}.sortmerna.1.stats
        clean_fq2_stats=$stat_dir/${sample}.sortmerna.2.stats
        echo "#$seqkit stats $clean_fq1 > $clean_fq1_stats" >> $runshell
        echo "#$seqkit stats $clean_fq2 > $clean_fq2_stats" >> $runshell
        echo "#if [ -s $clean_fq1_stats -a -s $clean_fq2_stats ];then" >> $runshell
        echo "# rm -f $prinseq_fq1 $prinseq_fq2">> $runshell
        echo "#fi" >> $runshell
        ##组装
        assembly_dir=$sampledir/02.Assembly/metaSPAdes
        mkdir -p $assembly_dir
        #echo "$metaspades -1 $clean_fq1 -2 $clean_fq2 -t 6 -m 20 -o $assembly_dir" >> $runshell
		echo "rm -fr $assembly_dir/*" >> $runshell
        echo "$metaspades -1 $clean_fq1 -2 $clean_fq2 -t 6 -m 80 -o $assembly_dir" >> $runshell

        echo -e "echo \"end at:\`date\`\"" >> $runshell
        if [ $qsub_status == 'yes' ];then
            cd $shelldir
			rm -f ${runshell}.* core.*
            qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:6 -l vf=80G,num_proc=6 $runshell
			#qsub -clear -cwd -q st_supermem.q -P P20Z10200N0206 -binding linear:6 -l vf=100G,num_proc=6 $runshell
        else
			cd $shelldir
            ls $runshell
			rm -f ${runshell}.* core.*
        fi
    done
fi

