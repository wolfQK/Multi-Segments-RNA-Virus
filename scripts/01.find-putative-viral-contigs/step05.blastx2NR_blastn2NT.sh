diamond=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/diamond/2.0.9/bin/diamond
blastn=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/blast/2.10.1/bin/blastn
NR_db=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/NCBI-BLAST/2021-05/NR/Diamond/NR.with.taxonomy
BLASTDB=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/NCBI-BLAST/2021-05/NT/BLASTDB/download
NT_db=$BLASTDB/nt
assembly_failed_list=$PWD/Assembly.failed.sample_list
echo -e "project\tsample" > $assembly_failed_list

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
project=`basename $poutdir`
if [ -z $sample_list ] || [ -z $poutdir ] || [ -z $qsub_status ];then
    echo "Some or all of the parameters are empty";
    helpFunction
else
    cat $sample_list|while read sample
    do
        sampledir=$poutdir/samples/$sample
        shelldir=$sampledir/shell
        runshell=$shelldir/step04.${sample}.blastx2NR-blastn2NT.sh
       
        contigs_fa=$sampledir/02.Assembly/metaSPAdes/contigs.renamed.fasta

        if [ -s $contigs_fa ];then
            echo -e "echo \"start at:\`date '+%F %T'\`\"" > $runshell
            blastx2NR_dir=$sampledir/04.blastx2NR
            blastx2NR_out=$blastx2NR_dir/${sample}.blastx2NR.blastx
            mkdir -p $blastx2NR_dir
            echo " $diamond blastx --threads 4 --query $contigs_fa --db $NR_db --sensitive --max-target-seqs 80 --evalue 1E-5 --block-size 2.0 --index-chunks 1 --out $blastx2NR_out --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps staxids sscinames sskingdoms skingdoms sphylums stitle qcovhsp scovhsp" >> $runshell
            blastn2NT_dir=$sampledir/05.blastn2NT
            blastn2NT_out=$blastn2NT_dir/${sample}.blastn2NT.blastn
            mkdir -p $blastn2NT_dir
            echo "export BLASTDB=$BLASTDB" >> $runshell
            echo " $blastn -num_threads 4 -query $contigs_fa -db $NT_db -out $blastn2NT_out -evalue 1E-5 -max_target_seqs 80 -outfmt \"6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle qcovs qcovhsp qcovus\"" >> $runshell
            echo -e "echo \"end at:\`date '+%F %T'\`\"" >> $runshell
            
            if [ $qsub_status == 'yes' ];then
                cd $shelldir
                rm -f ${runshell}.* core.*
                qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:4 -l vf=80G,num_proc=4 $runshell
            else
                ls $runshell
                rm -f ${runshell}.* core.*
            fi
        else
            echo -e "$project\t$sample" >> $assembly_failed_list
        fi
    done
fi

