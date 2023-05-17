diamond=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/diamond/2.0.9/bin/diamond
RdRp_lib=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Multi-Segments-Virus/data/RdRp/RdRp.all
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
		mkdir -p $shelldir
		#runshell=$shelldir/step09.${sample}.author.blastx2RdRp.sh
        #runshell=$shelldir/step09.${sample}.rmHost.blastx2RdRp.sh
        runshell=$shelldir/step09.${sample}.withHost.blastx2RdRp.sh
        
        echo -e "echo \"start at:\`date\`\"" > $runshell

        # putative_viral_contigs=$sampledir/06.putative_virus/author/${sample}.author.putative_virus_contigs.renamed.fasta
        # blastx2RdRp_dir=$sampledir/09.blastx2RdRp/author
        # blastx2RdRp_out=$blastx2RdRp_dir/${sample}.author-contigs.blastx2RdRp.blastx
        
        # putative_viral_contigs=$sampledir/06.putative_virus/rmHost/${sample}.rmHost.spades.putative_virus_contigs.renamed.fasta
        # blastx2RdRp_dir=$sampledir/09.blastx2RdRp/rmHost
        # blastx2RdRp_out=$blastx2RdRp_dir/${sample}.rmHost-spades-contigs.blastx2RdRp.blastx
        
        putative_viral_contigs=$sampledir/06.putative_virus/withHost/${sample}.withHost.spades.putative_virus_contigs.renamed.fasta
        blastx2RdRp_dir=$sampledir/09.blastx2RdRp/withHost
        blastx2RdRp_out=$blastx2RdRp_dir/${sample}.withHost-spades-contigs.blastx2RdRp.blastx
        mkdir -p $blastx2RdRp_dir
        
        if [ -f $putative_viral_contigs ];then
            echo "$diamond blastx --threads 3 --query $putative_viral_contigs --db $RdRp_lib  --sensitive --max-target-seqs 50 --evalue 1E-5 --block-size 2.0 --index-chunks 1 --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps --out $blastx2RdRp_out" >> $runshell
            echo -e "echo \"end at:\`date\`\"" >> $runshell
            if [ $qsub_status == 'yes' ];then
                cd $shelldir
                rm -f ${runshell}.* core.*
                #qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:1 -l vf=1G,num_proc=1 $runshell
                qsub -clear -cwd -q st_supermem.q -P P20Z10200N0206 -binding linear:2 -l vf=2G,num_proc=2 $runshell
            else
                ls $runshell
                rm -f core.*
            fi
        else
            echo -e "$project\t$sample" >> $assembly_failed_list
	    fi
    done
fi