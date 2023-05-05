seqkit=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/seqkit/2.0.0/bin/seqkit
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
		runshell=$shelldir/step03.${sample}.rename.contigs.sh
        echo -e "echo \"start at:\`date\`\"" > $runshell
		#raw_contigs=$sampledir/02.Assembly/metaSPAdes/contigs.fasta
		#new_contigs=$sampledir/02.Assembly/metaSPAdes/contigs.renamed.fasta
		#raw_contigs=$sampledir/02.Assembly/rmHost/metaSPAdes/contigs.fasta
		#new_contigs=$sampledir/02.Assembly/rmHost/metaSPAdes/contigs.renamed.fasta
		tmp_contigs=$sampledir/06.putative_virus/${sample}.putative_virus_contigs.fasta
		raw_contigs=$sampledir/06.putative_virus/${sample}.withHost.putative_virus_contigs.fasta
		if [ -f $tmp_contigs -a ! -f $raw_contigs ];then
			mv $tmp_contigs $raw_contigs
		fi
		new_contigs=$sampledir/06.putative_virus/${sample}.withHost.putative_virus_contigs.renamed.fasta
		if [ -f $raw_contigs ];then
            ## 更改序列的名称
            echo "$seqkit replace -p \"(\w+\.\d+)\" -r \"${sample}~\\\${1}\" $raw_contigs > $new_contigs" >> $runshell
            echo -e "echo \"end at:\`date\`\"" >> $runshell

            if [ $qsub_status == 'yes' ];then
                cd $shelldir
                rm -f ${runshell}.* core.*
                qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:1 -l vf=2G,num_proc=1 $runshell
            else
                ls $runshell
                rm -f core.*
            fi
        else
            echo -e "$project\t$sample" >> $assembly_failed_list
		fi
    done
fi

