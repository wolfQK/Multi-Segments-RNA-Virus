interproscan=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/interproscan/5.56_89.0/interproscan.sh

helpFunction()
{
	echo ""
    echo "Usage: $0 -l sample_list -o outdir of sample -r qsub or not"
    echo -e "\t-l which 'sample_list' is used to run find_putative_virus_contigs-proteins"
    echo -e "\t-o which parent directory to write the results"
    echo -e "\t-r qsub or not (yes|no)"
    exit 1 # Exit script after printing help
}

while getopts "l:o:r:" opt
do
	case "$opt" in
		l) sample_list="$OPTARG";;
		o) poutdir="$OPTARG";;
		r) qsub_status=$"$OPTARG";;
		?) helpFunction;;
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
		runshell=$shelldir/step15.${sample}.protein.annotation.interproscan.sh
		echo -e "echo \"start at:\`date '+%F %T'\`\"" > $runshell
		putative_contigs_dir=$sampledir/06.putative_virus
		annotation_dir=$putative_contigs_dir/annotation-interproscan
		mkdir -p $annotation_dir
		#作者的contigs
		author_protein_faa=$putative_contigs_dir/${sample}.author-contigs.putative_virus_contigs.gene.faa
		tmp_author_protein_faa=$putative_contigs_dir/${sample}.author-contigs.putative_virus_contigs.gene.faa.tmp
		echo "sed 's/*//g' $author_protein_faa > $tmp_author_protein_faa" >> $runshell
		author_protein_interproscan_out=$annotation_dir/${sample}.author-contigs.putative_virus_contigs.gene
		echo "$interproscan -appl Pfam -i $tmp_author_protein_faa -b $author_protein_interproscan_out -t p -dp && rm -fr $annotation_dir/temp && rm -f $tmp_author_protein_faa" >> $runshell
		echo -e "echo \"end at:\`date '+%F %T'\`\"" >> $runshell
		if [ $qsub_status == 'yes' ];then
            cd $shelldir
            rm -f ${runshell}.*
            qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:4 -l vf=40G,num_proc=4 $runshell
        else
            ls $runshell
        fi
    done
fi
