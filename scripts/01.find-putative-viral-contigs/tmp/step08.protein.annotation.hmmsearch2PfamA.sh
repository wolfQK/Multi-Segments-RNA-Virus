hmmscan=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/hmmer/3.3.2/bin/hmmscan
pfam_db=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/Pfam_A/Pfam-A.hmm

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
		runshell=$shelldir/step14.${sample}.protein.annotation.hmmsarch2PfamA.sh
		echo -e "echo \"start at:\`date '+%F %T'\`\"" > $runshell
		putative_contigs_dir=$sampledir/06.putative_virus
		annotation_dir=$putative_contigs_dir/annotation-hmmscan
		mkdir -p $annotation_dir
		#作者的contigs
		author_protein_faa=$putative_contigs_dir/${sample}.author-contigs.putative_virus_contigs.gene.faa
		author_protein_hmmscan_out=$annotation_dir/${sample}.author-contigs.putative_virus_contigs.gene.hmmscan2PfamA
		author_protein_hmmscan_tblout=$annotation_dir/${sample}.author-contigs.putative_virus_contigs.gene.tblout
		author_protein_hmmscan_domtblout=$annotation_dir/${sample}.author-contigs.putative_virus_contigs.gene.domtblout
		author_protein_hmmscan_pfamtblout=$annotation_dir/${sample}.author-contigs.putative_virus_contigs.gene.pfamtblout
		echo "$hmmscan -o $author_protein_hmmscan_out --tblout $author_protein_hmmscan_tblout --domtblout $author_protein_hmmscan_domtblout --pfamtblout $author_protein_hmmscan_pfamtblout --acc -E 1E-5 --domE 1E-5 --incE 1E-5 --incdomE 1E-5 --cpu 4 $pfam_db $author_protein_faa" >> $runshell
		#组装的contigs
		rmHost_protein_faa=$putative_contigs_dir/${sample}.rmHost.spades.contigs.putative_virus_contigs.gene.faa		
		rmHost_protein_hmmscan_out=$putative_contigs_dir/${sample}.rmHost.spades.contigs.putative_virus_contigs.gene.hmmscan2PfamA
		rmHost_protein_hmmscan_out=$annotation_dir/${sample}.rmHost.spades.contigs.putative_virus_contigs.gene.hmmscan2PfamA
        rmHost_protein_hmmscan_tblout=$annotation_dir/${sample}.rmHost.spades.contigs.putative_virus_contigs.gene.tblout
        rmHost_protein_hmmscan_domtblout=$annotation_dir/${sample}.rmHost.spades.contigs.putative_virus_contigs.gene.domtblout
        rmHost_protein_hmmscan_pfamtblout=$annotation_dir/${sample}.rmHost.spades.contigs.putative_virus_contigs.gene.pfamtblout
        echo "$hmmscan -o $rmHost_protein_hmmscan_out --tblout $rmHost_protein_hmmscan_tblout --domtblout $rmHost_protein_hmmscan_domtblout --pfamtblout $rmHost_protein_hmmscan_pfamtblout --acc -E 1E-5 --domE 1E-5 --incE 1E-5 --incdomE 1E-5 --cpu 4 $pfam_db $rmHost_protein_faa" >> $runshell

		echo -e "echo \"end at:\`date '+%F %T'\`\"" >> $runshell
		if [ $qsub_status == 'yes' ];then
            cd $shelldir
            rm -f ${runshell}.*
            qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:4 -l vf=8G,num_proc=4 $runshell
        else
            ls $runshell
        fi
    done
fi
