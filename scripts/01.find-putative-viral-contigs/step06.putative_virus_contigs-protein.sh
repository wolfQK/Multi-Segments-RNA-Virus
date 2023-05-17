#poutdir=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California
scriptsdir=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Multi-segmentsVirusFinder/scripts
python=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/gatk/3.7/bin/python
besthit_filter_py=$scriptsdir/python/besthit_filter.py
seqkit=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/seqkit/2.0.0/bin/seqkit
prodigal=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/prodigal/2.6.3/bin/prodigal

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
		runshell=$shelldir/step05.${sample}.putative_virus_contigs-protein.sh
		##开始
		echo -e "echo \"start at:\`date '+%F %T'\`\"" > $runshell	
		
		contigs_fna=$sampledir/02.Assembly/metaSPAdes/contigs.renamed.fasta
		contigs_id=$sampledir/02.Assembly/metaSPAdes/contigs.renamed.id
		blastx2NR_dir=$sampledir/04.blastx2NR
		blastx2NR_out=$blastx2NR_dir/${sample}.blastx2NR.blastx
		blastx2NR_besthit=$blastx2NR_dir/${sample}.blastx2NR.besthit
		blastx2NR_filter_id=$blastx2NR_dir/${sample}.blastx2NR.filtered.id
		blastn2NT_dir=$sampledir/05.blastn2NT
    	blastn2NT_out=$blastn2NT_dir/${sample}.blastn2NT.blastn
    	blastn2NT_besthit=$blastn2NT_dir/${sample}.blastn2NT.besthit
		blastn2NT_filter_id=$blastn2NT_dir/${sample}.blastn2NT.filtered.id
		putative_virus_dir=$sampledir/06.putative_virus
		filtered_id=$putative_virus_dir/${sample}.blastx2NR-blastn2NT.filtered.id
		putative_virus_contigs_id=$putative_virus_dir/${sample}.putative_virus_contigs.id
		putative_virus_contigs_fna=$putative_virus_dir/${sample}.putative_virus_contigs.fna
		putative_virus_contigs_gene_fna=$putative_virus_dir/${sample}.putative_virus_contigs.gene.fna
		putative_virus_contigs_gene_faa=$putative_virus_dir/${sample}.putative_virus_contigs.gene.faa
		putative_virus_contigs_gene_gff3=$putative_virus_dir/${sample}.putative_virus_contigs.gene.gff3
		putative_virus_contigs_all_potential_genes=$putative_virus_dir/${sample}.putative_virus_contigs.all_potential_genes
		mkdir -p $putative_virus_dir
	
		##获取组装序列id
		echo "$seqkit fx2tab -ni $contigs_fna > $contigs_id" >> $runshell	
		##筛选出最优比对以及最优比对不是病毒的序列id
		echo "$python $besthit_filter_py $blastx2NR_out $blastx2NR_besthit $blastx2NR_filter_id" >> $runshell
		echo "$python $besthit_filter_py $blastn2NT_out $blastn2NT_besthit $blastn2NT_filter_id" >> $runshell
		echo "cat $blastx2NR_filter_id $blastn2NT_filter_id|sort|uniq > $filtered_id" >> $runshell
		##在所有序列中删除最优比对不是病毒的序列;得到的便是潜在的病毒序列(包含dark matter)
		echo "sort $contigs_id $filtered_id|uniq -u > $putative_virus_contigs_id" >> $runshell
		echo "$seqkit grep -f $putative_virus_contigs_id $contigs_fna > $putative_virus_contigs_fna" >> $runshell
		##翻译潜在的病毒序列
		echo "$prodigal -p meta -i $putative_virus_contigs_fna -d $putative_virus_contigs_gene_fna -a $putative_virus_contigs_gene_faa -f gff -o $putative_virus_contigs_gene_gff3 -s $putative_virus_contigs_all_potential_genes" >> $runshell

		echo -e "echo \"end at:\`date '+%F %T'\`\"" >> $runshell

		##运行
		if [ $qsub_status == 'yes' ];then
            cd $shelldir
			rm -f ${runshell}.*
			#qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:2 -l vf=4G,num_proc=2 $runshell
			qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:3 -l vf=20G,num_proc=3 $runshell
        else
            ls $runshell
        fi
	done
fi
