RdRp_check_py=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Multi-Segments-Virus/scripts/02.related-virome/RdRp.segmented.check.py
related=100
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
		runshell=$shelldir/step12.${sample}.blastx2RdRp.sh
		echo -e "echo \"start at:\`date '+ %F %T'\`\"" > $runshell
		blastx2RdRp_dir=$sampledir/09.blastx2RdRp
		#mkdir -p $blastx2RdRp_dir
		author_contigs_blastx2RdRp=$blastx2RdRp_dir/author/${sample}.author-contigs.blastx2RdRp.blastx
		parser_author_contigs_blastx2RdRp=$blastx2RdRp_dir/author/${sample}.author-contigs.blastx2RdRp.parser.r${related}.txt
		#echo "# $diamond blastx --threads 4 --query $author_contigs --db $RdRp_db --sensitive --max-target-seqs 50 --evalue 1E-5 --block-size 2.0 --index-chunks 1 --out $author_contigs_blastx2RdRp --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps" >> $runshell
		echo "python $RdRp_check_py --blastx2RdRp_out $author_contigs_blastx2RdRp --sample_out $parser_author_contigs_blastx2RdRp --related_p $related" >> $runshell
		
        rmHost_spades_contigs_blastx2RdRp=$blastx2RdRp_dir/rmHost/${sample}.rmHost-spades-contigs.blastx2RdRp.blastx
		parser_rmHost_spades_contigs_blastx2RdRp=$blastx2RdRp_dir/rmHost/${sample}.rmHost-spades-contigs.blastx2RdRp.parser.r${related}.txt
        #echo "# $diamond blastx --threads 4 --query $rmHost_spades_contigs --db $RdRp_db --sensitive --max-target-seqs 50 --evalue 1E-5 --block-size 2.0 --index-chunks 1 --out $rmHost_spades_contigs_blastx2RdRp --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps" >> $runshell 
		echo "python $RdRp_check_py --blastx2RdRp_out $rmHost_spades_contigs_blastx2RdRp --sample_out $parser_rmHost_spades_contigs_blastx2RdRp --related_p $related" >> $runshell

        withHost_spades_contigs_blastx2RdRp=$blastx2RdRp_dir/withHost/${sample}.withHost-spades-contigs.blastx2RdRp.blastx
		parser_withHost_spades_contigs_blastx2RdRp=$blastx2RdRp_dir/withHost/${sample}.withHost-spades-contigs.blastx2RdRp.parser.r${related}.txt
        #echo "# $diamond blastx --threads 4 --query $withHost_spades_contigs --db $RdRp_db --sensitive --max-target-seqs 50 --evalue 1E-5 --block-size 2.0 --index-chunks 1 --out $withHost_spades_contigs_blastx2RdRp --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps" >> $runshell 
		echo "python $RdRp_check_py --blastx2RdRp_out $withHost_spades_contigs_blastx2RdRp --sample_out $parser_withHost_spades_contigs_blastx2RdRp --related_p $related" >> $runshell


		echo -e "echo \"end at:\`date '+ %F %T'\`\"" >> $runshell
		##运行
        if [ $qsub_status == 'yes' ];then
            cd $shelldir
            rm -f ${runshell}.*
            # qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:1 -l vf=2G,num_proc=1 $runshell
            sh $runshell
            echo $sample
        else
            ls $runshell
        fi
	done
fi
#tmp_author_contigs_blastx2RdRp=$blastx2RdRp_dir/${sample}.author-contigs.blastx2RdRp.blastx.tmp
#mv $author_contigs_blastx2RdRp $tmp_author_contigs_blastx2RdRp
#sed '/-like/d' $tmp_author_contigs_blastx2RdRp > $author_contigs_blastx2RdRp
#rm -f $tmp_author_contigs_blastx2RdRp
#tmp_rmHost_spades_contigs_blastx2RdRp=$blastx2RdRp_dir/${sample}.rmHost-spades-contigs.blastx2RdRp.blastx.tmp
#mv $rmHost_spades_contigs_blastx2RdRp $tmp_rmHost_spades_contigs_blastx2RdRp
#sed '/-like/d' $tmp_rmHost_spades_contigs_blastx2RdRp > $rmHost_spades_contigs_blastx2RdRp
#rm -f $tmp_rmHost_spades_contigs_blastx2RdRp
