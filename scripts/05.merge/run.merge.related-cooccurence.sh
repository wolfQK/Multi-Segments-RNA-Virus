samplelist=$1
poutdir=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California
related2clustering=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California/Analysis/clustering/shell/parser.author.ge500bp.cdhit-V2.1.py
cluster_file=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California/Analysis/clustering/Author/cd-hit/all_putative.author.ge500bp.contigs.cluster.clstr
overlap_cutoff=0.9
RdRp_list_file=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California/Analysis/results/new.author.putative.RdRp.list

cat $samplelist|while read sample
do
	sampledir=$poutdir/samples/$sample
	shelldir=$sampledir/shell
	runshell=$shelldir/step16.${sample}.related2clustering.sh
	echo -e "echo \"start at:\`date\` \"" > $runshell
	related_found_file=$sampledir/08.segments/${sample}.author-putative-viral-contigs.tblastx2ICTV_VH.virome.json
	outdir=$sampledir/08.segments/ICTV_VH
	echo "python $related2clustering --outdir $outdir --related_found_file $related_found_file --cluster_file $cluster_file --RdRp_list_file $RdRp_list_file --overlap_cutoff $overlap_cutoff" >> $runshell
	related_found_file=$sampledir/08.segments/${sample}.author-putative-viral-contigs.tblastx2NCBI_Virus.virome.json
	outdir=$sampledir/08.segments/NCBI_Virus
	echo "python $related2clustering --outdir $outdir --related_found_file $related_found_file --cluster_file $cluster_file --RdRp_list_file $RdRp_list_file --overlap_cutoff $overlap_cutoff" >> $runshell
	echo -e "echo \"end at:\`date\` \"" >> $runshell
	cd $shelldir
	#ls $runshell
	rm -fr ${runshell}.*
	qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:1 -l vf=2G,num_proc=1 $runshell
done
