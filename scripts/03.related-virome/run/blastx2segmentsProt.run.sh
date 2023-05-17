parser_blastx_py=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Multi-Segments-Virus/scripts/02.related-virome/blastx2segmentsProt.py
samplelist=$1
poutdir=$2
run=$3

cat $samplelist|while read sample
do
	sampledir=$poutdir/samples/$sample
	shelldir=$sampledir/shell
	runshell=$shelldir/step11.${sample}.parser.blastx2segments.sh
	echo -e "echo \"start at:\`date '+%F %T'\`\"" > $runshell
	
	segmentsdir=$sampledir/08.segments/author
    related_p=80
	blastx_out_file=$segmentsdir/${sample}.author.putative_viral_contigs.blastx2ICTV_VH.blastx
    RdRp_check_file=$sampledir/09.blastx2RdRp/author/${sample}.author-contigs.blastx2RdRp.parser.r100.txt
    out_file=$segmentsdir/${sample}.author.putative_viral_contigs.blastx2ICTV_VH.parser.txt
    virome_out_anyhit=$segmentsdir/${sample}.author.putative_viral_contigs.blastx2ICTV_VH.anyhit.json
    virome_out_besthit=$segmentsdir/${sample}.author.putative_viral_contigs.blastx2ICTV_VH.besthit.json
    
	echo "python $parser_blastx_py --related_p 80 --blastx_out_file $blastx_out_file --RdRp_check_file $RdRp_check_file --out_file $out_file  --virome_out_anyhit $virome_out_anyhit --virome_out_besthit $virome_out_besthit " >> $runshell

	echo -e "echo \"end at:\`date '+%F %T'\`\"" >> $runshell
    if [ $run == 'yes' ];then
	    cd $shelldir
	    rm -fr ${runshell}.*
	    #qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:1 -l vf=2G,num_proc=1 $runshell
        sh $runshell
    else
	    ls $runshell
    fi

done
