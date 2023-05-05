#poutdir=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California
python=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/gatk/3.7/bin/python
besthit_filter_py=$PWD/besthit_filter.py
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
		run2sample=$poutdir/../list/Run2libraryName.list	
		sampleName=`grep -w $sample $run2sample|awk '{print $2}'`
		sampledir=$poutdir/$sample
		shelldir=$sampledir/shell
		runshell=$shelldir/step10.${sample}.find.dark-contigs.sh
		##开始
		echo -e "echo \"start at:\`date '+%F %T'\`\"" > $runshell	
		
		## --- 自己组装的contigs --- ## 
		#未去宿主-metaSPAdes
		contigs_id=$sampledir/02.Assembly/metaSPAdes/contigs.id
		blastx2NR_dir=$sampledir/04.blastx2NR
		blastx2NR_besthit=$blastx2NR_dir/${sample}.blastx2NR.besthit
		blastx2NR_id=$blastx2NR_dir/${sample}.blastx2NR.id
		echo "awk '{print \$1}' $blastx2NR_besthit |sort|uniq > $blastx2NR_id" >> $runshell
		blastn2NT_dir=$sampledir/05.blastn2NT
    	blastn2NT_besthit=$blastn2NT_dir/${sample}.blastn2NT.besthit
		blastn2NT_id=$blastn2NT_dir/${sample}.blastn2NT.id
		echo "awk '{print \$1}' $blastn2NT_besthit |sort|uniq > $blastn2NT_id" >> $runshell
		putative_virus_dir=$sampledir/06.putative_virus
		blastId=$putative_virus_dir/blastId
		echo "sort $blastx2NR_id $blastn2NT_id|uniq > $blastId" >> $runshell
		darkcontigsId=$putative_virus_dir/${sample}.dark-contigs.id
		if [ -f $darkcontigsId ];then
            rm -f $darkcontigsId
        fi
        echo "sort $contigs_id $blastId|uniq -u|while read id" >> $runshell
        echo "do" >> $runshell
        echo -e "\techo \"$sampleName~\$id\" >> $darkcontigsId" >> $runshell
        echo "done" >> $runshell
        echo "rm -f $blastx2NR_id $blastn2NT_id $blastId" >> $runshell

		#去宿主-metaSPAdes
		contigs_id=$sampledir/02.Assembly/rmHost/metaSPAdes/contigs.id
        blastx2NR_besthit=$blastx2NR_dir/${sample}.rmHost.spades.contigs.blastx2NR.besthit
        blastx2NR_id=$blastx2NR_dir/${sample}.rmHost.spades.contigs.blastx2NR.id
		echo "awk '{print \$1}' $blastx2NR_besthit |sort|uniq > $blastx2NR_id" >> $runshell
        blastn2NT_besthit=$blastn2NT_dir/${sample}.rmHost.spades.contigs.blastn2NT.besthit
        blastn2NT_id=$blastn2NT_dir/${sample}.rmHost.spades.contigs.blastn2NT.id
		echo "awk '{print \$1}' $blastn2NT_besthit |sort|uniq > $blastn2NT_id" >> $runshell
        putative_virus_dir=$sampledir/06.putative_virus
        blastId=$putative_virus_dir/blastId
        echo "sort $blastx2NR_id $blastn2NT_id|uniq > $blastId" >> $runshell
        darkcontigsId=$putative_virus_dir/${sample}.rmHost.spades.dark-contigs.id
		if [ -f $darkcontigsId ];then
            rm -f $darkcontigsId
        fi
        echo "sort $contigs_id $blastId|uniq -u|while read id" >> $runshell
        echo "do" >> $runshell
        echo -e "\techo \"$sampleName~\$id\" >> $darkcontigsId" >> $runshell
        echo "done" >> $runshell
        echo "rm -f $blastx2NR_id $blastn2NT_id $blastId" >> $runshell

		## --- 作者的contigs --- ##
		author_contigs_dir=/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California/figshare/contigs
		contigs_id=$author_contigs_dir/$sampleName/contigs.id
		blastx2NR_besthit=$blastx2NR_dir/${sample}.author-contigs.blastx2NR.besthit
        blastx2NR_id=$blastx2NR_dir/${sample}.author-contigs.blastx2NR.id
		echo "awk '{print \$1}' $blastx2NR_besthit |sort|uniq > $blastx2NR_id" >> $runshell
		blastn2NT_besthit=$blastn2NT_dir/${sample}.author-contigs.blastn2NT.besthit
        blastn2NT_id=$blastn2NT_dir/${sample}.author-contigs.blastn2NT.id
		echo "awk '{print \$1}' $blastn2NT_besthit |sort|uniq > $blastn2NT_id" >> $runshell
        putative_virus_dir=$sampledir/06.putative_virus
        blastId=$putative_virus_dir/blastId
        echo "sort $blastx2NR_id $blastn2NT_id|uniq > $blastId" >> $runshell
        darkcontigsId=$putative_virus_dir/${sample}.author-contigs.dark-contigs.id
		if [ -f $darkcontigsId ];then
			rm -f $darkcontigsId
		fi
        echo "sort $contigs_id $blastId|uniq -u|while read id" >> $runshell
		echo "do" >> $runshell
		echo -e "\techo \"$sampleName~\$id\" >> $darkcontigsId" >> $runshell
		echo "done" >> $runshell
		echo "rm -f $blastx2NR_id $blastn2NT_id $blastId" >> $runshell


		##结束
		echo -e "echo \"end at:\`date '+%F %T'\`\"" >> $runshell
		##运行
		if [ $qsub_status == 'yes' ];then
            cd $shelldir
			rm -f ${runshell}.*
			qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:1 -l vf=2G,num_proc=1 $runshell
        else
            ls $runshell
        fi
	done
fi
