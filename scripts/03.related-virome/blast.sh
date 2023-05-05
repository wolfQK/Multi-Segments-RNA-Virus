# 输入contigs路径、输出路径，默认输出四种比对的结果；也可选择用哪种比对
BLASTBIN=/jdfssz1/ST_HEALTH/P20Z10200N0206/fengqikai/software/blast+2.10.1/blast-2.10.1+/bin
blastn=$BLASTBIN/blastn
tblastx=$BLASTBIN/tblastx
segmentsLIB=/jdfssz1/ST_HEALTH/P20Z10200N0206/fengqikai/biye/Multi-Segments-Virus/data/multi_segments_virus/ICTV_VH
ICTV_VH_segments_seq=$segmentsLIB/ICTV_VH.segments.nucl
diamond=/jdfssz1/ST_HEALTH/P20Z10200N0206/fengqikai/software/Diamond/2.0.9/diamond
ICTV_VH_segments_prot=$segmentsLIB/ICTV_VH.segments.prot

helpFunction()
{
    echo ""
    echo "Usage: $0 -l file_list -r qsub or not"
    echo -e "\t-l which 'file_list' is used to run BLAST2segments"
    echo -e "\t-r qsub or not (yes|no)"
    exit 1 # Exit script after printing help
}

while getopts "l:r:" opt
do
   case "$opt" in
        l ) file_list="$OPTARG" ;;
        r ) qsub_status="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

if [ -z $file_list ] || [ -z $qsub_status ];then
    echo "Some or all of the parameters are empty";
    helpFunction
else
    cat $file_list|while read file
    do
        sampledir=`dirname $file`
        shelldir=$sampledir/shell
		mkdir -p $shelldir
        sample=`basename $file|awk -F '.' '{print $1}'`
		runshell=$shelldir/step01.${sample}.blast2ICTV_VH.segments.nucl-prot.sh
        echo -e "echo \"start at:\`date\`\"" > $runshell
        blastn_out=$sampledir/${sample}.blastn2ICTV_VH.blastn
        blastp_out=$sampledir/${sample}.blastp2ICTV_VH.blastp
        blastx_out=$sampledir/${sample}.blastx2ICTV_VH.blastx
        tblastx_out=$sampledir/${sample}.tblastx2ICTV_VH.tblastx

        if [ -f $file ];then
            echo "$blastn -num_threads 3 -query $file -db $ICTV_VH_segments_seq -out $blastn_out -evalue 1E-5 -outfmt '6 qseqid qlen sseqid slen length pident mismatch gapopen gaps qcovs qstart qend sstart send evalue bitscore'" >> $runshell
            echo "$diamond blastp --threads 3 --query $file --db $ICTV_VH_segments_prot --out $blastp_out --sensitive --max-target-seqs 50 --evalue 1E-5 --block-size 2.0 --index-chunks 1  --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps" >> $runshell
            echo "$diamond blastx --threads 3 --query $file --db $ICTV_VH_segments_prot --out $blastx_out --sensitive --max-target-seqs 50 --evalue 1E-5 --block-size 2.0 --index-chunks 1 --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps" >> $runshell
            echo "$tblastx -num_threads 3 -query $file -db $ICTV_VH_segments_seq -out $tblastx_out -evalue 1E-5 -outfmt '6 qseqid qlen sseqid slen length pident mismatch gapopen gaps qcovs qstart qend sstart send evalue bitscore'" >> $runshell
            echo -e "echo \"end at:\`date\`\"" >> $runshell
            if [ $qsub_status == 'yes' ];then
                cd $shelldir
                rm -f ${runshell}.* core.*
                qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:3 -l vf=3G,num_proc=3 $runshell
                # qsub -clear -cwd -q st_supermem.q -P P20Z10200N0206 -binding linear:1 -l vf=1G,num_proc=1 $runshell
            else
                ls $runshell
                rm -f core.*
            fi
        else
            echo 'file do not exist !'
	    fi
    done
fi
