#!/bin/bash
#
# realign_varcall_by_chr.sh <runfile> <sample> <chr> <log.in> <log.ou> <qsub>
# 
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 6 ]
then
        MSG="Parameter mismatch. Rerun as: $0 <runfile> <sample> <chr> <log.in> <log.ou> <qsub> "
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi

echo -e "\n\n#####################################################################################"        
echo -e "#############             BEGIN ANALYSIS PROCEDURE                    ###############"
echo -e "#####################################################################################\n\n"        

echo -e "\n\n#####################################################################################"        
echo -e "#############             DECLARING VARIABLES                         ###############"
echo -e "#####################################################################################\n\n"        

umask 0027
set -x
echo `date`
scriptfile=$0
runfile=$1
SampleName=$2
chr=$3
elog=$4
olog=$5
qsubfile=$6
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


if [ ! -s $runfile ]
then
    MSG="$runfile configuration file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "Variant Calling Workflow failure message" "$redmine"
    exit 1;
fi

reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
sampledir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )        
aligner_mod=$( cat $runfile | grep -w ALIGNERMODULE | cut -d '=' -f2  )
aligner_parms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 )
picard_mod=$( cat $runfile | grep -w PICARDMODULE | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
samtools_mod=$( cat $runfile | grep -w SAMTOOLSMODULE | cut -d '=' -f2 )
samblaster_mod=$( cat $runfile | grep -w SAMBLASTERMODULE | cut -d '=' -f2 )
sorttool_mod=$( cat $runfile | grep -w SORTMODULE | cut -d '=' -f2 )
markduplicates=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
gatk_mod=$( cat $runfile | grep -w GATKMODULE | cut -d '=' -f2 ) 
gatk_dir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 ) 
sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
ref_local=${refdir}/$refgenome
dbsnp_local=${refdir}/$dbSNP
indel_local=${refdir}/$indeldir

outputdir=$rootdir/$SampleName

echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   we will need these guys throughout, let's take care of them now   ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"          

module load $samtools_mod
module load $sorttool_mod
module load $gatk_mod

SampleDir=$outputdir
AlignDir=$outputdir/align
RealignDir=$outputdir/realign
VarcallDir=$outputdir/variant
DeliveryDir=$rootdir/$deliverydir/$SampleName

qcfile=$rootdir/$deliverydir/docs/QC_test_results.txt            # name of the txt file with all QC test results
inputbam=$AlignDir/${SampleName}.wdups.sorted.bam                # name of the bam file that align-dedup produced
dedupsortedbam=${SampleName}.$chr.wdups.sorted.bam               # name of the aligned file -one chr out of inputbam
realignedbam=${SampleName}.$chr.realigned.bam                    # name of the realigned file
recalibratedbam=${SampleName}.$chr.recalibrated.bam              # name of the recalibrated file
rawvariant=${SampleName}.$chr.raw.g.vcf                          # name of the raw variant file


echo -e "\n\n##################################################################################" 
echo -e "##################################################################################"        
echo -e "#############                       SANITY CHECK                   ###############"
echo -e "##################################################################################"
echo -e "##################################################################################\n\n"

if [ ! -d $tmpdir ]
then
    mkdir -p $tmpdir
fi

if [ `expr ${#SampleName}` -lt 1]
then
    MSG="$SampleName sample undefined variable"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1     
else
    sID=$SampleName
    sPU=$SampleName
    sSM=$SampleName
fi
if [ `expr ${#sLB}` -lt 1 -o `expr ${#sPL}` -lt 1 -o `expr ${#sCN}` -lt 1 ] 
then
    MSG="SAMPLELB=$sLB SAMPLEPL=$sPL SAMPLECN=$sCN at least one of these fields has invalid values. "
    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )
rgheader=$( echo -n -e "@RG\t" )$( echo -e "${RGparms}"  | tr ":" "\t" | tr "=" ":" )

if [ ! -d $rootdir ]
then
    MSG="Invalid value specified for OUTPUTDIR=$rootdir in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -d $AlignDir ]
then
    MSG="$AlignDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -d $RealignDir ]
then
    MSG="$RealignDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi
if [ ! -d $VarcallDir ]
then
    MSG="$VarcallDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi
if [ ! -d $DeliveryDir ]
then
    MSG="$DeliveryDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -s $inputbam ]
then
    MSG="$inputbam aligned-dedup file not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi


echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   REALIGN-RECALIBRATION BLOCK STARTS HERE   $SampleName $chr        ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n" 



echo -e "\n\n##################################################################################"
echo -e "########### PREP WORK: split aligned.bam by chr and grab indels for chr            ###" 
echo -e "##################################################################################\n\n"

echo -e "\n### ploidy variable, its value is 2 for all chr except mitochondrial               ###\n"

if [ $chr == "M" ]
then
    ploidy=1
else
    ploidy=2
fi

echo -e "\n### split aligned.bam by region with  $inputbam and chr=$chr                       ###\n"     

cd $RealignDir

samtools view -bu -@ $thr -h $inputbam $chr > $dedupsortedbam

exitcode=$?

echo `date`
if [ $exitcode -ne 0 ]
then
	 MSG="samtools command failed to split $inputbam by chr exitcode=$exitcode"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	 exit $exitcode;
fi
if [ -s $dedupsortedbam ]
then     
    echo -e "### the file was created. But we are not done.     #############"
    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
    numAlignments=$( samtools view -c $dedupsortedbam ) 

    echo `date`
    if [ $numAlignments -eq 0 ]
    then
	echo -e "${SampleName}\tREALIGNMENT\tFAIL\tsamtools command did not produce alignments for $dedupsortedbam\n" >> $qcfile	    
	MSG="samtools command did not produce alignments for $dedupsortedbam"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
    else
	echo -e "####### $dedupsortedbam seems to be in order ###########"
    fi
else 
    MSG="samtools command did not produce a file $dedupsortedbam"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"        
    exit 1;          
fi	

echo -e "\n### index  $dedupsortedbam                                 ###\n"  

samtools index $dedupsortedbam

exitcode=$?

echo `date`
if [ $exitcode -ne 0 ]
then
	 MSG="samtools index command failed to split $inputbam by chr exitcode=$exitcode"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	 exit $exitcode;
fi

echo -e "\n\n### grab indels for this region. We need two variables one for realignment and one for recalibration"
echo -e "###  because they are specified differently for each GATK tool                            ##############\n\n"

cd $indel_local


recalparms2=$( find ${PWD} -name "${chr}.*.vcf" | sed "s/^/ --knownSites /g" | tr "\n" " " )
recalparms1=" -knownSites $dbsnp_local "
realparms=$( find ${PWD} -name "${chr}.*.vcf" | sed "s/^/ -known /g" | tr "\n" " " )

if [ `expr ${#realparms}` -lt 1 ]
then
    MSG="no indels were found for $chr in this folder $indel_local"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ `expr ${#recalparms2}` -lt 1 ]
then
    MSG="no indels were found for $chr in this folder $indel_local"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi


echo -e "\n\n##################################################################################"
echo -e "########### command one: executing GATK RealignerTargetCreator using known indels ####" 
echo -e "##################################################################################\n\n"

cd $RealignDir

java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatk_dir/GenomeAnalysisTK.jar \
	 -R $ref_local \
	 -I $dedupsortedbam \
	 -T RealignerTargetCreator \
	 -nt $thr \
	 $realparms \
	 -o ${SampleName}.$chr.realignTargetCreator.intervals

exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
	 MSG="RealignerTargetCreator command failed exitcode=$exitcode. realignment for sample $SampleName.$chr. stopped"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	 exit $exitcode;
fi
if [ ! -s ${SampleName}.$chr.realignTargetCreator.intervals ]
then
	 echo -e "${SampleName}\tREALIGNMENT\tWARN\t${SampleName}.$chr.RealignTargetCreator.intervals is an empty file. Skipping Indel realignment cmd\n" >> $qcfile
         ln -s $dedupsortedbam $RealignDir/$realignedbam
else 

	echo -e "\n\n##################################################################################" 
	echo -e "########### command two: executing GATK IndelRealigner and generating BAM    #####"
	echo -e "##################################################################################\n\n"

	java -Xmx8g -Djava.io.tmpdir=$tmpdir -jar $gatk_dir/GenomeAnalysisTK.jar \
                  -R $ref_local -I $dedupsortedbam -T IndelRealigner $realparms \
                  --targetIntervals ${SampleName}.$chr.realignTargetCreator.intervals \
                  -o $realignedbam

	exitcode=$?

	echo -e "\n\n##################################################################################" 
	echo -e "########### command three: sanity check for GATK IndelRealigner                  #####"
	echo -e "##################################################################################\n\n"

	echo `date`
	if [ $exitcode -ne 0 ]
	then
		 MSG="indelrealigner command failed exitcode=$exitcode. realignment for sample $SampleName stopped"
		 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		 exit $exitcode;
	fi

	if [ -s $realignedbam ]
	then     
	    echo -e "### the file was created. But we are not done.     #############"
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
	    numAlignments=$( samtools view -c $realignedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tREALIGNMENT\tFAIL\tGATK IndelRealigner command did not produce alignments for $realignedbam\n" >> $qcfile	    
		MSG="GATK IndelRealigner command did not produce alignments for $realignedbam"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    else
		echo -e "####### $realignedbam seems to be in order ###########"
	    fi
	else 
	    MSG="indelrealigner command did not produce a file $realignedbam"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"        
	    exit 1;          
	fi	
fi

echo `date`     

echo -e "\n\n##################################################################################" 
echo -e "########### command four: run GATK BaseRecalibrator using known indels and dbSNP    ##"
echo -e "##################################################################################\n\n"

java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatk_dir/GenomeAnalysisTK.jar \
         -T BaseRecalibrator \
         -R $ref_local \
         -I $realignedbam \
         $recalparms1 \
         --out $SampleName.$chr.recal_report.grp \
         -nct $thr 

exitcode=$?

echo `date`
if [ $exitcode -ne 0 ]
then
	MSG="BaseRecalibrator command failed exitcode=$exitcode. recalibration for sample $SampleName stopped"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	exit $exitcode;
fi
if [ ! -s $SampleName.$chr.recal_report.grp ]
then
	echo -e "${SampleName}\tRECALIBRATION\tWARN\t$SampleName.$chr.recal_report.grp is an empty file. Skipping recalibration cmd\n" >> $qcfile
        ln -s $dedupsortedbam $RealignDir/$recalibratedbam
else
	echo -e "\n\n##################################################################################" 
	echo -e "########### command five: GATK BQSR step                                         #####"
	echo -e "##################################################################################\n\n"


        java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatk_dir/GenomeAnalysisTK.jar \
                -R $ref_local \
                -I $realignedbam \
                -T PrintReads \
                -BQSR $SampleName.$chr.recal_report.grp \
                --out $recalibratedbam \
                -nct $thr

        exitcode=$?


	echo -e "\n\n##################################################################################" 
	echo -e "########### command six: sanity check for GATK BQSR step                        #####"
	echo -e "##################################################################################\n\n"


	if [ -s $recalibratedbam ]
	then     
	    echo -e "### the file was created. But we are not done.     #############"
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
	    numAlignments=$( samtools view -c $recalibratedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
                echo -e "${SampleName}\tRECALIBRATION\tFAIL\tBQSR Recalibrator command did not produce alignments for $recalibratedbam\n" >> $qcfile	    
		MSG="GATK BQSR Recalibrator command did not produce alignments for $recalibratedbam"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    else
		echo -e "####### $realignedbam seems to be in order ###########"
	    fi
	else 
	    MSG="GATK BQSR Recalibrator command did not produce a file $recalibratedbam"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"        
	    exit 1;          
	fi	

fi      
echo `date` 

echo -e "\n\n##################################################################################"
echo -e "#############     END REALIGN-RECALIBRATION BLOCK                         ############"
echo -e "##################################################################################\n\n"          


echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"	
echo -e "##################################################################################"        
echo -e "#############    GATK VARIANT CALLING   FOR SAMPLE $SampleName  $chr   ###########"
echo -e "##################################################################################"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"

echo `date`        

cd $VarcallDir

echo -e "\n\n##################################################################################"            
echo -e "########### command one: executing GATK HaplotypeCaller command         ##########" 
echo -e "##################################################################################\n\n"



java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatk_dir/GenomeAnalysisTK.jar \
	 -T HaplotypeCaller \
	 -R $ref_local \
	 --dbsnp $dbsnp_local \
	 -I $RealignDir/$recalibratedbam \
	 --emitRefConfidence GVCF \
	 -gt_mode DISCOVERY \
	 -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
	 -stand_call_conf 30 \
	 -stand_emit_conf 30 \
	 --sample_ploidy $ploidy \
	 -nt 1 -nct $thr \
	 -o $rawvariant


exitcode=$?
echo `date`

if [ $exitcode -ne 0 ]
then
	MSG="GATK HaplotypeCaller command failed exitcode=$exitcode for $rawvariant"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	exit $exitcode;
fi
if [ ! -s $rawvariant ]
then
	echo -e "${SampleName}\tVCALL\tFAIL\tHaplotypeCaller command did not produce results for $rawvariant\n" >> $qcfile	    
	MSG="GATK HaplotypeCaller command did not produce results for $rawvariant"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

echo -e "\n\n##################################################################################"
echo -e "#############       END VARIANT CALLING BLOCK                         ############"        
echo -e "##################################################################################\n\n"


echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"		
echo -e "##################################################################################"        
echo -e "#############   WRAP UP                                               ############"        
echo -e "##################################################################################"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"	

module unload $samtools_mod
module unload $sorttol_mod
module unload $gatk_mod
echo `date`
 
#cp $RealignDir/${SampleName}.$chr.recalibrated.bam   $DeliveryDir

# we will merge all variant files for this sample and copy that file to delivery
#cp $VarcallDir/rawvariant=${SampleName}.$chr.raw.vcf $DeliveryDir  

echo `date`
echo -e "\n\n##################################################################################"
echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############"
echo -e "##################################################################################\n\n"


