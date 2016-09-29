#!/bin/bash
#
# merge_vcf.sh <runfile> <sample> <log.in> <log.ou> <qsub>
# 
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch. Rerun as: $0 <runfile> <sample> <log.in> <log.ou> <qsub> "
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi
set -x
echo -e "\n\n#####################################################################################"        
echo -e "#############             BEGIN ANALYSIS PROCEDURE                    ###############"
echo -e "#####################################################################################\n\n"        

echo -e "\n\n#####################################################################################"        
echo -e "#############             DECLARING VARIABLES                         ###############"
echo -e "#####################################################################################\n\n"        
set -x
echo `date`
scriptfile=$0
runfile=$1
SampleName=$2
elog=$3
olog=$4
qsubfile=$5
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


if [ ! -s $runfile ]
then
    MSG="$runfile runfile not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "Variant Calling Workflow failure message" "$redmine"
    exit 1;
fi

reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 ) 
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
gatkdir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
novocraftdir=$( cat $runfile | grep -w NOVOCRAFTDIR | cut -d '=' -f2  )
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
ref_local=${refdir}/$refgenome
dbsnp_local=${refdir}/$dbSNP
indel_local=${refdir}/$indeldir
indices=$( cat $runfile | grep -w CHRNAMES | cut -d '=' -f2 | tr ':' ' ' )
outputdir=$rootdir/$SampleName

set +x
echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   we will need these guys throughout, let's take care of them now   ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"          
set -x 

SampleDir=$outputdir
RealignDir=$outputdir/realign
VarcallDir=$outputdir/variant
DeliveryDir=$rootdir/$deliverydir/$SampleName

qcfile=$rootdir/$deliverydir/docs/QC_test_results.txt            # name of the txt file with all QC test results
tmpvariant=${SampleName}.raw.vcf                                 # name of raw variant file pre-sorting
rawvariant=${SampleName}.GATKCombineGVCF.raw.vcf                 # name of the raw variant file
outbam=${SampleName}.recalibrated.bam                            # name of the ready-for-analysis bam file

set +x
echo -e "\n\n##################################################################################" 
echo -e "##################################################################################"        
echo -e "#############                       SANITY CHECK                   ###############"
echo -e "##################################################################################"
echo -e "##################################################################################\n\n"
set -x

if [ ! -d $tmpdir ]
then
    mkdir -p $tmpdir
fi

if [ `expr ${#SampleName}` -lt 1 ]
then
    MSG="$SampleName sample undefined variable"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1     
fi


if [ ! -d $rootdir ]
then
    MSG="Invalid value specified for OUTPUTDIR=$rootdir in the runfile."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -d $VarcallDir ]
then
    MSG="$VarcallDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi
if [ ! -d $RealignDir ]
then
    MSG="$RealignDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi
if [ ! -d $DeliveryDir ]
then
    MSG="$DeliveryDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi


set +x
echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   MERGE BAMS BLOCK STARTS HERE  FOR              $SampleName        ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n" 
set -x
cd $RealignDir
set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command one: gather all files to be merged             #####"
echo -e "##################################################################################\n\n"
set -x 
chr_bamList=$( ls -1 ${SampleName}.*.recalibrated.bam | tr "\n" " " )

if [ `expr ${#chr_bamList}` -lt 1 ]
then
    MSG="no bam files were found at $RealignDir"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi
set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command two: merge files with novosort              #####"
echo -e "##################################################################################\n\n"
set -x
$novocraftdir/novosort --index --threads $thr --tmpdir $tmpdir -o $outbam  ${SampleName}.*.recalibrated.bam 

exitcode=$?
set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command three: sanity check for novosort                        #####"
echo -e "##################################################################################\n\n"
set -x 

echo `date`
if [ $exitcode -ne 0 ]
then
	 MSG="novosort command failed exitcode=$exitcode. merge for sample $SampleName stopped"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	 exit $exitcode;
fi
if [ -s "$outbam" ]
then     
    echo -e "### the file was created. But we are not done.     #############"
    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
    numAlignments=$( $samtoolsdir/samtools view -c $outbam ) 

    echo `date`
    if [ $numAlignments -eq 0 ]
    then
	echo -e "${SampleName}\tMERGE\tFAIL\tnovosort command did not produce alignments for $outbam\n" >> $qcfile	    
	MSG="novosort command did not produce alignments for $outbam"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
    else
	echo -e "####### $outbam seems to be in order ###########"
    fi
else 
    MSG="novosort command did not produce a file $outbam"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"        
    exit 1;          
fi	

set +x
echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   MERGE VCFs BLOCK STARTS HERE  FOR              $SampleName        ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n" 
set -x
cd $VarcallDir

set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command one: checking that there are vcf files to merge          #####"
echo -e "##################################################################################\n\n"
set -x 
chr_vcfList=$( ls -1 ${SampleName}.*.vcf | tr "\n" " " )

if [ `expr ${#chr_vcfList}` -lt 1 ]
then
    MSG="no variant files were found at $VarcallDir"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi
set +x
echo -e "\n\n##################################################################################" 
echo -e "########### STAGE two: prepare vcfs and create ordered list of vcfs to merge   #####"
echo -e "##################################################################################\n\n"
set -x

ordered_vcfs=""

for chr in $indices
do
	set +x
	echo -e "\n\n##################################################################################"  
	echo -e "##################################################################################"          	
	echo -e "#######   Processing vcf for chr=$chr                                       ######"
	echo -e "##################################################################################"  
	echo -e "##################################################################################\n\n"
	set -x
        ### the vcf file for this chr and  GATK-CombineGVCFs
        
        thisvcf=$( ls -1 ${SampleName}.$chr.*.vcf | sed 's/^/ --variant /' | sed 's/\n/ /' )
        
        ### the vcf file for this chr and  picard-SortVcf
        
        #thisvcf=$( ls -1 ${SampleName}.$chr.*.vcf | sed 's/^/ INPUT=/' | sed 's/\n/ /' )

        ### now we append this name to the ordered list
        
        ordered_vcfs=${ordered_vcfs}$thisvcf
done

set +x
echo -e "\n\n#### outside of loop over indices   \n\n" 
set -x


if [ `expr ${#ordered_vcfs}` -lt 1 ]
then
    MSG="no variant files were found at $VarcallDir"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi
set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command three: merge with picard-SortVcf using ordered list of vcfs   #####"
echo -e "##################################################################################\n\n"
set -x
#module load $java_mod
#module load picard_mod

#java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $picardir/picard.jar SortVcf  OUTPUT=$tmpvariant  $ordered_vcfs

#exitcode=$?
#echo `date`
#if [ $exitcode -ne 0 ]
#then
#	 MSG="picard SortVcf command failed exitcode=$exitcode."
#	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
#	 exit $exitcode;
#fi
#if [ ! -s $tmpvariant ]
#then
#	 echo -e "${SampleName}\tVARCALLING\tFAIL\tpicard SortVcf produced an empty file $tmpvariant\n" >> $qcfile
#	 MSG="picard SortVcf produced an empty file $tmpvariant exitcode=$exitcode."
#	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
#	 exit $exitcode;
#
#fi

#vcf-sort $tmpvariant > $rawvariant
    
#exitcode=$?
#echo `date`
#if [ $exitcode -ne 0 ]
#then
#	 MSG="Cvcf-sort  command failed exitcode=$exitcode."
#	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
#	 exit $exitcode;
#fi

set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command three: merge with GATK-CombineGVCFs                          #####"
echo -e "##################################################################################\n\n"
set -x 

$javadir/java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar \
	 -R $ref_local \
	 --dbsnp $dbsnp_local \
	 $ordered_vcfs  \
	 -T CombineGVCFs \
	 -o $rawvariant 
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
	 MSG="GATK CombineGVCFs  command failed exitcode=$exitcode."
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	 exit $exitcode;
fi
if [ ! -s $rawvariant ]
then
	 echo -e "${SampleName}\tVARCALLING\tFAIL\tGATKCombineGVCFs produced an empty file $rawvariant\n" >> $qcfile
	 MSG="vcf-sort produced an empty file $rawvariant exitcode=$exitcode."
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	 exit $exitcode;

fi 

set +x
echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"		
echo -e "##################################################################################"        
echo -e "#############   COPY RESULTS TO DELIVERY and also to                              ############"        
echo -e "##################################################################################"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"	
set -x 

echo `date`

cp $RealignDir/$outbam       $DeliveryDir

if [ ! -s $DeliveryDir/$outbam ]
then

    MSG="copy failed from $DeliveryDir/$outbam to $DeliveryDir/$outbam"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi

cp $VarcallDir/$rawvariant   $DeliveryDir

if [ ! -s $DeliveryDir/$rawvariant ]
then

    MSG="copy failed from $DeliveryDir/$rawvariant to $DeliveryDir/$rawvariant"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi

echo -e "${SampleName}\tVARCALLING\tPASS\tAll analyses completed successfully for this sample" >> $qcfile
set +x
echo `date`
echo -e "\n\n##################################################################################"
echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############"
echo -e "##################################################################################\n\n"
set -x 

