#! /bin/bash
#
# joint_vcf.sh <runfile> <log.in> <log.ou> <qsub>  #fix this in start.sh
# 
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 4 ]  #fix this depending on the actually needed parameters
then
        MSG="Parameter mismatch. Rerun as: $0 <runfile> <log.in> <log.ou> <qsub> "
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi
set +x
echo -e "\n\n#####################################################################################"        
echo -e "#############             BEGIN ANALYSIS PROCEDURE                    ###############"
echo -e "#####################################################################################\n\n"        

echo -e "\n\n#####################################################################################"        
echo -e "#############             DECLARING VARIABLES                         ###############"
echo -e "#####################################################################################\n\n"        
set -x

set -x
echo `date`
scriptfile=$0
runfile=$1
elog=$2
olog=$3
qsubfile=$4
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
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
gatkdir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
ref_local=${refdir}/$refgenome

outputdir=$rootdir
set +x
echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"            
echo -e "#######   we will need these guys throughout, let's take care of them now   ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"   
set -x 
DeliveryDir=$rootdir/$deliverydir/jointVCFs

qcfile=$rootdir/$deliverydir/docs/QC_test_results.txt            # name of the txt file with all QC test results
jointVCF=jointVCFcalled.vcf					# name of the resulting jointly called variants file
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

if [ ! -d $rootdir ]
then
    MSG="Invalid value specified for OUTPUTDIR=$rootdir in the runfile."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -d $DeliveryDir ]
then
    MSG="$DeliveryDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

############# This is actually what i need to do:
set +x
echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"            
echo -e "########### Joint genotyping script for all SAMPLEs: each 200 together          #####"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n" 
set -x 
cd $DeliveryDir

variantFiles=$( find ${outputdir} -name "*.GATKCombineGVCF.raw.vcf" | sed "s/^/ --variant /g" | tr "\n" " " )

if [ `expr ${#variantFiles}` -lt 1 ]
then
    MSG="no gvcf variant files were found at $rootdir"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

$javadir/java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar \
         -R $ref_local \
         -T  GenotypeGVCFs \
         $variantFiles  \
         -o $jointVCF

exitcode=$?
echo `date`

if [ $exitcode -ne 0 ]
then
         MSG="GATK GenotypeGVCFs  command failed exitcode=$exitcode."
         echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
         exit $exitcode;
fi

if [ ! -s $jointVCF ]
then
         echo -e "${SampleName}\tJointVARCALLING\tFAIL\tGATKGenotypeGVCFs produced an empty file $jointVCF\n" >> $qcfile
         MSG="joint-vcf produced an empty file $jointVCF exitcode=$exitcode."
         echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
         exit $exitcode;

fi

echo -e "Joint VARCALLING\tPASS\tAll analyses completed successfully for all samples" >> $qcfile
set +x
echo `date`
echo -e "\n\n##################################################################################"
echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############"
echo -e "##################################################################################\n\n"
set -x 

