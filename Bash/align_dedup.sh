#!/bin/bash
#
# align_dedup.sh <runfile> <SampleName> <read1> <read2> <reportticket> 
# 
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu

set -x
if [ $# != 4 ]
then
        MSG="Parameter mismatch. Rerun as: $0 <runfile> <SampleName> <read1> <read2> "
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi

set +x
echo -e "\n\n#####################################################################################" >&2       
echo -e "#############             BEGIN ANALYSIS PROCEDURE                    ###############" >&2
echo -e "#####################################################################################\n\n" >&2

echo -e "\n\n#####################################################################################" >&2        
echo -e "#############             DECLARING VARIABLES                         ###############" >&2
echo -e "#####################################################################################\n\n" >&2
set -x        

umask 0027
echo `date`
scriptfile=$0
runfile=$1
SampleName=$2
R1=$3
R2=$4

if [ ! -s $runfile ]
then
    MSG="$runfile runfile not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "Variant Calling Workflow failure message" "$redmine"
    exit 1;
fi

set +x
echo -e "\n\n##################################################################################" >&2 
echo -e "##################################################################################" >&2          	
echo -e "#######   we will need these guys throughout, let's take care of them now   ######" >&2
echo -e "##################################################################################" >&2 
echo -e "##################################################################################\n\n" >&2
set -x          

#######################################################
# remove these parts, as they are not needed. just for testing!
# just create the corresponding inputs!

outputdir=$rootdir/$SampleName
deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )

refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
alignertool=$( cat $runfile | grep -w ALIGNERTOOL | cut -d '=' -f2  )
bwamemdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2  )
bwamem_parms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 )
bwa_index=$( cat $runfile | grep -w BWAINDEX | cut -d '=' -f2 )
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
dup_cutoff=$( cat $runfile | grep -w  DUP_CUTOFF | cut -d '=' -f2 )
map_cutoff=$( cat $runfile | grep -w  MAP_CUTOFF | cut -d '=' -f2 )
outputdir=$rootdir/$SampleName

######################################################
SampleDir=$outputdir
AlignDir=$outputdir/align
RealignDir=$outputdir/realign
VarcallDir=$outputdir/variant
DeliveryDir=$rootdir/$deliverydir/$SampleName

qcfile=$rootdir/$deliverydir/docs/QC_test_results.txt            # name of the txt file with all QC test results
alignedbam=${SampleName}.nodups.bam                              # name of the aligned bam
alignedsortedbam=${SampleName}.nodups.sorted.bam                 # name of the aligned-sorted bam
dedupbam=${SampleName}.wdups.bam                                 # name of the deduplicated bam
dedupsortedbam=${SampleName}.wdups.sorted.bam                    # name of the dedup-sorted bam (output of this module)

set +x
echo -e "\n\n##################################################################################" >&2      
echo -e "#############                       SANITY CHECK                   ###############" >&2
echo -e "##################################################################################\n\n" >&2
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

if [ ! -d $outputdir ]
then
    MSG="$outputdir outputdir not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -d $DeliveryDir ]
then
    mkdir -p $DeliveryDir
fi
if [ `expr ${#R1}` -lt 1]
then
    MSG="$R1 read one file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1
elif [ ! -s $R1 ]
then
    MSG="$R1 read one file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1                
fi

if [ `expr ${#R2}` -lt 1]
then
    MSG="$R2 read two file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1
elif [ ! -s $R2 ]
then
    MSG="$R2 read two  file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1                
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


if [ `expr ${#markduplicates}` -lt 1 ]
then
    markduplicates="NOVOSORT"
fi

        
set +x
echo -e "\n\n##################################################################################" >&2 
echo -e "##################################################################################" >&2          
echo -e "##################################################################################" >&2        
echo -e "#############   ALIGN-DEDUPPLICATION  FOR SAMPLE $SampleName       ###############" >&2
echo -e "##################################################################################" >&2
echo -e "##################################################################################" >&2 
echo -e "##################################################################################\n\n" >&2          
set -x 

echo `date` 

cd  $AlignDir

set +x
echo -e "\n\n##################################################################################" >&2
echo -e "#############   select the dedup tool and then run the command        ############" >&2
echo -e "#############   choices:  SAMBLASTER or NOVOSORT                      ############" >&2
echo -e "##################################################################################\n\n" >&2          
set -x


if [ $markduplicates == "SAMBLASTER" ]
then
        set +x
	echo -e "\n\n##################################################################################" >&2
	echo -e "##CASE1: dedup tool is $markduplciates we use a single command for align-deduplication ##" >&2 
	echo -e "##################################################################################\n\n" >&2

	echo -e "\n\n##################################################################################" >&2
	echo -e "############# step one: alignment and deduplication                ############" >&2	     
	echo -e "##################################################################################\n\n" >&2
	set -x
	

        $bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index $R1 $R2 | $samblaster | $samtoolsdir/samtools view -@ $thr -bSu -> $dedupbam 
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment-dedup step  failed for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"| mail -s "[Task #${reportticket}]" "$redmine,$email"
	    exit $exitcode;
	fi
	        
	set +x
	echo -e "\n\n##################################################################################" >&2	     
	echo -e "#############  step two: making sure that a file was produced with alignments    #####" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x

	if [ -s $AlignDir/$dedupbam ]
	then 
	    set +x 		    
	    echo -e "### the file was created. But we are not done.     #############" >&2
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
	    set -x
	    numAlignments=$( $samtoolsdir/samtools view -c $AlignDir/$dedupbam ) 	
	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for $AlignDir/$dedupbam\n" >> $qcfile	    
		MSG="bwa mem command did not produce alignments for $AlignDir/$dedupbam alignment failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    else
		set +x
		echo -e "####### $AlignDir/$dedupbam seems to be in order ###########" >&2
		set -x
	    fi
	else 
	    MSG="bwa mem command did not produce a file $AlignDir/$dedupbam alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi       

	set +x
	echo -e "\n\n##################################################################################" >&2
	echo -e "#############  step three: sort                                      ############" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x

	$novocraftdir/novosort --index --tmpdir $tmpdir --threads $thr --compression 1 -o $dedupsortedbam $dedupbam
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="align-sorting step failed for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"       
	    exit $exitcode;
	fi
	
	set +x
	echo -e "\n\n##################################################################################" >&2	     
	echo -e "#############  step four: making sure that a file was produced with alignments    #####" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x

	if [ -s $AlignDir/$dedupsortedbam ]
	then 
	    set +x	    
	    echo -e "### the file was created. But we are not done.     #############" >&2
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
	    set -x 
	    numAlignments=$( $samtoolsdir/samtools view -c $AlignDir/$dedupsortedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tnovosort command did not produce a file for $AlignDir/$dedupsortedbam\n" >> $qcfile	    
		MSG="novosort command did not produce a file for $AlignDir/$dedupsortedbam"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    else
		set +x
		echo -e "####### $AlignDir/$dedupbam seems to be in order ###########" >&2
	        set -x
       
        fi
	else 
	    MSG="novosort command did not produce a file $AlignDir/$dedupsortedbam"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi       	

	set +x		
	echo -e "\n\n##################################################################################" >&2
	echo -e "#############      END SAMBLASTER BLOCK                               ############" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x             

elif  [ $markduplicates == "NOVOSORT" ]
then
	set +x
	echo -e "\n\n##################################################################################" >&2
	echo -e "CASE2: dedup tool is NOVOSORT. one cmd for align and one for dedup-sort   ########" >&2  
	echo -e "##################################################################################\n\n" >&2

	echo -e "\n\n##################################################################################" >&2	     
	echo -e "#############  step one: alignment                                 ############" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x 
        

        if [ $alignertool== "BWA" ]
        then


	   $bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index $R1 $R2 | $samtoolsdir/samtools view -@ $thr -bSu -> $alignedbam 
	   exitcode=$?
	   echo `date`
	   if [ $exitcode -ne 0 ]
	   then
	       MSG="alignment step  failed for sample $SampleName exitcode=$exitcode."
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	       exit $exitcode;
	   fi
        elif [ $alignertool == "NOVOALIGN" ]
	then
           $novocraftdir/novoalign $novoalign_parms  -c $thr -d ${refdir}/${novoalign_index} -f $R1 $R2 | $samtoolsdir/samtools view -@ $thr -bS - > $alignedbam
           exitcode=$?
           echo `date`
           if [ $exitcode -ne 0 ]
           then
               MSG="alignment step  failed for sample $SampleName exitcode=$exitcode."
               echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
               exit $exitcode;
           fi
        fi   
	
	set +x 
	echo -e "\n\n##################################################################################" >&2	     
	echo -e "#############  step two:making sure that a file was produced with alignments     #####" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x 

	if [ -s $AlignDir/$alignedbam ]
	then
	    set +x            
	    echo -e "### the file was created. But we are not done.     #############" >&2
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
 	    set -x 		

	    numAlignments=$( $samtoolsdir/samtools view -c $AlignDir/$alignedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\taligner command did not produce alignments for $AlignDir/$alignedbam\n" >> $qcfile	    
		MSG="aligner command did not produce alignments for $AlignDir/$alignedbam alignment failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"      
		exit 1;
	    else
		set +x
		echo -e "####### $AlignDir/$alignedbam seems to be in order ###########" >&2
		set -x 
	    fi
	else 
	    MSG="aligner command did not produce a file $AlignDir/$alignedbam alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi 
	set +x      
	echo -e "\n\n##################################################################################" >&2	     
	echo -e "#############  step three: sort + dedup + indexing                       ############" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x

	$novocraftdir/novosort -r "${rgheader}" --markDuplicates  -t $tmpdir -c $thr -i -o $dedupsortedbam $alignedbam

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed during sorting-deduplication for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit $exitcode
	fi
 
	set +x	
	echo -e "\n\n##################################################################################" >&2	     
	echo -e "#############  step four: making sure that a file was produced with alignments #######" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x
	
	if [ -s $AlignDir/$dedupsortedbam ]
	then
	    set +x	
	    echo -e "### the file was created. But we are not done.     #############" >&2
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
	    set -x 
	    
	    numAlignments=$( $samtoolsdir/samtools view -c $AlignDir/$dedupsortedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tnovosort command did not produce a file for $AlignDir/$dedupsortedbam\n" >> $qcfile	    
		MSG="novosort command did not produce a file for $AlignDir/$dedupsortedbam"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		exit 1;
	    else
		echo -e "####### $AlignDir/$dedupbam seems to be in order ###########"
	    fi
	else 
	    MSG="novosort command did not produce a file $AlignDir/$dedupsortedbam"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi   
	
	echo -e "\n\n##################################################################################"
	echo -e "#############      END NOVOSRT  BLOCK                                 ############"
	echo -e "##################################################################################\n\n"             
	

elif  [ $markduplicates == "PICARD" ]
then
	set +x
	echo -e "\n\n##################################################################################" >&2
	echo -e "CASE2: dedup tool is PICARD. one cmd for align, one for sort, one for dedup   ########" >&2 
	echo -e "##################################################################################\n\n" >&2

	echo -e "\n\n##################################################################################" >&2	     
	echo -e "#############  step one: alignment                                 ############" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x 	

	$bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index $R1 $R2 | $samtoolsdir/samtools view -@ $thr -bSu -> $alignedbam 
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	    exit $exitcode;
	fi   
	set +x
	echo -e "\n\n##################################################################################" >&2	     
	echo -e "#############  step two:making sure that a file was produced with alignments     #####" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x

	if [ -s $AlignDir/$alignedbam ]
	then 
	    set +x		           
	    echo -e "### the file was created. But we are not done.     #############" >&2
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
	    set -x

	    numAlignments=$( $samtoolsdir/samtools view -c $AlignDir/$alignedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for $AlignDir/$alignedbam\n" >> $qcfile	    
		MSG="bwa mem command did not produce alignments for $AlignDir/$alignedbam alignment failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"      
		exit 1;
	    else
		set +x
		echo -e "####### $AlignDir/$alignedbam seems to be in order ###########" >&2
		set -x 
	    fi
	else 
	    MSG="bwa mem command did not produce a file $AlignDir/$alignedbam alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi      
	set +x 
	echo -e "\n\n##################################################################################" >&2
	echo -e "#############  step three: sort                                        ############" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x

	$novocraftdir/novosort -t $tmpdir -c ${thr} -i -o $alignedsortedbam $alignedbam

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed during sorting for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit $exitcode
	fi 

	set +x
	echo -e "\n\n##################################################################################" >&2
	echo -e "#############  step four: dedup                                        ############" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x 


$javadir/java -Xmx8g -Djava.io.tmpdir=$tmpdir -jar $picardir/picard.jar  MarkDuplicates \
INPUT=$alignedsortedbam OUTPUT=$dedupsortedbam TMP_DIR=$tmpdir \
ASSUME_SORTED=true MAX_RECORDS_IN_RAM=null CREATE_INDEX=true \
METRICS_FILE=${SampleName}.picard.metrics \
VALIDATION_STRINGENCY=SILENT

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed during deduplication for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit $exitcode
	fi 

	set +x	
	echo -e "\n\n##################################################################################" >&2	     
	echo -e "#############  step five: making sure that a file was produced with alignments #######" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x
	
	if [ -s $AlignDir/$dedupsortedbam ]
	then
	    set +x			
	    echo -e "### the file was created. But we are not done.     #############" >&2
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
	    set -x 
	    
	    numAlignments=$( $samtoolsdir/samtools view -c $AlignDir/$dedupsortedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tpicard command did not produce a file for $AlignDir/$dedupsortedbam\n" >> $qcfile	    
		MSG="novosort command did not produce a file for $AlignDir/$dedupsortedbam"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		exit 1;
	    else
		set +x
		echo -e "####### $AlignDir/$dedupsortedbam seems to be in order ###########" >&2
		set -x
	    fi
	else 
	    MSG="picard command did not produce a file $AlignDir/$dedupsortedbam"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi   


        set +x 
	echo -e "\n\n##################################################################################" >&2
	echo -e "#############      END PICARD  BLOCK                                 ############" >&2
	echo -e "##################################################################################\n\n" >&2
	set -x             

	
else
	MSG="unrecognized deduplication tool $markduplicates"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;        

fi


set +x 
echo -e "\n\n##################################################################################" >&2
echo -e "#############     END ALIGNMENT-DEDUPLICATION BLOCK                   ############" >&2
echo -e "##################################################################################\n\n" >&2

echo `date`

echo -e "\n\n##################################################################################" >&2
echo -e "##################################################################################" >&2          
echo -e "##################################################################################" >&2        
echo -e "########   ALIGNMENT QC TEST   FOR SAMPLE $SampleName              ###############" >&2
echo -e "########   QC rule1: duplication cutoff <= $dup_cutoff             ###############" >&2
echo -e "########   QC rule2: mapped_reads cutoff >= $map_cutoff            ###############" >&2
echo -e "##################################################################################" >&2
echo -e "##################################################################################" >&2 
echo -e "##################################################################################\n\n" >&2
        
     

echo -e "\n\n##################################################################################" >&2	     
echo -e "#############  step one: generating the relevant file with flagstat       ############" >&2
echo -e "##################################################################################\n\n" >&2
set -x

flagstats=${dedupsortedbam}.flagstats

echo `date`             
$samtoolsdir/samtools flagstat $dedupsortedbam > $flagstats
echo `date`

set +x
echo -e "\n\n##################################################################################" >&2	     
echo -e "#############  step two: sanity check                                 ############" >&2
echo -e "##################################################################################\n\n" >&2
set -x  

if [ ! -s $flagstats ]
then
	 echo -e "${SampleName}\tQCT/EST\tFAIL\tsamtools/samtools flagstat command produced an empty file $flagstats\n" >> $qcfile
	 MSG="samtools flagstat command produced an empty file  $flagstats"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	 exit $exitcode;
fi

set +x
echo -e "\n\n##################################################################################" >&2	     
echo -e "#############  step three: parsing the file and grabbing stats for the QC test     ###" >&2
echo -e "##################################################################################\n\n" >&2
set -x           


tot_mapped=$( cat $flagstats | grep "mapped (" | cut -d ' ' -f1 )
tot_reads=$( cat $flagstats | grep "in total" | cut -d ' ' -f1 )
tot_dups=$( cat $flagstats | grep "duplicates" | cut -d ' ' -f1 )

#now testing if these variables are numeric and have numbers

if [ $tot_dups -eq $tot_dups 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi
if [ $tot_reads -eq $tot_reads 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi

if [ $tot_mapped -eq $tot_mapped 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi

set +x           
echo -e "\n\n##################################################################################" >&2	     
echo -e "#############  step four: calculating stats according to QC rules                  ###" >&2
echo -e "##################################################################################\n\n" >&2
set -x        


perc_dup=$(( tot_dups * 100 / tot_reads ))
perc_mapped=$(( tot_mapped * 100 / tot_reads ))

#now testing if these variables have numbers

if [ $perc_dup -eq $perc_dup 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi

if [ $perc_mapped -eq $perc_mapped 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi
set +x
echo -e "\n\n##################################################################################" >&2	     
echo -e "#############  step five: applying the  QC rules                                  ###" >&2
echo -e "##################################################################################\n\n" >&2
set -x         

if [ $perc_dup -lt $dup_cutoff ]
then
	echo -e "$#####  sample passed first filter percent_duplicates with value $perc_dup, maximum cutoff is $dup_cutoff"
	
	if [ $perc_mapped -gt $map_cutoff ]
	then
	        echo -e "##### $sample passed second filter map_cutoff with value $perc_mapped, minimum cutoff is $map_cutoff"	
	        echo -e "${SampleName}\tQCTEST\tPASS\trule1 ok: percent_duplication=$perc_dup <? duplication_cutoff=$dup_cutoff\trule2 ok: percent_mapped=$perc_mapped >? mapping_cutoff=$map_cutoff" >> $qcfile
	else
	        echo -e "##### $sample DID NOT pass second filter map_cutoff with value $perc_mapped, minimum cutoff is $map_cutoff"	
	        echo -e "${SampleName}\tQCTEST\tFAIL\trule1 ok: percent_duplication=$perc_dup <? duplication_cutoff=$dup_cutoff\trule2 notok: percent_mapped=$perc_mapped >? mapping_cutoff=$map_cutoff" >> $qcfile	          
	fi
else
	echo -e "$#####  sample DID NOT pass first filter percent_duplicates with value $perc_dup, maximum cutoff is $dup_cutoff"
	echo -e "${SampleName}\tQCTEST\tFAIL\trule1 not ok: percent_duplication=$perc_dup <? duplication_cutoff=$dup_cutoff\trule2 not evaluated: percent_mapped=$perc_mapped >? mapping_cutoff=$map_cutoff" >> $qcfile
fi

set +x
echo -e "\n\n##################################################################################" >&2
echo -e "#############       END QC TEST                                       ############" >&2        
echo -e "##################################################################################\n\n" >&2

echo `date`

echo -e "\n\n##################################################################################"  >&2
echo -e "##################################################################################" >&2
echo -e "##################################################################################" >&2        
echo -e "#############   WRAP UP                                               ############" >&2       
echo -e "##################################################################################" >&2
echo -e "##################################################################################" >&2 
echo -e "##################################################################################\n\n" >&2	

echo `date`
set -x

### perhaps this bam file is not necessary in the delivery folder           
### cp $AlignDir/${SampleName}.wdups.sorted.bam          $DeliveryDir   


MSG="ALIGNMENT-DEDUPLICATION for $SampleName finished successfully"
echo -e "program=$scriptfile at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"

echo `date`

set +x
echo -e "\n\n##################################################################################" >&2
echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############" >&2
echo -e "##################################################################################\n\n" >&2
set -x
