@dispatch=WORKER
/*
void is set as an output so that dependency chaining can take place
(i.e., other functions can wait for this one to complete if necessary)
*/
app (void signal) samtools_index(string samtoolsdir, file inputFilename) {
	samtoolsdir "index" inputFilename;
}

@dispatch=WORKER
app (file outputfile, file outLog) RealignerTargetCreator (string javaexe, string java_heap, string gatkjar, string reference,
   file inputFile, int thr, string known[]
  ) {
	javaexe java_heap "-jar" gatkjar
	"-T" "RealignerTargetCreator"
	"-R" reference
	"-I" inputFile known
	"-nt" thr
	"-o" outputfile
	@stderr=outLog;
}

@dispatch=WORKER
app (file outputfile, file outLog) IndelRealigner (string javaexe, string java_heap, string gatkjar, string reference,
				   file inputFile, string known[], file intervals
				  ) {
	javaexe java_heap "-jar" gatkjar
	"-T" "IndelRealigner"
	"-R" reference
	"-I" inputFile known 
	"--targetIntervals" intervals
	"-o" outputfile
	@stderr=outLog;
}

@dispatch=WORKER
app (file outputfile, file outLog) BaseRecalibrator (string javaexe, string java_heap, string gatkjar, string reference,
				     file inputFile, int thr, string knownindels[], string dbsnp
				    ) {
	javaexe java_heap "-jar" gatkjar
	"-T" "BaseRecalibrator"
	"-R" reference
	"-I" inputFile
	"-knownSites" dbsnp knownindels
	"-nct" thr
	"-o" outputfile
	@stderr=outLog;
}

@dispatch=WORKER
app (file outputfile, file outLog) PrintReads (string javaexe, string java_heap, string gatkjar, string reference,
			       file inputFile, int thr, file recalreport
			      ) {
	javaexe java_heap "-jar" gatkjar
	"-T" "PrintReads"
	"-R" reference
	"-I" inputFile
	"-BQSR" recalreport
	"-nct" thr
	"--out" outputfile
	@stderr=outLog;
}


// The chromosome splitting version
@dispatch=WORKER
app (file outputfile, file outLog) HaplotypeCaller(string javaexe, string java_heap, string gatkjar, string reference,
						   file inputFile, string dbsnp, int thr, int ploidy, string chr
						  ) {
	javaexe java_heap "-jar" gatkjar
	"-T" "HaplotypeCaller"
	"-R" reference
	"-I" inputFile
	"--dbsnp" dbsnp
	"-o" outputfile
	"--emitRefConfidence" "GVCF"
	"-gt_mode" "DISCOVERY"
	"-A" "Coverage"
	"-A" "FisherStrand"
	"-A" "StrandOddsRatio"
	"-A" "HaplotypeScore"
	"-A" "MappingQualityRankSumTest"
	"-A" "QualByDepth"
	"-A" "RMSMappingQuality"
	"-A" "ReadPosRankSumTest"
	"--sample_ploidy" ploidy
	"-nt" 1
	"-nct" thr
	"-L" chr
	@stderr=outLog;
}

// The whole genome version
@dispatch=WORKER	   
app (file outputfile, file outLog) HaplotypeCaller(string javaexe, string java_heap, string gatkjar, string reference,      
						   file inputFile, string dbsnp, int thr	  
						  ) {			     
	javaexe java_heap "-jar" gatkjar			    
	"-T" "HaplotypeCaller"				     
	"-R" reference     
	"-I" inputFile     
	"--dbsnp" dbsnp    
	"-o" outputfile    
	"--emitRefConfidence" "GVCF"			       
	"-gt_mode" "DISCOVERY"				     
	"-A" "Coverage"    
	"-A" "FisherStrand"
	"-A" "StrandOddsRatio"				     
	"-A" "HaplotypeScore"				      
	"-A" "MappingQualityRankSumTest"			   
	"-A" "QualByDepth" 
	"-A" "RMSMappingQuality"				   
	"-A" "ReadPosRankSumTest"				  			      
	"-nt" 1	    
	"-nct" thr				    
	@stderr=outLog;    
}
