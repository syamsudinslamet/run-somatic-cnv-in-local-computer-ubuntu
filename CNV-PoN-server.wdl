version 1.0

workflow PanelofNormal {
	#>>> 1. set global input
	input{
		Array[File] inputBams = [] 
		Array[File] inputBamIndexes = []
		File referenceFasta
		File referenceFastaDict
		File referenceFastaFai
		File? regions
		String outputDir
		String PONname = "cnv.pon"
		Boolean performExplicitGcCorrection = true
		String gatk
		

	}
	
	#>>> 2. Call preprocess interval for an input in collect read count
	call PreprocessIntervals {
		input:
		    referenceFasta = referenceFasta,
		    referenceFastaDict = referenceFastaDict,
		    referenceFastaFai = referenceFastaFai,
		    intervals = regions,
		    outputIntervalList = outputDir + "/preprocessed.interval_list",
		    gatk = gatk
    	}
	
	#>>> 3. Additional step: GC Correction, not necessary for matched normal-case
	if (performExplicitGcCorrection) {
		call AnnotateIntervals {
		    input:
		        referenceFasta = referenceFasta,
		        referenceFastaDict = referenceFastaDict,
		        referenceFastaFai = referenceFastaFai,
		        annotatedIntervalsPath = outputDir + "/annotated_intervals.tsv",
		        intervals = PreprocessIntervals.intervalList,
		        gatk = gatk
		}
	}
	
	#>>> 4. Collect coverage count for normal - do it parallelly
	scatter (bamAndIndex in zip(inputBams, inputBamIndexes)) {
		call CollectReadCounts {
		    input:
		        countsPath = outputDir + "/" + basename(bamAndIndex.left) + ".readcounts.hdf5",
		        intervals =  PreprocessIntervals.intervalList,
		        inputBam = bamAndIndex.left,
		        inputBamIndex = bamAndIndex.right,
		        referenceFasta = referenceFasta,
		        referenceFastaDict = referenceFastaDict,
		        referenceFastaFai = referenceFastaFai,
		        gatk = gatk
		}
	    }
	    
	
	#>>> 5. Create Panel of Normal
	if (length(CollectReadCounts.counts) > 0){
		call CreateReadCountPanelOfNormals {
		    input:
		        PONpath = outputDir + "/" + PONname + ".hdf5",
		        readCountsFiles = CollectReadCounts.counts,
		        annotatedIntervals = AnnotateIntervals.annotatedIntervals,
		        gatk = gatk
        	}
    	}

	output {
		File preprocessedIntervals = PreprocessIntervals.intervalList
		File? PON = CreateReadCountPanelOfNormals.PON
		File? annotatedIntervals = AnnotateIntervals.annotatedIntervals
	}
	
}	

#Declare task here:
#>>> 2. preprocess interval task	
task PreprocessIntervals {
    input {
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String outputIntervalList 
        Int binLength = if defined(intervals) then 0 else 1000
        Int padding = if defined(intervals) then 250 else 0
        String intervalMergingRule = "OVERLAPPING_ONLY"

        File? intervals
	String gatk
    }

    command {
        set -e #terminates the execution when the error occurs
        mkdir -p "$(dirname ~{outputIntervalList})"
        ~{gatk} \
        PreprocessIntervals \
        -R ~{referenceFasta} \
        --sequence-dictionary ~{referenceFastaDict} \
        --bin-length ~{binLength} \
        --padding ~{padding} \
        ~{"-L " + intervals} \
        --interval-merging-rule ~{intervalMergingRule} \
        -O ~{outputIntervalList}
    }

    output {
        File intervalList = outputIntervalList
    }

}


#>>> 3. Annotate Interval
task AnnotateIntervals {
    input {
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String annotatedIntervalsPath
        File intervals
        String intervalMergingRule = "OVERLAPPING_ONLY"
        Int featureQueryLookahead = 1000000

        File? mappabilityTrack
        File? segmentalDuplicationTrack

        String gatk
    }

    command {
        set -e
        mkdir -p "$(dirname ~{annotatedIntervalsPath})"
        ~{gatk}  \
        AnnotateIntervals \
        -R ~{referenceFasta} \
        -L ~{intervals} \
        ~{"--mappability-track  " + mappabilityTrack} \
        ~{"--segmental-duplication-track " + segmentalDuplicationTrack} \
        --feature-query-lookahead ~{featureQueryLookahead} \
        --interval-merging-rule ~{intervalMergingRule} \
        -O ~{annotatedIntervalsPath}
    }

    output {
        File annotatedIntervals = annotatedIntervalsPath
    }

}


#>>> 4. Collect coverage count
task CollectReadCounts {
    input {
        String countsPath
        File intervals
        File inputBam
        File inputBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String intervalMergingRule = "OVERLAPPING_ONLY"
      
        String gatk
    }

    command {
        set -e
        mkdir -p "$(dirname ~{countsPath})"
        ~{gatk} \
        CollectReadCounts \
        -L ~{intervals} \
        -I ~{inputBam} \
        -R ~{referenceFasta} \
        --format HDF5 \
        --interval-merging-rule ~{intervalMergingRule} \
        -O ~{countsPath}
    }

    output {
        File counts = countsPath
    }

}

#>>> 5. Create PoN
task CreateReadCountPanelOfNormals {
    input {
        String PONpath
        Array[File]+ readCountsFiles

        File? annotatedIntervals #adding this does not work with 1 sample

        String gatk
    }

    command {
        set -e
        mkdir -p "$(dirname ~{PONpath})"
        ~{gatk} \
        CreateReadCountPanelOfNormals \
        -I ~{sep=" -I " readCountsFiles} \
        --minimum-interval-median-percentile 5.0 \
        -O ~{PONpath}
    }

    output {
        File PON = PONpath
    }

  
}

































