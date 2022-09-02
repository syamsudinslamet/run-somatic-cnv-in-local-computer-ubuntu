version 1.0

workflow Sample {
    input {
        File preprocessedIntervals
        File PON
        File? annotatedIntervals
        
        File? matchedNormalAllelicCounts
        
        File inputBam
        File inputBamIndex
        
        File commonVariantSites
        File? commonVariantSitesIndex
        
        String sampleName
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        
        String outputDir 
        
        Int? minimumContigLength

        String gatk
    }
    
    
    #>>> Task 1
    call CollectAllelicCounts {
    	input:
		allelicCountsPath = outputDir + "/" + sampleName + ".allelic_counts.tsv",
		commonVariantSites = commonVariantSites,
		commonVariantSitesIndex = commonVariantSitesIndex,
		inputBam = inputBam,
		inputBamIndex = inputBamIndex,
		referenceFasta = referenceFasta,
		referenceFastaDict = referenceFastaDict,
            	referenceFastaFai = referenceFastaFai,
		gatk = gatk
    }
    
    
    #>>> Task 2
    call CollectReadCounts {
        input:
            countsPath = outputDir + "/" + sampleName + ".readcounts.hdf5",
            intervals = preprocessedIntervals,
            inputBam = inputBam,
            inputBamIndex = inputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            gatk = gatk
    }
    

    #>>> Task 3
    call DenoiseReadCounts {
        input:
            PON = PON,
            annotatedIntervals = annotatedIntervals,
            readCounts = CollectReadCounts.counts,
            outputPrefix = outputDir + "/" + sampleName,
            gatk = gatk
    }
    
    #>>> Task 4
    call ModelSegments {
        input:
            outputDir = outputDir,
            outputPrefix = sampleName,
            denoisedCopyRatios = DenoiseReadCounts.denoisedCopyRatios,
            allelicCounts = CollectAllelicCounts.allelicCounts,
            normalAllelicCounts = matchedNormalAllelicCounts,
            gatk = gatk
    }

    #>>> Task 5
    call CallCopyRatioSegments {
        input:
            outputPrefix = outputDir + "/" + sampleName,
            copyRatioSegments = ModelSegments.copyRatioSegments,
            gatk = gatk
    }
    

    output {
        File allelicCounts = CollectAllelicCounts.allelicCounts
        File readCounts = CollectReadCounts.counts
        File standardizedCopyRatios = DenoiseReadCounts.standardizedCopyRatios
        File denoisedCopyRatios = DenoiseReadCounts.denoisedCopyRatios
        File hetrozygousAllelicCounts = ModelSegments.hetrozygousAllelicCounts
        File? normalHetrozygousAllelicCounts = ModelSegments.normalHetrozygousAllelicCounts
        File copyRatioSegments = ModelSegments.copyRatioSegments
        File copyRatioCBS = ModelSegments.copyRatioCBS
        File alleleFractionCBS = ModelSegments.alleleFractionCBS
        File unsmoothedModeledSegments = ModelSegments.unsmoothedModeledSegments
        File unsmoothedCopyRatioParameters = ModelSegments.unsmoothedCopyRatioParameters
        File unsmoothedAlleleFractionParameters = ModelSegments.unsmoothedAlleleFractionParameters
        File modeledSegments = ModelSegments.modeledSegments
        File copyRatioParameters = ModelSegments.copyRatioParameters
        File alleleFractionParameters = ModelSegments.alleleFractionParameters
        File calledSegments = CallCopyRatioSegments.calledSegments
        File calledSegmentsIgv = CallCopyRatioSegments.calledSegmentsIgv
    }

}

#>>> Define task here

#>>> Task 1
task CollectAllelicCounts {
    input {
    
        String allelicCountsPath
        File commonVariantSites
        File? commonVariantSitesIndex
        File inputBam
        File inputBamIndex
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
	
	String gatk
    }

    command {
        set -e
        mkdir -p "$(dirname ~{allelicCountsPath})"
        ~{gatk} \
        CollectAllelicCounts \
        -L ~{commonVariantSites} \
        -I ~{inputBam} \
        --read-index ~{inputBamIndex} \
        -R ~{referenceFasta} \
        -O ~{allelicCountsPath}
    }

    output {
        File allelicCounts = allelicCountsPath
    }

}

#>>> Task 2
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

#>>> Task 3
task DenoiseReadCounts {
    input {
        File readCounts
        String outputPrefix
        
        String gatk

        File PON
        File? annotatedIntervals

    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        ~{gatk} \
        DenoiseReadCounts \
        -I ~{readCounts} \
        --count-panel-of-normals ~{PON} \
        ~{"--annotated-intervals " + annotatedIntervals} \
        --standardized-copy-ratios ~{outputPrefix}.standardizedCR.tsv \
        --denoised-copy-ratios ~{outputPrefix}.denoisedCR.tsv
    }

    output {
        File standardizedCopyRatios = outputPrefix + ".standardizedCR.tsv"
        File denoisedCopyRatios = outputPrefix + ".denoisedCR.tsv"
    }

}


#>>> Task 4
task ModelSegments {
    input {
        String outputDir
        String outputPrefix
        File denoisedCopyRatios
        File allelicCounts
        Int minimumTotalAlleleCountCase = if defined(normalAllelicCounts) then 0 else 30
        Int maximumNumberOfSmoothingIterations = 10
        String gatk

        File? normalAllelicCounts

    }

    command {
        set -e
        mkdir -p ~{outputDir}
        ~{gatk} \
        ModelSegments \
        --denoised-copy-ratios ~{denoisedCopyRatios} \
        --allelic-counts ~{allelicCounts} \
        ~{"--normal-allelic-counts " + normalAllelicCounts} \
        --minimum-total-allele-count-case ~{minimumTotalAlleleCountCase} \
        --maximum-number-of-smoothing-iterations ~{maximumNumberOfSmoothingIterations} \
        --output ~{outputDir} \
        --output-prefix ~{outputPrefix}
    }

    output {
        File hetrozygousAllelicCounts = outputDir + "/" + outputPrefix + ".hets.tsv"
        File copyRatioSegments = outputDir + "/" + outputPrefix + ".cr.seg"
        File copyRatioCBS = outputDir + "/" + outputPrefix + ".cr.igv.seg"
        File alleleFractionCBS = outputDir + "/" + outputPrefix + ".af.igv.seg"
        File unsmoothedModeledSegments = outputDir + "/" + outputPrefix + ".modelBegin.seg"
        File unsmoothedCopyRatioParameters = outputDir + "/" + outputPrefix + ".modelBegin.cr.param"
        File unsmoothedAlleleFractionParameters = outputDir + "/" + outputPrefix + ".modelBegin.af.param"
        File modeledSegments = outputDir + "/" + outputPrefix + ".modelFinal.seg"
        File copyRatioParameters = outputDir + "/" + outputPrefix + ".modelFinal.cr.param"
        File alleleFractionParameters = outputDir + "/" + outputPrefix + ".modelFinal.af.param"
        File? normalHetrozygousAllelicCounts = outputDir + "/" + outputPrefix + ".hets.normal.tsv"
    }

}


#>>> Task 5
task CallCopyRatioSegments {
    input {
        String outputPrefix
        File copyRatioSegments

        String gatk
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        ~{gatk} \
        CallCopyRatioSegments \
        -I ~{copyRatioSegments} \
        -O ~{outputPrefix}.called.seg
    }

    output {
        File calledSegments = outputPrefix + ".called.seg"
        File calledSegmentsIgv = outputPrefix + ".called.igv.seg"
    }

}

