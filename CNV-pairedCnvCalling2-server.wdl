version 1.0

import "CNV-Task2-server.wdl" as sampleWorkflow

workflow PairedCnvCalling {
    input {
        File caseBam
        File caseBamIndex
        File controlBam
        File controlBamIndex
        File PON
        File? annotatedIntervals
        File preprocessedIntervals
        File commonVariantSites
        File? commonVariantSitesIndex
        String outputDir 
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        Int? minimumContigLength

        String gatk
    }
    
    String caseSampleName = basename(caseBam, ".bam") + "_T"
    String controlSampleName = basename(controlBam, ".bam") + "_N"

    call sampleWorkflow.Sample as controlSample {
        input:
            preprocessedIntervals = preprocessedIntervals,
            PON = PON,
            annotatedIntervals = annotatedIntervals,
            inputBam = controlBam,
            inputBamIndex = controlBamIndex,
            commonVariantSites = commonVariantSites,
            commonVariantSitesIndex = commonVariantSitesIndex,
            sampleName = controlSampleName,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            minimumContigLength = minimumContigLength,
            outputDir = outputDir + "/" + controlSampleName + "/",
            gatk = gatk
    }

    call sampleWorkflow.Sample as caseSample {
        input:
            preprocessedIntervals = preprocessedIntervals,
            PON = PON,
            annotatedIntervals = annotatedIntervals,
            matchedNormalAllelicCounts = controlSample.allelicCounts,
            inputBam = caseBam,
            inputBamIndex = caseBamIndex,
            commonVariantSites = commonVariantSites,
            commonVariantSitesIndex = commonVariantSitesIndex,
            sampleName = caseSampleName,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            minimumContigLength = minimumContigLength,
            outputDir = outputDir + "/" + caseSampleName + "/",
            gatk = gatk
    }

    output {
        File caseAllelicCounts = caseSample.allelicCounts
        File caseReadCounts = caseSample.readCounts
        File caseStandardizedCopyRatios = caseSample.standardizedCopyRatios
        File caseDenoisedCopyRatios = caseSample.denoisedCopyRatios
        File caseHetrozygousAllelicCounts = caseSample.hetrozygousAllelicCounts
        File caseNormalHetrozygousAllelicCounts = select_first([caseSample.normalHetrozygousAllelicCounts])
        File caseCopyRatioSegments = caseSample.copyRatioSegments
        File caseCopyRatioCBS = caseSample.copyRatioCBS
        File caseAlleleFractionCBS = caseSample.alleleFractionCBS
        File caseUnsmoothedModeledSegments = caseSample.unsmoothedModeledSegments
        File caseUnsmoothedCopyRatioParameters = caseSample.unsmoothedCopyRatioParameters
        File caseUnsmoothedAlleleFractionParameters = caseSample.unsmoothedAlleleFractionParameters
        File caseModeledSegments = caseSample.modeledSegments
        File caseCopyRatioParameters = caseSample.copyRatioParameters
        File caseAlleleFractionParameters = caseSample.alleleFractionParameters
        File caseCalledSegments = caseSample.calledSegments
        File caseCalledSegmentsIgv = caseSample.calledSegmentsIgv
        

        File controlAllelicCounts = controlSample.allelicCounts
        File controlReadCounts = controlSample.readCounts
        File controlStandardizedCopyRatios = controlSample.standardizedCopyRatios
        File controlDenoisedCopyRatios = controlSample.denoisedCopyRatios
        File controlHetrozygousAllelicCounts = controlSample.hetrozygousAllelicCounts
        File controlCopyRatioSegments = controlSample.copyRatioSegments
        File controlCopyRatioCBS = controlSample.copyRatioCBS
        File controlAlleleFractionCBS = controlSample.alleleFractionCBS
        File controlUnsmoothedModeledSegments = controlSample.unsmoothedModeledSegments
        File controlUnsmoothedCopyRatioParameters = controlSample.unsmoothedCopyRatioParameters
        File controlUnsmoothedAlleleFractionParameters = controlSample.unsmoothedAlleleFractionParameters
        File controlModeledSegments = controlSample.modeledSegments
        File controlCopyRatioParameters = controlSample.copyRatioParameters
        File controlAlleleFractionParameters = controlSample.alleleFractionParameters
        File controlCalledSegments = controlSample.calledSegments
        File controlCalledSegmentsIgv = controlSample.calledSegmentsIgv
    }

}
