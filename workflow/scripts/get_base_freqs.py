#################################################################################################################################################################
###     This program extracts position specific nucleotide frequencies from all vcf files in the given file hierarchy.
###     *  Can exclude deletion read counts by setting "noDel=True".
###     *  Columns in vcf file: CHROM    POS	ID	REF	ALT	QUAL	FILTER	INFO
###     *  Example line:
###        NC_045512.2	241	.	T	C	4.77121	PASS	Freq1=1;Freq2=1;Freq3=0;Post1=1.324;Post2=1.0286;Post3=0;Fvar=2793;Rvar=2818;Ftot=2794;Rtot=2819;Pval=1;Qval=1
#################################################################################################################################################################

import sys
import os
import glob


collectionName = snakemake.input['vcf_dir']
delCoverageThreshold = 0     # by default deletions with any coverage are used
delThreshold = 0.0           # by default deletions at all frequencies are used
delThreshold = snakemake.config['params']['deletion_threshold']
delCoverageThreshold = snakemake.config['params']['deletion_coverage_threshold']
#list = sys.argv[4]
#usedSamples = []


#f = open(list,'r')
#lines = f.readlines()[1:]
#f.close()

#for line in lines:
#    usedSamples.append(line.split(',')[0])
#print len(usedSamples)


delFilter = ""
if delThreshold > 0.0 or delCoverageThreshold > 0:
    delFilter = "_delFilter_" + str(delThreshold) + "_" + str(delCoverageThreshold)


outFileName = snakemake.output['fname']
outFile = open(outFileName, "w")
outFile.write("SAMPLE" + "\t" + "POS" + "\t" + "REF_BASE" + "\t" + "(ADJUSTED_)READ_COUNT" + "\t" + "A_freq" + "\t" + "C_freq" + "\t" + "G_freq" + "\t" + "T_freq" + "\t" + "DEL_freq"+ "\n")
print(outFileName)
#pattern = "samples/*/*/variants/SNVs/snvs.vcf"
pattern = collectionName + "/*.vcf"       # all files matching this pattern are processed
fileList = glob.glob(pattern)

usedFiles = 0
for file in fileList:
    countMuts = 0
    sampleName = os.path.basename(file)
    if sampleName.endswith('.vcf'):
        sampleName  = sampleName[:-4]   # remove file extension to get sample name
    if sampleName.startswith('snvs_'):
        sampleName = sampleName[5:]
    #if sampleName not in usedSamples:
    #    print sampleName
    #    continue
    usedFiles += 1
    #sampleDate = file.split('/')[2]
    posBaseCounts = {}      # pos -> base counts
    refBaseAtPos = {}       # pos -> refBase
    postProbs = {}
    prevLine = ""

    with open(file) as f:
        #print file
        lines = f.read().splitlines()
        for line in lines:

            if line.startswith('#'):   #skip lines not referring to a mutation
                continue

             ## get the base read counts from the line
            elems = line.split('\t')
            pos= int(elems[1])
            altBase = elems[4]
            refBase = elems[3]
            refBaseAtPos[pos] = refBase
            info = elems[7].split(';')
            Fvar = elems[7].split("Fvar=")[1].split(";")[0]
            Rvar = elems[7].split("Rvar=")[1].split(";")[0]
            Ftot = elems[7].split("Ftot=")[1].split(";")[0]
            Rtot = elems[7].split("Rtot=")[1].split(";")[0]
            total = int(Ftot) + int(Rtot)
            altTotal = int(Fvar) + int(Rvar)
            mutFreq = float(altTotal)/float(total)


            if pos not in posBaseCounts:
                baseCounts = dict([('A', 0), ('C', 0), ('G', 0), ('T', 0), ('-', 0), ('Total', 0)])
                posBaseCounts[pos] = baseCounts


            ## some checks that read counts in multiple rows refering to same position make sense
            if posBaseCounts[pos][altBase] != 0:
                if posBaseCounts[pos][altBase] != altTotal:
                    print("Warning: different alt base count: " + str(posBaseCounts[pos][altBase]) + " = " + str(altTotal))

            if posBaseCounts[pos]['Total'] != 0 and posBaseCounts[pos]['Total'] != total:
                print("Warning: different total base count: " + str(posBaseCounts[pos]['Total']) + " = " + str(total))

            ## set alt base and total read count for this position
            posBaseCounts[pos][altBase] = altTotal
            posBaseCounts[pos]['Total'] = total

            prevLine = line


    ## assign difference between total read count and sum of alt read counts as ref count
    for pos in posBaseCounts:
        refCount = posBaseCounts[pos]['Total'] - (posBaseCounts[pos]['A'] + posBaseCounts[pos]['C'] + posBaseCounts[pos]['G'] + posBaseCounts[pos]['T'] + posBaseCounts[pos]['-'])
        posBaseCounts[pos][refBaseAtPos[int(pos)]] = refCount

        ## another check that read counts make sense
        if refCount < 0:
            print("ERROR: NEGATIVE REF COUNT AT POSITION " + pos + " in " + file + ":")
            print(posBaseCounts[pos])


    ## write base frequencies to output file
    for pos, refBase in sorted(refBaseAtPos.items()):

        ## check if deletion makes the threshold at this postion
        useDel = True
        if float(posBaseCounts[pos]['-'])/float(posBaseCounts[pos]['Total']) <= delThreshold:
            useDel = False
        if posBaseCounts[pos]['Total'] <= delCoverageThreshold:
            useDel = False

        ## if we don't consider the deletion, update total considered read counts
        totalConsidered = posBaseCounts[pos]['Total']
        if not useDel:
            totalConsidered = posBaseCounts[pos]['Total']-posBaseCounts[pos]['-']

        ## check if we should skip this position
        if totalConsidered == 0:                            # no reads left after removing deletion
            continue
        if posBaseCounts[pos][refBase] == totalConsidered:  # no mutation left after removing deletions
            continue

        # compute variant frequencies from base counts
        baseFreqStats = str(float(posBaseCounts[pos]['A'])/float(totalConsidered)) + "\t"
        baseFreqStats += str(float(posBaseCounts[pos]['C'])/float(totalConsidered)) + "\t"
        baseFreqStats += str(float(posBaseCounts[pos]['G'])/float(totalConsidered)) + "\t"
        baseFreqStats += str(float(posBaseCounts[pos]['T'])/float(totalConsidered))
        if useDel:
            baseFreqStats += "\t" + str(float(posBaseCounts[pos]['-'])/float(totalConsidered))
        else:
            baseFreqStats += "\t" + "0"

        posInfo = sampleName + "\t" + str(pos) + "\t" + refBaseAtPos[pos]
        outFile.write( posInfo + "\t" + str(totalConsidered) + "\t" + baseFreqStats + "\n")
        countMuts += 1
    #print str(countMuts) + "\t" + sampleName
print(usedFiles)

outFile.close()
