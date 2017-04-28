#!/usr/bin/python



##add your name below if you are editing this
'''
pipeline contributors:
Eric Karlins


This is to run a SnakeMake pipeline doing GWAS QC and generating several QC metrics files to return to give to the lab
'''

import glob
import sys
import os
import shutil
import codecs
from snakemake.utils import R


configfile: "config.yaml"

input_sample_sheet = config['sample_sheet']
(outName, x, sampSheetDate) = os.path.basename(input_sample_sheet)[:-4].split('_')
baseProjDir = '/DCEG/CGF/Infinium/ScanData/CGF/ByProject'
snp_cr_1 = config['snp_cr_1']
samp_cr_1 = config['samp_cr_1']
snp_cr_2 = config['snp_cr_2']
samp_cr_2 = config['samp_cr_2']
ld_prune_r2 = config['ld_prune_r2']
maf_for_ibd = config['maf_for_ibd']
sample_sheet = 'IlluminaSampleSheet.csv'
subject_id_to_use = config['subject_id_to_use']
ibd_pi_hat_cutoff = config['ibd_pi_hat_cutoff']
dup_concordance_cutoff = config['dup_concordance_cutoff']
illumina_manifest_file = config['illumina_manifest_file']
expected_sex_col_name = config['expected_sex_col_name']
lims_output_dir = config['lims_output_dir']
contam_threshold = config['contam_threshold']


def makeRepDiscordantDict(knownConcordanceFile, thresh = dup_concordance_cutoff):
    RepDiscordantDict = {}
    with open(knownConcordanceFile) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            (Subject_ID,Sample_ID1,Sample_ID2,Concordance,PI_HAT) = line.rstrip().split(',')
            if Concordance != 'NA' and float(Concordance) < thresh:
                RepDiscordantDict[Sample_ID1] = 1
                RepDiscordantDict[Sample_ID2] = 1
            line = f.readline()
    return RepDiscordantDict


def makeUnexpectedRepDict(unknownRepFile):
    UnexpectedRepDict = {}
    with open(unknownRepFile) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            (Subject_ID1,Subject_ID2,Sample_ID1,Sample_ID2,Concordance,PI_HAT) = line.rstrip().split(',')
            UnexpectedRepDict[Sample_ID1] = 1
            UnexpectedRepDict[Sample_ID2] = 1
            line = f.readline()
    return UnexpectedRepDict



def makeSubjectToSampListDict(SampleSheet):
    subToSampListDict = {}
    sampToSubIdDict = {}
    with codecs.open(SampleSheet,"r",encoding='utf-8', errors='ignore') as f:
        head = f.readline()
        while 'SentrixBarcode_A' not in head and head != '':
            head = f.readline()
        if 'SentrixBarcode_A' not in head:
            print('Sample sheet not formatted correctly')
            sys.exit(1)
        head_list = head.rstrip().split(',')
        subjectIdCol = None
        for i in range(len(head_list)):
            if head_list[i] == subject_id_to_use:
                subjectIdCol = i
        if subjectIdCol == None:
            print('Subject ID not found in sample sheet')
            sys.exit(1)
        line = f.readline()
        while line != '':
            if line.strip():
                line_list = line.rstrip().split(',')
                sampId = line_list[0]
                subId = line_list[subjectIdCol]
                sampToSubIdDict[sampId] = subId
                if not subToSampListDict.get(subId):
                    subToSampListDict[subId] = []
                subToSampListDict[subId].append(sampId)
            line = f.readline()
    return (subToSampListDict, sampToSubIdDict)


def makeSampToExpectedSexDict(SampleSheet):
    sampToExpectedSexDict = {}
    with codecs.open(SampleSheet,"r",encoding='utf-8', errors='ignore') as f:
        head = f.readline()
        while 'SentrixBarcode_A' not in head and head != '':
            head = f.readline()
        if 'SentrixBarcode_A' not in head:
            print('Sample sheet not formatted correctly')
            sys.exit(1)
        head_list = head.rstrip().split(',')
        expectedSexCol = None
        for i in range(len(head_list)):
            if head_list[i] == expected_sex_col_name:
                expectedSexCol = i
        if expectedSexCol == None:
            print('Expected Sex not found in sample sheet')
            sys.exit(1)
        line = f.readline()
        while line != '':
            if line.strip():
                line_list = line.rstrip().split(',')
                sampId = line_list[0]
                expectedSex = line_list[expectedSexCol][0].upper()
                sampToExpectedSexDict[sampId] = expectedSex
            line = f.readline()
    return sampToExpectedSexDict





def getConcIBS(ibs0, ibs1, ibs2):
    '''
    (int, int, int) -> float
    return the concordance when given the ibs values from the PLINK IBD file
    '''
    return float(ibs2)/(ibs0 + ibs1 + ibs2)




def makeCallRateDict(imissFile):
    '''
    (str) -> dict
    return a dict where the key is sample id and the value is call rate
    '''
    crDict = {}
    with open(imissFile) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            line_list = line.split()
            samp = line_list[1]
            miss = float(line_list[-1])
            cr = 1.0 - miss
            crDict[samp] = cr
            line = f.readline()
    return crDict



def makeSexDict(sexCheckFile):
    sexDict = {}
    with open(sexCheckFile) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            line_list = line.split()
            samp = line_list[1]
            snpsex = line_list[3]
            if snpsex == '1':
                snpsex = 'M'
            elif snpsex == '2':
                snpsex = 'F'
            else:
                snpsex = 'U'
            homozygosity = float(line_list[-1])
            het = 1.0 - homozygosity
            sexDict[samp] = (snpsex, het)
            line = f.readline()
    return sexDict



def makeContamDict(ContamFile):
    contamDict = {}
    with open(ContamFile) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            line_list = line.split(',')
            samp = line_list[0]
            contam = line_list[1]
            contamDict[samp] = contam
            line = f.readline()
    return contamDict



def makeAncestryDict(snpWeightsFile):
    ancestryDict = {}
    with open(snpWeightsFile) as f:
        for line in f:
            line_list = line.split()
            samp = line_list[0]
            (AFR, EUR, ASN) = line_list[-3:]
            ancestryDict[samp] = (float(AFR), float(EUR), float(ASN)) 
    return ancestryDict


def makeChipIdToSampDict(SampleSheet):
    allSampleIds = []
    idats = []
    SAMPLE_IDS = []
    gtcFiles = []
    chipIdToSampDict = {}
    sampIdToGtcDict = {}
    sampIdToProjDict = {}
    with codecs.open(SampleSheet,"r",encoding='utf-8', errors='ignore') as f:
        head = f.readline()
        while 'SentrixBarcode_A' not in head and head != '':
            head = f.readline()
        if 'SentrixBarcode_A' not in head:
            print('Sample sheet not formatted correctly')
            sys.exit(1)
        projCol = None
        head_list = head.rstrip().split(',')
        for i in range(len(head_list)):
            if head_list[i] == 'Project':
                projCol = i
        if not projCol:
            print('No column named "Project" found')
            sys.exit(1)
        line = f.readline()
        while line != '':
            if line.strip():
                line_list = line.rstrip().split(',')
                sampId = line_list[0]
                allSampleIds.append(sampId)
                chipId = line_list[1] + '_' + line_list[2]
                proj = line_list[projCol]
                sampIdToProjDict[sampId] = proj
                idatBase = baseProjDir + '/' + proj + '/' + line_list[1] + '/' + chipId
                if os.path.isfile(idatBase + '.gtc'):
                    gtcFiles.append(idatBase + '.gtc')
                    SAMPLE_IDS.append(sampId)
                    sampIdToGtcDict[sampId] = idatBase + '.gtc'
                if os.path.isfile(idatBase + '_Red.idat') and os.path.isfile(idatBase + '_Grn.idat'):
                    idats.append(idatBase)
                chipIdToSampDict[chipId] = sampId
            line = f.readline()
    return (chipIdToSampDict, sampIdToGtcDict, sampIdToProjDict, allSampleIds, idats, SAMPLE_IDS, gtcFiles, chipIdToSampDict)

(chipIdToSampDict, sampIdToGtcDict, allSampleIds, sampIdToProjDict, idats, SAMPLE_IDS, gtcFiles, chipIdToSampDict) = makeChipIdToSampDict(sample_sheet)


def getGtc(wildcards):
    sampId = wildcards.sample_id
    return sampIdToGtcDict[sampId]



if config['plink_genotype_file'] == 'None':
    include: 'modules/Snakefile_gtc_preprocess'

elif config['plink_genotype_file'][-4:] == '.ped':
    include: 'modules/Snakefile_ped_preprocess'

elif config['plink_genotype_file'][-4:] == '.bed':
    include: 'modules/Snakefile_bed_preprocess'


include: 'modules/Snakefile_plink_stats_filters'
include: 'modules/Snakefile_replicate_concordance'
include: 'modules/Snakefile_sample_qc_report'


rule all:
    input:
        'all_sample_qc.csv',
        'concordance/KnownReplicates.csv',
        'concordance/UnknownReplicates.csv',
        'snpweights/samples.snpweights',
        'all_contam/contam.csv'


rule plink_sample_qc_stats:
    input:
        'plink_start/samples.bed',
        'plink_start/samples.bim',
        'plink_start/samples.fam'
    params:
        inProj = 'plink_start/samples',
        outProj = 'sample_qc_stats/samples'
    output:
        'sample_qc_stats/samples.imiss',
        'sample_qc_stats/samples.lmiss',
        'sample_qc_stats/samples.sexcheck',
        'sample_qc_stats/samples.frq',
        'sample_qc_stats/samples.hwe'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --freq --missing --hardy --het --check-sex --out {params.outProj}'


rule filter_call_rate_1:
    input:
        'plink_start/samples.bed',
        'plink_start/samples.bim',
        'plink_start/samples.fam'
    params:
        inProj = 'plink_start/samples',
        outProj = 'plink_filter_call_rate_1/samples',
        mind = str(1 - float(samp_cr_1)),
        geno = str(1 - float(snp_cr_1))
    output:
        'plink_filter_call_rate_1/samples.bed',
        'plink_filter_call_rate_1/samples.bim',
        'plink_filter_call_rate_1/samples.fam'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --mind {params.mind} --geno {params.geno} --make-bed --out {params.outProj}'

rule qc_stats_after_filter_1:
    input:
        'plink_filter_call_rate_1/samples.bed',
        'plink_filter_call_rate_1/samples.bim',
        'plink_filter_call_rate_1/samples.fam'
    params:
        inProj = 'plink_filter_call_rate_1/samples',
        outProj = 'plink_filter_call_rate_1/samples_filter1'
    output:
        'plink_filter_call_rate_1/samples_filter1.imiss',
        'plink_filter_call_rate_1/samples_filter1.lmiss',
        'plink_filter_call_rate_1/samples_filter1.sexcheck',
        'plink_filter_call_rate_1/samples_filter1.frq',
        'plink_filter_call_rate_1/samples_filter1.hwe'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --freq --missing --hardy --het --check-sex --out {params.outProj}'




rule filter_call_rate_2:
    input:
        'plink_filter_call_rate_1/samples.bed',
        'plink_filter_call_rate_1/samples.bim',
        'plink_filter_call_rate_1/samples.fam'
    params:
        inProj = 'plink_filter_call_rate_1/samples',
        outProj = 'plink_filter_call_rate_2/samples',
        mind = str(1 - float(samp_cr_2)),
        geno = str(1 - float(snp_cr_2))
    output:
        'plink_filter_call_rate_2/samples.bed',
        'plink_filter_call_rate_2/samples.bim',
        'plink_filter_call_rate_2/samples.fam'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --mind {params.mind} --geno {params.geno} --make-bed --out {params.outProj}'



rule qc_stats_after_filter_2:
    input:
        'plink_filter_call_rate_2/samples.bed',
        'plink_filter_call_rate_2/samples.bim',
        'plink_filter_call_rate_2/samples.fam'
    params:
        inProj = 'plink_filter_call_rate_2/samples',
        outProj = 'plink_filter_call_rate_2/samples_filter2'
    output:
        'plink_filter_call_rate_2/samples_filter2.imiss',
        'plink_filter_call_rate_2/samples_filter2.lmiss',
        'plink_filter_call_rate_2/samples_filter2.sexcheck',
        'plink_filter_call_rate_2/samples_filter2.frq',
        'plink_filter_call_rate_2/samples_filter2.hwe'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --freq --missing --hardy --het --check-sex --out {params.outProj}'




rule plink_ld_prune:
    input:
        'plink_filter_call_rate_2/samples.bed',
        'plink_filter_call_rate_2/samples.bim',
        'plink_filter_call_rate_2/samples.fam'
    params:
        inProj = 'plink_filter_call_rate_2/samples',
        outProj = 'ld_prune/ldPruneList',
        r2 = ld_prune_r2,
        maf = maf_for_ibd
    output:
        'ld_prune/ldPruneList.prune.in',
        'ld_prune/ldPruneList.prune.out'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --indep-pairwise 50 5 {params.r2}  --maf {params.maf} --out {params.outProj}'


rule extract_ld_prune:
    input:
        bed = 'plink_filter_call_rate_2/samples.bed',
        bim = 'plink_filter_call_rate_2/samples.bim',
        fam = 'plink_filter_call_rate_2/samples.fam',
        prune = 'ld_prune/ldPruneList.prune.in'
    params:
        inProj = 'plink_filter_call_rate_2/samples',
        outProj = 'ld_prune/samples'
    output:
        'ld_prune/samples.bed',
        'ld_prune/samples.bim',
        'ld_prune/samples.fam'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --extract {input.prune} --make-bed --out {params.outProj}'


rule plink_ibd:
    input:
        'ld_prune/samples.bed',
        'ld_prune/samples.bim',
        'ld_prune/samples.fam'
    params:
        inProj = 'ld_prune/samples',
        outProj = 'ibd/samples'
    output:
        'ibd/samples.genome'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --genome full --min 0.05 --out {params.outProj}'

'''
# I only need this rule if I'm writing a script to calculate concordance.
# I'll calculate concordance from IBD output for now

rule plink_recode_A:
    input:
        'ld_prune/samples.bed',
        'ld_prune/samples.bim',
        'ld_prune/samples.fam'
    params:
        inProj = 'ld_prune/samples',
        outProj = 'snp_concordance/samples'
    output:
        'snp_concordance/samples.raw'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --recode A --out {params.outProj}'
'''



rule output_replicates:
    input:
        sampSheet = sample_sheet,
        imiss3 = 'plink_filter_call_rate_2/samples_filter2.imiss',
        ibd = 'ibd/samples.genome'
    output:
        known = 'concordance/KnownReplicates.csv',
        unknown = 'concordance/UnknownReplicates.csv'
    run:
        (SubToSampListDict, sampToSubIdDict) = makeSubjectToSampListDict(sample_sheet)
        crDict = makeCallRateDict(input.imiss3)
        piHatDict = {}
        minConcordance = 1.0
        with open(input.ibd) as f:
            head = f.readline()
            (piHatCol, ibs0col, ibs1col, ibs2col) = [None, None, None, None]
            head_list = head.split()
            for i in range(len(head_list)):
                if head_list[i] == 'PI_HAT':
                    piHatCol = i
                elif head_list[i] == 'IBS0':
                    ibs0col = i
                elif head_list[i] == 'IBS1':
                    ibs1col = i
                elif head_list[i] == 'IBS2':
                    ibs2col = i
            if not piHatCol or not ibs0col or not ibs1col or not ibs2col:
                print('Necessary column headers not in IBD output.')
                sys.exit(1)
            line = f.readline()
            while line != '':
                line_list = line.split()
                samp1 = line_list[1]
                samp2 = line_list[3]
                piHat = float(line_list[piHatCol])
                ibs0 = float(line_list[ibs0col])
                ibs1 = float(line_list[ibs1col])
                ibs2 = float(line_list[ibs2col])
                concordance = getConcIBS(ibs0, ibs1, ibs2)
                if concordance < minConcordance:
                    minConcordance = concordance
                sampList = sorted([samp1, samp2])
                piHatDict[(sampList[0], sampList[1])] = (piHat, concordance)
                line = f.readline()
        with open(output.known, 'w') as outKnown:
            outKnown.write('Subject_ID,Sample_ID1,Sample_ID2,Concordance,PI_HAT\n')
            for subId in SubToSampListDict:
                sampList = sorted(set(SubToSampListDict[subId]))
                if len(sampList) > 1:
                    for i in range(len(sampList)):
                        for j in range(i + 1, len(sampList)):
                            samp1 = sampList[i]
                            samp2 = sampList[j]
                            if not piHatDict.get((samp1, samp2)):
                                piHat = .05
                                concordance = minConcordance
                            else:
                                (piHat, concordance) = piHatDict[(samp1, samp2)]
                            if not crDict.get(samp1) or not crDict.get(samp2):
                                concordance = 'NA'
                                piHat = 'NA'
                            outKnown.write(','.join([subId, samp1, samp2, str(concordance), str(piHat)]) + '\n')
        with open(output.unknown, 'w') as outUn:
            outUn.write('Subject_ID1,Subject_ID2,Sample_ID1,Sample_ID2,Concordance,PI_HAT\n')
            for (samp1, samp2) in piHatDict.keys():
                (piHat, concordance) = piHatDict[(samp1,samp2)]
                if concordance > dup_concordance_cutoff:
                    subId1 = sampToSubIdDict[samp1]
                    subId2 = sampToSubIdDict[samp2]
                    if subId1 != subId2 and crDict.get(samp1) and crDict.get(samp2):
                        outUn.write(','.join([subId1, subId2, samp1, samp2, str(concordance), str(piHat)]) + '\n')





rule sample_qc_report:
    input:
        sampSheet = sample_sheet,
        imiss1 = 'sample_qc_stats/samples.imiss',
        imiss2 = 'plink_filter_call_rate_1/samples_filter1.imiss',
        imiss3 = 'plink_filter_call_rate_2/samples_filter2.imiss',
        sexcheck = 'sample_qc_stats/samples.sexcheck',
        freq = 'sample_qc_stats/samples.frq',
        hwe = 'sample_qc_stats/samples.hwe',
        contam = 'all_contam/contam.csv',
        ancestry = 'snpweights/samples.snpweights',
        knownConcordance = 'concordance/KnownReplicates.csv',
        unknownRep = 'concordance/UnknownReplicates.csv'
    output:
        allQc = 'all_sample_qc.csv',
        lims = outName + '_LimsUpload_' + sampSheetDate + '.csv',
        limsDirOut = lims_output_dir + '/' + outName + '_LimsUpload_' + sampSheetDate + '.csv'
    run:
        (SubToSampListDict, sampToSubIdDict) = makeSubjectToSampListDict(sample_sheet)
        ExpectedSexDict = makeSampToExpectedSexDict(sample_sheet)
        crDict1 = makeCallRateDict(input.imiss1)
        crDict2 = makeCallRateDict(input.imiss2)
        crDict3 = makeCallRateDict(input.imiss3)
        sexDict = makeSexDict(input.sexcheck)
        contamDict = makeContamDict(input.contam)
        ancestryDict = makeAncestryDict(input.ancestry)
        repDiscordantDict = makeRepDiscordantDict(input.knownConcordance)
        unexpectedRepDict = makeUnexpectedRepDict(input.unknownRep)
        with open(input.sampSheet) as f, open(output.allQc, 'w') as out, open(output.lims, 'w') as limsOut:
            out.write('Sample_ID,CGR_ID,Project_ID,Subject_ID,Expected_Sex,Predicted_Sex,SexMatch,ChrX_Het,AFR,EUR,ASN,Ancestry,Contamination_Rate,Call_Rate_Initial,Call_Rate_1_filter,Call_Rate_1,Call_Rate_2_filter,Call_Rate_2\n')
            limsOut.write(subject_id_to_use + ',Sample ID,Project-Sample ID,Call Rate,Low Call Rate,Contaminated,Sex Discordant,Expected Replicate Discordance,Unexpected Replicate\n')
            head = f.readline()
            while 'SentrixBarcode_A' not in head and head != '':
                head = f.readline()
            if 'SentrixBarcode_A' not in head:
                print('Sample sheet not formatted correctly')
                sys.exit(1)
            head_list = head.rstrip().split(',')
            expectedSexCol = None
            subjectIdCol = None
            projCol = None
            projSampIdCol = None
            for i in range(len(head_list)):
                if head_list[i] == expected_sex_col_name:
                    expectedSexCol = i
                elif head_list[i] == subject_id_to_use:
                    subjectIdCol = i
                elif head_list[i] == 'Project':
                    projCol = i
                elif head_list[i] == 'Project-Sample ID':
                    projSampIdCol = i
            if expectedSexCol == None:
                print('Expected Sex not found in sample sheet')
                sys.exit(1)
            if subjectIdCol == None:
                print('Subject ID not found in sample sheet')
                sys.exit(1)
            if not projCol:
                print('No column named "Project" found')
                sys.exit(1)
            if not projSampIdCol:
                print('No column named "Project-Sample ID" found')
            line = f.readline()
            while line != '':
                if line.strip():
                    line_list = line.rstrip().split(',')
                    samp = line_list[0]
                    chipId = line_list[1] + '_' + line_list[2]
                    proj = line_list[projCol]
                    subId = line_list[subjectIdCol]
                    expectedSex = line_list[expectedSexCol]
                    projSampId = line_list[projSampIdCol]
                    cgrId = samp.split('_')[0]
                    if not crDict1.get(samp):
                        cr1 = 'NA'
                    else:
                        cr1 = crDict1[samp]
                    if not crDict2.get(samp):
                        filt1 = 'Y'
                        cr2 = 'NA'
                    else:
                        filt1 = 'N'
                        cr2 = crDict2[samp]
                    if not crDict3.get(samp):
                        filt2 = 'Y'
                        cr3 = 'NA'
                    else:
                        filt2 = 'N'
                        cr3 = crDict3[samp]
                    if not sexDict.get(samp):
                        (snpsex, het) = ['NA', 'NA']
                    else:
                        (snpsex, het) = sexDict[samp]
                    if snpsex == 'U' or expectedSex == 'U':
                        sexMatch = 'U'
                    elif snpsex != expectedSex:
                        sexMatch = 'N'
                    else:
                        sexMatch = 'Y'
                    if not ancestryDict.get(samp):
                        (AFR, EUR, ASN, ancestry) = ['NA', 'NA', 'NA', 'NA']
                    else:
                        (AFR, EUR, ASN) = ancestryDict[samp]
                        if AFR > 0.8:
                            ancestry = 'AFR'
                        elif EUR > 0.8:
                            ancestry = 'EUR'
                        elif ASN > 0.8:
                            ancestry = 'ASN'
                        else:
                            ancestry = 'ADMIXED'
                    if cr1 == 'NA' or cr2 == 'NA' or cr3 == 'NA':
                        lowCR = 'TRUE'
                    else:
                        lowCR = 'FALSE'
                    if not contamDict.get(samp):
                        contam = 'NA'
                    else:
                        contam = contamDict[samp]
                    if contam != 'NA' and float(contam) > contam_threshold:
                        isContaminated = 'TRUE'
                    else:
                        isContaminated = 'FALSE'
                    if sexMatch == 'U':
                        sexDiscordant = 'NA'
                    elif sexMatch == 'N':
                        sexDiscordant = 'TRUE'
                    else:
                        sexDiscordant = 'FALSE'
                    if repDiscordantDict.get(samp):
                        expectedRepDisc = 'TRUE'
                    else:
                        expectedRepDisc= 'FALSE'
                    if unexpectedRepDict.get(samp):
                        unexpectedRep = 'TRUE'
                    else:
                        unexpectedRep = 'FALSE'
                    limsOut.write(','.join([subId, samp, projSampId, str(cr1), lowCR, isContaminated, sexDiscordant, expectedRepDisc, unexpectedRep]) + '\n')
                    out.write(','.join([samp, cgrId, proj, subId, expectedSex, snpsex, sexMatch, str(het), str(AFR), str(EUR), str(ASN), ancestry, contam, str(cr1), filt1, str(cr2), filt2, str(cr3)]) + '\n')
                    line = f.readline()
        shell('cp {output.lims} {output.limsDirOut}')





rule thousG_match:
    input:
        'ld_prune/samples.bim'
    output:
        bim = 'ld_prune/samples.vcfStrand.bim',
        remove = 'ld_prune/thousG_rename.remove.txt'
    shell:
        'module load python/2.7.5;python /DCEG/CGF/Bioinformatics/Production/Eric/scripts/thousGmatchAndRemoveDups.py {input} {output.remove}'

rule remove_not_thousG:
    input:
        bim = 'ld_prune/samples.vcfStrand.bim',
        bed = 'ld_prune/samples.bed',
        fam = 'ld_prune/samples.fam',
        remove = 'ld_prune/thousG_rename.remove.txt'
    params:
        outProj = 'plink_thousG_match/samples'
    output:
        'plink_thousG_match/samples.bed',
        'plink_thousG_match/samples.bim',
        'plink_thousG_match/samples.fam'
    shell:
        'module load plink2/1.90b3.32;plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --exclude {input.remove} --make-bed --out {params.outProj}'



rule ped_for_snpweights:
    input:
        'plink_thousG_match/samples.bed',
        'plink_thousG_match/samples.bim',
        'plink_thousG_match/samples.fam'
    params:
        inProj = 'plink_thousG_match/samples',
        outProj = 'snpweights/samples'
    output:
        'snpweights/samples.ped',
        'snpweights/samples.map'
    shell:
        'module load plink2/1.90b3.32;plink --bfile {params.inProj} --maf 0.05 --recode --out {params.outProj}'

rule convert_eigen:
    input:
        ped = 'snpweights/samples.ped',
        map = 'snpweights/samples.map'
    output:
        par = 'snpweights/convertEigen.par',
        gen = 'snpweights/samples.eigenstratgeno',
        snp = 'snpweights/samples.snp',
        ind = 'snpweights/samples.ind'
    run:
        parTxt = 'genotypename: ' + input.ped + '\n'
        parTxt += 'snpname: ' + input.map + '\n'
        parTxt += 'indivname: ' + input.ped + '\n'
        parTxt += 'outputformat: EIGENSTRAT\n'
        parTxt += 'genooutfilename: ' + output.gen + '\n'
        parTxt += 'snpoutfilename: ' + output.snp + '\n'
        parTxt += 'indoutfilename: ' + output.ind + '\n'
        parTxt += 'familynames: NO\n'
        with open(output.par, 'w') as out:
            out.write(parTxt)
        shell('module load eigensoft/6.0.1;convertf -p {output.par}')


rule snpweights:
    input:
        gen = 'snpweights/samples.eigenstratgeno',
        snp = 'snpweights/samples.snp',
        ind = 'snpweights/samples.ind'
    params:
        weight = '/DCEG/CGF/Bioinformatics/Production/Eric/software/SNPweights2.1/snpwt.CO',
        software = '/DCEG/CGF/Bioinformatics/Production/Eric/software/SNPweights2.1/bin/inferanc'
    output:
        par = 'snpweights/SNPWEIGHTS.par',
        predpcoutput = 'snpweights/samples.snpweights'
    run:
        parTxt = 'geno: ' + input.gen + '\n'
        parTxt += 'snp: ' + input.snp + '\n'
        parTxt += 'ind: ' + input.ind + '\n'
        parTxt += 'snpwt: ' + params.weight + '\n'
        parTxt += 'predpcoutput: ' + output.predpcoutput + '\n'
        with open(output.par, 'w') as out:
            out.write(parTxt)
        shell('{params.software} -p {output.par}')


rule gtc_to_adpc:
    input:
        gtc = getGtc,
        manifest = illumina_manifest_file
    output:
        adpc = 'contam/{sample_id}.adpc.bin',
        txt = 'contam/{sample_id}.adpc.bin.numSnps.txt'
    shell:
        'module load python/2.7.5;python /DCEG/CGF/Bioinformatics/Production/Eric/scripts/gtc2adpc.py {input.gtc} {input.manifest} {output.adpc}'


rule contam_test:
    input:
        adpc = 'contam/{sample_id}.adpc.bin',
        txt = 'contam/{sample_id}.adpc.bin.numSnps.txt'
    params:
        software = '/DCEG/CGF/Bioinformatics/Production/BZhu/verifyIDintensity/verifyIDintensity/verifyIDintensity'
    output:
        'contam/{sample_id}.contam.out'
    run:
        with open(input.txt) as f:
            line = f.readline()
        markers = line.strip()
        shell('{params.software} -m ' + markers + ' -n 1 -v -p -i {input.adpc} > {output}')

rule combine_contam:
    input:
        expand('contam/{sample_id}.contam.out', sample_id = SAMPLE_IDS)
    output:
        'all_contam/contam.csv'
    run:
        with open(output[0], 'w') as out:
            out.write('ID,%Mix,LLK,LLK0\n')
            for i in input:
                samp = os.path.basename(i).split('.')[0]
                with open(i) as f:
                    head = f.readline()
                    head = f.readline()
                    line = f.readline()
                    line_list = line.split()
                    line_list[0] = samp
                    out.write(','.join(line_list) + '\n')


