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


plinkIn = config['plink_genotype_file'][:-4]

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
adpc_file = config['adpc_file']


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
            inbreed = line_list[-1]
            sexDict[samp] = (snpsex, inbreed)
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
    include: 'modules/Snakefile_gtc_contam'

elif config['plink_genotype_file'][-4:] == '.ped':
    include: 'modules/Snakefile_ped_preprocess'
    include: 'modules/Snakefile_plink_contam'

elif config['plink_genotype_file'][-4:] == '.bed':
    include: 'modules/Snakefile_bed_preprocess'
    include: 'modules/Snakefile_plink_contam'


include: 'modules/Snakefile_plink_stats_filters'
include: 'modules/Snakefile_replicate_concordance'
include: 'modules/Snakefile_sample_qc_report'
include: 'modules/Snakefile_ancestry'


rule all:
    input:
        'all_sample_qc.csv',
        'concordance/KnownReplicates.csv',
        'concordance/UnknownReplicates.csv',
        'snpweights/samples.snpweights',
        'all_contam/contam.csv'
