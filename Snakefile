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
gtc_dir = config['gtc_dir']
input_sample_sheet = config['sample_sheet']
(outName, sampSheetDate) = os.path.basename(input_sample_sheet)[:-4].split('_AnalysisManifest_')
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



POPS = ['ADMIXED_EUR', 'ADMIXED_ASN', 'ADMIXED_AFR', 'ASN_EUR', 'AFR_EUR', 'AFR_ASN', 'AFR_ASN_EUR', 'EUR', 'ASN', 'AFR']
PCs = ['PC1_PC2', 'PC2_PC3', 'PC3_PC4', 'PC4_PC5', 'PC5_PC6', 'PC6_PC7', 'PC7_PC8', 'PC8_PC9', 'PC9_PC10']

D = ['plink_start', 'plink_filter_call_rate_1', 'plink_filter_call_rate_2']
FILT = ['_start', '_filter1', '_filter2']

ylimDict = {'_start':('0', '100'), '_filter1':(str(samp_cr_1 * 100), '100'), '_filter2':(str(samp_cr_2 * 100), '100')}



def DictDiff(dict1, dict2):
    diffDict = {}
    for key in dict1.keys():
        if not dict2.get(key):
            diffDict[key] = 1
    return diffDict


def makeSampDict(fam_file):
    sampDict = {}
    with open(fam_file) as f:
        for line in f:
            samp = line.split()[1]
            sampDict[samp] = 1
    return sampDict


def getCountsByCaCo(sampList, CaCoDict):
    countList = [0, 0, 0, 0]
    for s in sampList:
        CaCo = int(CaCoDict[s])
        countList[CaCo] += 1
    tot = sum(countList)
    (controls, cases, qc, other) = countList
    return [str(controls), str(cases), str(qc), str(other), str(tot)]


def makeSampleToCaCoDict(SampleSheet):
    sampToCaCoDict = {}
    with codecs.open(SampleSheet,"r",encoding='utf-8', errors='ignore') as f:
        head = f.readline()
        while 'SentrixBarcode_A' not in head and head != '':
            head = f.readline()
        if 'SentrixBarcode_A' not in head:
            print('Sample sheet not formatted correctly')
            sys.exit(1)
        head_list = head.rstrip().split(',')
        CaCoCol = None
        sampGroupCol = None
        for i in range(len(head_list)):
            if head_list[i] == 'Case/Control_Status':
                CaCoCol = i
        for i in range(len(head_list)):
            if head_list[i] == 'Sample_Group':
                sampGroupCol = i
        if CaCoCol == None:
            print('Case/Control_Status not found in sample sheet')
            sys.exit(1)
        if sampGroupCol == None:
            print('Sample_Group not found in sample sheet')
            sys.exit(1)
        line = f.readline()
        while line != '':
            if line.strip():
                line_list = line.rstrip().split(',')
                samp = line_list[0]
                CaCo = line_list[CaCoCol]
                sampGroup = line_list[sampGroupCol]
                if not CaCo.strip():
                    if sampGroup == 'sVALD-001':
                        CaCo = 'QC'
                    else:
                        CaCo = 'NA'
                if CaCo == 'Control':
                    CaCo = '0'
                elif CaCo == 'Case':
                    CaCo = '1'
                elif CaCo == 'QC':
                    CaCo = '2'
                else:
                    CaCo = '3'
                sampToCaCoDict[samp] = CaCo
            line = f.readline()
    return sampToCaCoDict





def makeSubjectToCaCoDict(SampleSheet):
    subToCaCoDict = {}
    with codecs.open(SampleSheet,"r",encoding='utf-8', errors='ignore') as f:
        head = f.readline()
        while 'SentrixBarcode_A' not in head and head != '':
            head = f.readline()
        if 'SentrixBarcode_A' not in head:
            print('Sample sheet not formatted correctly')
            sys.exit(1)
        head_list = head.rstrip().split(',')
        subjectIdCol = None
        CaCoCol = None
        sampGroupCol = None
        for i in range(len(head_list)):
            if head_list[i] == subject_id_to_use:
                subjectIdCol = i
        for i in range(len(head_list)):
            if head_list[i] == 'Sample_Group':
                sampGroupCol = i
        if subjectIdCol == None:
            print('Subject ID not found in sample sheet')
            sys.exit(1)
        for i in range(len(head_list)):
            if head_list[i] == 'Case/Control_Status':
                CaCoCol = i
        if CaCoCol == None:
            print('Case/Control_Status not found in sample sheet')
            sys.exit(1)
        if sampGroupCol == None:
            print('Sample_Group not found in sample sheet')
            sys.exit(1)
        line = f.readline()
        while line != '':
            if line.strip():
                line_list = line.rstrip().split(',')
                subId = line_list[subjectIdCol]
                CaCo = line_list[CaCoCol]
                sampGroup = line_list[sampGroupCol]
                if not CaCo.strip():
                    if sampGroup == 'sVALD-001':
                        CaCo = 'QC'
                    else:
                        CaCo = 'NA'
                if CaCo == 'Control':
                    CaCo = '0'
                elif CaCo == 'Case':
                    CaCo = '1'
                elif CaCo == 'QC':
                    CaCo = '2'
                else:
                    CaCo = '3'
                subToCaCoDict[subId] = CaCo
            line = f.readline()
    return subToCaCoDict



def MakeSampToSubDict(sample_to_sub_file):
    sampToSubDict = {}
    with open(sample_to_sub_file) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            (sub, samp) = line.rstrip().split(',')
            if samp != 'NA':
                sampToSubDict[samp] = sub
            line = f.readline()
    return sampToSubDict

def MakeRelatedDict(ibd_file, sample_to_sub_file, sub_fam_file, relatedThresh = float(config['pi_hat_threshold'])):
    sampToSubDict = MakeSampToSubDict(sample_to_sub_file)
    subToKeepDict = {}
    with open(sub_fam_file) as f:
        for line in f:
            sub = line.split()[1]
            subToKeepDict[sub] = 1
    relatedDict = {}
    with open(ibd_file) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            line_list = line.split()
            samp1 = line_list[1]
            samp2 = line_list[3]
            piHat = float(line_list[9])
            if sampToSubDict.get(samp1) and sampToSubDict.get(samp2) and piHat > relatedThresh:
                sub1 = sampToSubDict[samp1]
                sub2 = sampToSubDict[samp2]
                if subToKeepDict.get(sub1) and subToKeepDict.get(sub2):
                    if not relatedDict.get(sub1):
                        relatedDict[sub1] = [sub2]
                    else:
                        relatedDict[sub1].append(sub2)
                    if not relatedDict.get(sub2):
                        relatedDict[sub2] = [sub1]
                    else:
                        relatedDict[sub2].append(sub1)
            line = f.readline()
    return relatedDict






def makeSampSheetDict(SampleSheet, headers):
    '''
    (str, list) -> dict
    headers is a list of header names in the SampleSheet
    headers[0] will be the key to the dict and the rest will be values
    headers[0] is usually "Sample_ID"
    '''
    def getValList(sampSheetLine, headers, headerToColDict, colDict):
        colToValDict = {}
        line_list = sampSheetLine.rstrip().split(',')
        for i in range(len(line_list)):
            if colDict.get(i):
                colToValDict[i] = line_list[i]
        valList = []
        for h in headers:
            col = headerToColDict[h]
            valList.append(colToValDict[col])
        return valList

    sampSheetDict = {}
    headerToColDict = {}
    colDict = {}
    with codecs.open(SampleSheet,"r",encoding='utf-8', errors='ignore') as f:
        head = f.readline()
        while 'SentrixBarcode_A' not in head and head != '':
            head = f.readline()
        if 'SentrixBarcode_A' not in head:
            print('Sample sheet not formatted correctly')
            sys.exit(1)
        head_list = head.rstrip().split(',')
        for i in range(len(head_list)):
            if head_list[i] in headers:
                myHead = head_list[i]
                headerToColDict[myHead] = i
                colDict[i] = 1
        for h in headers:
            if headerToColDict.get(h) == None:
                print(h + ' not found in sample sheet')
                sys.exit(1)
        line = f.readline()
        while line != '':
            if line.strip():
                valList = getValList(line, headers, headerToColDict, colDict)
                sampSheetDict[valList[0]] = valList
            line = f.readline()
    return sampSheetDict


def getBestSamp(CrSampList, crList):
    if len(CrSampList) == 1:
        return CrSampList[0]
    maxCr = 0.0
    maxSamp = 'NA'
    for i in range(len(CrSampList)):
        samp = CrSampList[i]
        CR = float(crList[i])
        if CR > maxCr:
            maxCr = CR
            maxSamp = samp
    return maxSamp


def makeControlDict(SampleSheet):
    controlDict = {}
    with codecs.open(SampleSheet,"r",encoding='utf-8', errors='ignore') as f:
        head = f.readline()
        while 'SentrixBarcode_A' not in head and head != '':
            head = f.readline()
        if 'SentrixBarcode_A' not in head:
            print('Sample sheet not formatted correctly')
            sys.exit(1)
        head_list = head.rstrip().split(',')
        sampGroupCol = None
        subjectIdCol = None
        for i in range(len(head_list)):
            if head_list[i] == 'Sample_Group':
                sampGroupCol = i
            elif head_list[i] == subject_id_to_use:
                subjectIdCol = i
        if sampGroupCol == None:
            print('Sample_Group not found in sample sheet')
            sys.exit(1)
        if subjectIdCol == None:
            print('subject ID col not found in sample sheet')
            sys.exit(1)
        line = f.readline()
        while line != '':
            if line.strip():
                line_list = line.rstrip().split(',')
                sampId = line_list[0]
                subId = line_list[subjectIdCol]
                sampGroup = line_list[sampGroupCol]
                if sampGroup == 'sVALD-001':
                    controlDict[subId] = 1
            line = f.readline()
    return controlDict


def classify_ancestry(labels,x,threshold):
    '''
    This code is modified from GLU admix.py
    https://github.com/bioinformed/glu-genetics/blob/master/glu/modules/struct/admix.py#L345
    An individual is considered of a given ancestry based on the supplied
    labels and estimated admixture coefficients if their coefficient is
    greater than a given threshold.
    Otherwise, an individual who has no single estimated admixture coefficient
    that meets the specified threshold then one of two behaviors result.  If
    only one population group exceeds 1-threshold then the ancestry is deemed
    'ADMIXED' for that population.  Otherwise, a list of populations with
    estimated admixture above 1-threshold is returned.
    '''
    if len(labels) != len(x):
        print('for classify_ancestry labels and x need to be same length')
        sys.exit(1)
    popset = set()
    cmax = -1
    for i in range(len(labels)):
        pop = labels[i]
        coeff = float(x[i])
        if coeff >= 1 - threshold:
            popset.add(pop)
            cmax = max(cmax,coeff)
    if len(popset) == 1 and cmax < threshold:
        ipop = 'ADMIXED_%s' % popset.pop()
    else:
        ipop = '_'.join(sorted(popset))
    return ipop


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



def makeIdatIntensDict(intensFile):
    intensDict = {}
    with open(intensFile) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            line_list = line.rstrip().split(',')
            samp = line_list[0]
            intens = line_list[2]
            intensDict[samp] = intens
            line = f.readline()
    return intensDict



def makeAncestryDict(snpWeightsCsvFile):
    ancestryDict = {}
    with open(snpWeightsCsvFile) as f:
        head = f.readline()
        line = f.readline()
        while line != '':
            line_list = line.rstrip().split(',')
            samp = line_list[0]
            (AFR, EUR, ASN, ancestry) = line_list[-4:]
            ancestryDict[samp] = (AFR, EUR, ASN, ancestry)
            line = f.readline()
    return ancestryDict


def makeChipIdToSampDict(SampleSheet, gtc_dir):
    allSampleIds = []
    idats = []
    noIdats = []
    SAMPLE_IDS = []
    chipIdToGtcDict= {}
    if gtc_dir == 'None':
        gtcFiles = []
    else:
        gtcFiles = []
        for filename in glob.iglob(gtc_dir + '/**/*.gtc', recursive=True):
            gtcFiles.append(filename)
        for gtc in gtcFiles:
            chipId = os.path.basename(gtc)[:-4]
            chipIdToGtcDict[chipId] = gtc
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
                
                if gtc_dir == 'None' and os.path.isfile(idatBase + '.gtc'):
                    gtcFiles.append(idatBase + '.gtc')
                    SAMPLE_IDS.append(sampId)
                    sampIdToGtcDict[sampId] = idatBase + '.gtc'
                elif chipIdToGtcDict.get(chipId):
                    SAMPLE_IDS.append(sampId)
                    sampIdToGtcDict[sampId] = chipIdToGtcDict[chipId]
                if os.path.isfile(idatBase + '_Red.idat') and os.path.isfile(idatBase + '_Grn.idat'):
                    idats.append(idatBase)
                else:
                    noIdats.append(idatBase)
                chipIdToSampDict[chipId] = sampId
            line = f.readline()
    return (chipIdToSampDict, sampIdToGtcDict, sampIdToProjDict, allSampleIds, idats, noIdats, SAMPLE_IDS, gtcFiles, chipIdToSampDict)

(chipIdToSampDict, sampIdToGtcDict, allSampleIds, sampIdToProjDict, idats, noIdats, SAMPLE_IDS, gtcFiles, chipIdToSampDict) = makeChipIdToSampDict(sample_sheet, gtc_dir)


def getGtc(wildcards):
    sampId = wildcards.sample_id
    return sampIdToGtcDict[sampId]


idatBaseDict = {}
for i in idats:
    idatBase = os.path.basename(i)
    idatBaseDict[idatBase] = i

def getRedIdat(wildcards):
    idatBase = wildcards.idatBase
    return idatBaseDict[idatBase] + '_Red.idat'


def getGreenIdat(wildcards):
    idatBase = wildcards.idatBase
    return idatBaseDict[idatBase] + '_Grn.idat'


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
include: 'modules/Snakefile_for_lab'
include: 'modules/Snakefile_idat_intensity'
include: 'modules/Snakefile_identifiler'
include: 'modules/Snakefile_remove_qc_failures'
include: 'modules/Snakefile_subject_level'
include: 'modules/Snakefile_remove_related'
include: 'modules/Snakefile_split_pop'
include: 'modules/Snakefile_autosomal_het'
include: 'modules/Snakefile_subject_ancestry'
include: 'modules/Snakefile_pca'
include: 'modules/Snakefile_HWP'
include: 'modules/Snakefile_plot_completion'
include: 'modules/Snakefile_plot_sex'
include: 'modules/Snakefile_subject_qc_fail'
include: 'modules/Snakefile_count_exclusions'

localrules: summary_stats

rule all:
    input:
        'all_sample_qc.csv',
        'summary_stats.txt',
        'all_sample_idat_intensity/idat_intensity.csv',
        'concordance/KnownReplicates.csv',
        'concordance/UnknownReplicates.csv',
        'snpweights/samples.snpweights.csv',
        'all_contam/contam.csv',
        'files_for_lab/' + outName + '_all_sample_qc_' + sampSheetDate + '.csv',
        'files_for_lab/' + outName + '_KnownReplicates_' + sampSheetDate + '.csv',
        'files_for_lab/' + outName + '_UnknownReplicates_' + sampSheetDate + '.csv',
        'files_for_lab/' + outName + '_LimsUpload_' + sampSheetDate + '.csv',
        'files_for_lab/' + outName + '_Identifiler_' + sampSheetDate + '.csv',
        'subject_level/subjects_qc.imiss',
        'ibd/unrelated_subjects.genome',
        'ancestry/subjects.ancestry.png',
        expand('pca/{pop}_subjects.{pc}.png', pop = POPS, pc = PCs),
        expand('HWP/{pop}_subjects_qc.hwe', pop = POPS),
        expand('{d}/samples{filt}.completion.png', zip, d = D, filt = FILT),
        expand('autosomal_heterozygosity/{pop}_subjects_qc.het.png', pop = POPS),
        'sex_plot/sex.png',
        'final_qc_subject_level/subjects.bed',
        'counts/exclusion_counts.csv'

