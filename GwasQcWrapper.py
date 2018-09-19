#!/usr/bin/env python2.7
'''
Gwas QC Pipeline written by Eric Karlins
karlinser@mail.nih.gov
This script is a wrapper for Snakemake.  The pipeline is run through rules in the Snakefile.
'''


import sys
import math
import os
import subprocess
import shutil
import argparse
import time
import glob
from os import stat
from pwd import getpwuid

def makeQsub(qsubFile, qsubText):
    '''
    Make a qsub file in the qsubDir
    '''
    with open(qsubFile, 'w') as output:
        output.write('#!/bin/bash\n\n')
        output.write(qsubText)


def runQsub(qsubFile, proj, queue):
    '''
    (str, int) -> None
    '''
    qsubHeader = qsubFile[:-3]
    user = subprocess.check_output('whoami').rstrip('\n')
    qsubCall = ['qsub', '-M', user + '@mail.nih.gov', '-m', 'beas', '-q', queue, '-o', qsubHeader + '.stdout', '-e', qsubHeader + '.stderr', '-N', 'GwasQC.' + proj, '-S', '/bin/sh', qsubFile]
    retcode = subprocess.call(qsubCall)




def getNumSamps(sampleSheet):
    samps = 0
    with open(sampleSheet) as f:
        head = f.readline()
        while 'SentrixBarcode_A' not in head and head != '':
            head = f.readline()
        if 'SentrixBarcode_A' not in head:
            print('Sample sheet not formatted correctly')
            sys.exit(1)
        line = f.readline()
        while line != '':
            samps += 1
            line = f.readline()
    return samps




def makeClusterConfig(outDir, queue):
    q_list = queue.split(',')
    new_q_list = []
    for q in q_list:
        if q != 'all.q':
            new_q_list.append(q)
    if not new_q_list:
        print('There must be a queue other than all.q in your queue list')
        sys.exit(1)
    new_q = ','.join(new_q_list)
    with open(outDir + '/cluster.yaml', 'w') as output:
        output.write('__default__:\n')
        output.write('    q: ' + queue + '\n')
        output.write('contam_test:\n')
        output.write('    q: ' + new_q + '\n')
        output.write('merge_sample_peds:\n')
        output.write('    q: ' + new_q + '\n')
        output.write('plink_ibd:\n')
        output.write('    q: ' + new_q + '\n')


def makeConfig(outDir, plink_genotype_file, snp_cr_1, samp_cr_1, snp_cr_2, samp_cr_2, ld_prune_r2, maf_for_ibd, sample_sheet,
               subject_id_to_use, ibd_pi_hat_cutoff, dup_concordance_cutoff, illumina_manifest_file, expected_sex_col_name, numSamps, lims_output_dir, 
               contam_threshold, adpc_file, gtc_dir, remove_contam, remove_sex_discordant, remove_rep_discordant, remove_unexpected_rep, pi_hat_threshold,
               autosomal_het_thresh, minimum_pop_subjects, control_hwp_thresh, word_doc_template, contam_pop, strand):
    '''
    (str, str, str) -> None
    '''
    paths = os.listdir(outDir)
    shutil.copyfile(sample_sheet, outDir + '/IlluminaSampleSheet.csv')
    if 'config.yaml' in paths:
        start = getStartTime(outDir + '/config.yaml')
    else:
        start = time.ctime()
    with open(outDir + '/config.yaml', 'w') as output:
        output.write('plink_genotype_file: ' + str(plink_genotype_file) + '\n')
        output.write('illumina_manifest_file: ' + str(illumina_manifest_file) + '\n')
        output.write('snp_cr_1: ' + str(snp_cr_1) + '\n')
        output.write('samp_cr_1: ' + str(samp_cr_1) + '\n')
        output.write('snp_cr_2: ' + str(snp_cr_2) + '\n')
        output.write('samp_cr_2: ' + str(samp_cr_2) + '\n')
        output.write('ld_prune_r2: ' + str(ld_prune_r2) + '\n')
        output.write('maf_for_ibd: ' + str(maf_for_ibd) + '\n')
        output.write('sample_sheet: ' + sample_sheet + '\n')
        output.write('subject_id_to_use: ' + subject_id_to_use + '\n')
        output.write('ibd_pi_hat_cutoff: ' + str(ibd_pi_hat_cutoff) + '\n')
        output.write('dup_concordance_cutoff: ' + str(dup_concordance_cutoff) + '\n')
        output.write('expected_sex_col_name: ' + expected_sex_col_name + '\n')
        output.write('num_samples: ' + str(numSamps) + '\n')
        output.write('lims_output_dir: ' + lims_output_dir + '\n')
        output.write('contam_threshold: ' + str(contam_threshold) + '\n')
        output.write('adpc_file: ' + str(adpc_file) + '\n')
        output.write('gtc_dir: ' + str(gtc_dir) + '\n')
        output.write('remove_contam: "' + remove_contam + '"\n')
        output.write('remove_sex_discordant: "' + remove_sex_discordant + '"\n')
        output.write('remove_rep_discordant: "' + remove_rep_discordant + '"\n')
        output.write('remove_unexpected_rep: "' + remove_unexpected_rep + '"\n')
        output.write('pi_hat_threshold: ' + str(pi_hat_threshold) + '\n')
        output.write('autosomal_het_thresh: ' + str(autosomal_het_thresh) + '\n')
        output.write('minimum_pop_subjects: ' + str(minimum_pop_subjects) + '\n')
        output.write('control_hwp_thresh: ' + str(control_hwp_thresh) + '\n')
        output.write('doc_template: ' + word_doc_template + '\n')
        output.write('contam_pop: ' + contam_pop + '\n')
        output.write('strand: ' + strand + '\n')
        output.write('start_time: ' + start + '\n')


def getStartTime(configFile):
    with open(configFile) as f:
        line = f.readline()
        while not line.startswith('start_time'):
            line = f.readline()
        return line.rstrip('\n').split(': ')[1]
    return time.ctime()


def get_args():
    '''
    return the arguments from parser
    '''
    parser = argparse.ArgumentParser()
    requiredArgs = parser.add_argument_group('Required Arguments')
    oneOfTheseRequired = parser.add_argument_group('Exactly one of these arguments is required')
    requiredWithDefaults = parser.add_argument_group('Required arguments with default settings')
    oneOfTheseRequired.add_argument('-p', '--path_to_plink_file', type=str, required=False, help='Full path to either PLINK ped or bed to use as input.\n\
                        need either this or gtc file project directory -g')
    requiredArgs.add_argument('-d', '--directory_for_output', type=str, help='REQUIRED. Full path to the base directory for the Gwas QC pipeline output.  Defaults to /DCEG/CGF/GWAS/Scans/GSA_Lab_QC/Manifest#/builds/QC_v#_date')
    requiredWithDefaults.add_argument('--snp_cr_1', type=float, default= 0.80, help='REQUIRED. SNP call rate filter 1.  default= 0.80')
    requiredWithDefaults.add_argument('--samp_cr_1', type=float, default= 0.80, help='REQUIRED. Sample call rate filter 1.  default= 0.80')
    requiredWithDefaults.add_argument('--snp_cr_2', type=float, default= 0.95, help='REQUIRED. SNP call rate filter 2.  default= 0.95')
    requiredWithDefaults.add_argument('--samp_cr_2', type=float, default= 0.95, help='REQUIRED. Sample call rate filter 2.  default= 0.95')
    requiredWithDefaults.add_argument('--ld_prune_r2', type=float, default= 0.10, help='REQUIRED. r-squared cutoff for ld pruning of SNPs to use for IBD and concordance.  default= 0.10')
    requiredWithDefaults.add_argument('--maf_for_ibd', type=float, default= 0.20, help='REQUIRED. MAF cutoff of SNPs to use for IBD and concordance.  default= 0.20')
    requiredArgs.add_argument('-s', '--sample_sheet', type=str, required=True, help='Full path to illimina style sample sheet csv file.')
    requiredArgs.add_argument('--subject_id_to_use', type=str, default= 'LIMS_Individual_ID', help='Name of column in sample sheet that corresponds to subject ID to use.')##I should be able to add a default once this is available
    requiredWithDefaults.add_argument('--ibd_pi_hat_cutoff', type=float, default= 0.95, help='REQUIRED. PI_HAT cutoff to call samples replicates.  default= 0.95')##this can be deleted if just using SNP concordance
    requiredWithDefaults.add_argument('--dup_concordance_cutoff', type=float, default= 0.95, help='REQUIRED. SNP concordance cutoff to call samples replicates.  default= 0.95')
    requiredWithDefaults.add_argument('--lims_output_dir', type = str, default = '/DCEG/CGF/Laboratory/LIMS/drop-box-prod/gwas_primaryqc', help='Directory to copy QC file to upload to LIMS')
    requiredWithDefaults.add_argument('--contam_threshold', type=float, default= 0.10, help='REQUIRED. Cutoff to call a sample contaminated.  default= 0.10')
    requiredWithDefaults.add_argument('--contam_pop', type=str, default= 'AF', choices=['AF','EAS_AF','AMR_AF','AFR_AF','EUR_AF','SAS_AF'], help='REQUIRED. 1KG pop to use for freq for contam test. Must be one of:  AF,EAS_AF,AMR_AF,AFR_AF,EUR_AF,SAS_AF.  default= AF')
    parser.add_argument('-i', '--illumina_manifest_file',type=str, default='/DCEG/CGF/Infinium/Resources/Manifests/GSAMD-Files/build37/GSAMD-24v1-0_20011747_A1.bpm', help='Full path to illimina .bpm manifest file. Required for gtc files.')
    requiredWithDefaults.add_argument('--strand', type=str, default= 'TOP', choices=['TOP','FWD'], help='REQUIRED. strand to use for genotypes. Must be one of:  TOP, FWD.  default= TOP')
    parser.add_argument('-a', '--adpc_file', type=str, help='Full path to adpc.bin file. Required for PLINK input.')
    parser.add_argument('-g', '--gtc_dir', type=str, help='Full path to gtc directory to use instead of project directory, which is the default.  Will recursively find gtc files in this directory.')
    requiredArgs.add_argument('--expected_sex_col_name', type=str, default='Expected_Sex', help='Name of column in sample sheet that corresponds to expected sex of sample.')##I should be able to add a default once this is available
    requiredWithDefaults.add_argument('-q', '--queue', type=str, default='all.q,seq-alignment.q,seq-calling.q,seq-calling2.q,seq-gvcf.q,bigmem.q', help='OPTIONAL. Queue on cgemsiii to use to submit jobs.  Defaults to all of the seq queues and all.q if not supplied.  default="all.q,seq-alignment.q,seq-calling.q,seq-calling2.q,seq-gvcf.q,bigmem.q"')
    requiredWithDefaults.add_argument('--remove_contam', type=str, default='YES', help='REQUIRED. If "YES" contaminated samples will be removed prior to sample to subject transformation.  Defaults to "YES"')
    requiredWithDefaults.add_argument('--remove_sex_discordant', type=str, default='YES', help='REQUIRED. If "YES" sex discordant samples will be removed prior to sample to subject transformation.  Defaults to "YES"')
    requiredWithDefaults.add_argument('--remove_rep_discordant', type=str, default='YES', help='REQUIRED. If "YES" known replicate discordant samples will be removed prior to sample to subject transformation.  Defaults to "YES"')
    requiredWithDefaults.add_argument('--remove_unexpected_rep', type=str, default='YES', help='REQUIRED. If "YES" all unexpected replicate samples will be removed prior to sample to subject transformation.  Defaults to "YES"')
    requiredWithDefaults.add_argument('--pi_hat_threshold', type=float, default= 0.20, help='REQUIRED. PI_HAT cutoff to call subjects related.  default= 0.20, to remove 2nd degree relatives or higher.')
    requiredWithDefaults.add_argument('--minimum_pop_subjects', type=int, default= 100, help='REQUIRED. Number of subjects needed in order to analyze a population.  default= 100.')
    requiredWithDefaults.add_argument('--control_hwp_thresh', type=int, default= 50, help='REQUIRED. Number of controls needed in order to just use controls for HWP.  default= 50.')
    requiredWithDefaults.add_argument('--autosomal_het_thresh', type=float, default= 0.10, help='REQUIRED. F coefficient from autosomal heterozygosity check cutoff for subject removal.  default= 0.10')
    requiredWithDefaults.add_argument('-w', '--word_doc_template', type=str, default='/DCEG/CGF/Bioinformatics/Production/Eric/refFiles/GwasQcPipeline_Template_interim.docx', help='REQUIRED. Template to use to generate QC report word docx')
    parser.add_argument('-u', '--unlock_snakemake', action='store_true', help='OPTIONAL. If pipeline was killed unexpectedly you may need this flag to rerun')
    parser.add_argument('-f', '--finish', action='store_true', help='OPTIONAL. Use with -d option to restart pipeline or update with new features without making new config file, etc.')
    args = parser.parse_args()
    return args





def getOutDir(sampleSheet, baseDir = '/DCEG/CGF/GWAS/Scans/GSA_Lab_QC/'):
    '''
    return the name of the output directory and make it if it doesn't exist
    '''
    (sr, sampSheetDate) = os.path.basename(sampleSheet)[:-4].split('_AnalysisManifest_')
    srDir = baseDir + sr + '_' + sampSheetDate
    if not os.path.isdir(srDir):
        os.mkdir(srDir)
    if not os.path.isdir(srDir + '/builds'):
        os.mkdir(srDir + '/builds')
    prevBuilds = glob.glob(srDir + '/builds/QC_v*')
    buildVersions = []
    for build in prevBuilds:
        d = os.path.basename(build)
        v = d.split('_')[1][1:]
        if v.isdigit():
            buildVersions.append(int(v))
    if buildVersions:
        vers = sorted(buildVersions)[-1] + 1
    else:
        vers = 1
    date = time.strftime("%m%d%Y")
    outDir = srDir + '/builds/QC_v' + str(vers) + '_' + date
    if not os.path.isdir(outDir):
        os.mkdir(outDir)
    return outDir


def main():
    scriptDir = os.path.dirname(os.path.abspath(__file__))
    args = get_args()
    sampSheetOwner = getpwuid(stat(args.sample_sheet).st_uid).pw_name
#I want to add this check that the file was created by LIMS, but we can't roll it out just yet
#    if sampSheetOwner != 'cgflims':
#        print('Pipeline requires AnalysisManifest to be generated by LIMS.  File owner is ' + sampSheetOwner)
#        print('Exiting.  Generate a file using LIMS and rerun.')
#        sys.exit(1)
    outDir = args.directory_for_output
    if not outDir:
        if args.finish:
            print('you need to use the -d option if using the -f option.')
            sys.exit(1)
        outDir = getOutDir(args.sample_sheet)
    if outDir[0] != '/':
        print('-d argument must be full path to working directory.  Relative paths will not work.')
        sys.exit(1)
    paths = os.listdir(outDir)
    if 'logs' not in paths:
        os.mkdir(outDir + '/logs')
    if 'modules' not in paths:
        os.mkdir(outDir + '/modules')
    if 'figures' not in paths:
        os.mkdir(outDir + '/figures')
    if not args.path_to_plink_file:
        if not args.illumina_manifest_file:
            print('--illumina_manifest_file is required for gtc files.')
            sys.exit(1)
        plinkPedOrFam = None
    else:
        plinkFile = args.path_to_plink_file
        plinkDir = os.path.dirname(plinkFile)
        if not args.adpc_file:
            args.adpc_file = plinkDir + '/Data/adpc.bin'
        if plinkFile[-4:] == '.ped':
            plinkPedOrFam = plinkFile
        elif plinkFile[-4:] == '.bed':
            plinkPedOrFam = plinkFile[:-4] + '.fam'
        else:
            print('Unrecognized PLINK file format.')
            sys.exit(1)
    shutil.copy2(scriptDir + '/Snakefile', outDir + '/Snakefile')
    moduleFiles = glob.glob(scriptDir + '/modules/*')
    for f in moduleFiles:
        shutil.copy2(f, outDir + '/modules')
    figFiles = glob.glob(scriptDir + '/figures/*')
    for f in figFiles:
        shutil.copy2(f, outDir + '/figures')
    if args.finish:
        numSamps = getNumSamps(outDir + '/IlluminaSampleSheet.csv')
    if not args.finish:
        numSamps = getNumSamps(args.sample_sheet)
        makeConfig(outDir, args.path_to_plink_file, args.snp_cr_1, args.samp_cr_1, args.snp_cr_2, args.samp_cr_2, args.ld_prune_r2, args.maf_for_ibd, args.sample_sheet,
               args.subject_id_to_use, args.ibd_pi_hat_cutoff, args.dup_concordance_cutoff, args.illumina_manifest_file, args.expected_sex_col_name, numSamps, args.lims_output_dir, args.contam_threshold,
               args.adpc_file, args.gtc_dir, args.remove_contam, args.remove_sex_discordant, args.remove_rep_discordant, args.remove_unexpected_rep, args.pi_hat_threshold, args.autosomal_het_thresh,
               args.minimum_pop_subjects, args.control_hwp_thresh, args.word_doc_template, args.contam_pop, args.strand)
        makeClusterConfig(outDir, args.queue)
    qsubTxt = 'cd ' + outDir + '\n'
    qsubTxt += 'module load sge\n'
    qsubTxt += 'module load python3/3.5.1\n'
    qsubTxt += 'module load R/3.4.0\n'
    qsubTxt += 'module load plink2/1.90b5\n'
    qsubTxt += 'module load python/2.7.5\n'
    qsubTxt += 'module load eigensoft/7.2.1\n'
    qsubTxt += 'module load R/3.4.0\n'
    if args.unlock_snakemake:
        qsubTxt += 'snakemake --unlock\n'
    qsubTxt += 'snakemake --rerun-incomplete --cluster "qsub -V -q {cluster.q} -pe by_node {threads} '
    qsubTxt += '-o ' + outDir + '/logs/ -e ' + outDir + '/logs/" --jobs 4000 --cluster-config cluster.yaml --latency-wait 300\n'
    makeQsub(outDir + '/GwasQcPipeline.sh', qsubTxt)
    runQsub(outDir + '/GwasQcPipeline.sh', os.path.basename(outDir), 'seq-alignment.q,seq-calling.q,seq-calling2.q,seq-gvcf.q')
    print('GWAS QC Pipeline submitted.  You should receive an email when the pipeline starts and when it completes.')
    print('Your input project has ' + str(numSamps) + ' samples.')


if __name__ == "__main__":
    main()
