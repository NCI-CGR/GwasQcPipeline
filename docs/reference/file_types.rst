File Type Reference
===================

Here I briefly describe the large number of filetypes used throughout the GWAS QC Pipeline.


.. A

adpc.bin
--------

An Illumina file format that can be output by GenomeStudio.
I cannot find a spec file or any details about this format beyond a small reader function used by Picard_.

.. rubric:: Fields

:A allele intensity:
:B allele intensity:
:A normalized allele intensity:
:B normalized allele intensity:
:cluster confidence score:
:called genotype: 0 (AA), 1 (AB), 2 (BB), and 3 (NN)

.. _Picard: https://javadoc.io/static/org.broadinstitute/gatk/4.1.4.1/picard/arrays/illumina/IlluminaAdpcFileWriter.html

abf.txt
-------

This is a simple text file with the ``B`` allele frequencies for each snp.

.. rubric:: Fields

:SNP_ID:
:BAF: (a.k.a ABF)

.. B

bed (``plink --make-bed``)
--------------------------

The ``bed`` file is a binary genotype table with genotype calls at biallelic variants.
It is always accompanied by a ``fam`` and ``bim`` file.

.. epigraph::

    The rest of the file is a sequence of V blocks of N/4 (rounded up) bytes each, where V is the number of variants and N is the number of samples.
    The first block corresponds to the first marker in the .bim file, etc.

    The low-order two bits of a block's first byte store the first sample's genotype code.
    ("First sample" here means the first sample listed in the accompanying .fam file.)
    The next two bits store the second sample's genotype code, and so on for the 3rd and 4th samples.
    The second byte stores genotype codes for the 5th-8th samples, the third byte stores codes for the 9th-12th, etc.

    The two-bit genotype codes have the following meanings:

    00 Homozygous for first allele in .bim file
    01 Missing genotype
    10 Heterozygous
    11 Homozygous for second allele in .bim file

    If N is not divisible by four, the extra high-order bits in the last byte of each block are always zero.

    For example, 0x6c 0x1b 0x01 0xdc 0x0f 0xe7 0x0f 0x6b 0x01

    1. The first three bytes are the magic number (0x6c 0x1b 0x01).
    2. Since there are six samples, each marker block has size 2 bytes (six divided by four, rounded up).
       Thus genotype data for the first marker ('snp1') is stored in the 4th and 5th bytes (0xdc 0x0f).
    3. The 4th byte value of 0xdc is [11][01][11][00] in binary.
       Since the low-order two bits are '00' (right most), the first sample is homozygous for the first allele for this marker listed in the .bim file, which is 'G'.
       The second sample has genotype code '11', which means she's homozygous for the second allele ('A').
       The third sample's code of '01' designates a missing genotype call, and the fourth code of '11' indicates another AA.
    4. The 5th byte value of 0x0f is [00][00][11][11] in binary.
       This indicates that the fifth and sixth samples also have the AA genotype at snp1.
       There is no sample #7 or #8, so the high-order 4 bits of this byte are zero (far left).
    5. The 6th and 7th bytes store genotype data for the second marker ('snp2').
       The 6th byte value of 0xe7 is 11100111 in binary.
       The '11' code for the first sample means that he's homozygous for the second snp2 allele ('2'), the '01' code for the second sample indicates a missing call, the '10' code for the third indicates a heterozygous genotype, and '11' for the fourth indicates another homozygous '2'.
       The 7th byte value of 0x0f indicates that the fifth and sixth samples also have homozygous '2' genotypes.
    6. Finally, the 8th and 9th bytes store genotype data for the third marker ('snp3').
       You can test your understanding of the file format by interpreting this by hand and then comparing to the .ped file above.

    -- https://www.cog-genomics.org/plink2/formats#bed

bim (``plink --make-bed``)
--------------------------

The ``bim`` file describes each marker SNP. It is usually accompanied by a ``fam`` and ``bed`` file.

.. epigraph::
  A text file with no header line, and one line per variant with the following six fields.
  Allele codes can contain more than one character. Variants with negative bp coordinates are ignored by PLINK.

.. rubric:: Fields

:Chromosome code: (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
:Variant identifier:
:Position: in morgans or centimorgans (safe to use dummy value of '0')
:Base-pair coordinate: (1-based; limited to 2^{31}-2)
:Allele 1: (corresponding to clear bits in .bed; usually minor)
:Allele 2: (corresponding to set bits in .bed; usually major)

bpm (Illumina Manifest File)
----------------------------

The Illumina manifest file is a binary file that describes the specific array.
It has information about each probe set and their corresponding SNPs.
This version ``GSAMD-24v1-0_20011747_A1.bpm`` is linked to in the config.
There is also a ``csv`` version that can be downloaded from Illumina.

.. rubric:: Fields

:IlmnID:
:Name:
:IlmnStrand: ``BOT|TOP|PLUS|MINUs``
:SNP: ``[A/B]`` where ``A`` and ``B`` can be ``[ACTGID]``
:AddressA_ID: The location on the array for ``A``
:AlleleA_ProbeSeq: The probe sequence for ``A``
:AddressB_ID: The location on the array for ``B``
:AlleleB_ProbeSeq: The probe sequence for ``B``
:GenomeBuild: Reference genome build number (i.e. 37)
:Chr: Chromosome
:MapInfo:
:Ploidy: Ploidy of the species (i.e. diploid)
:Species: The species (i.e. Homo sapiens)
:Source: The source for the SNP (e.g. 1000genomes, PAGE, ClinVar_ACMG).
:SourceVersion: Version of the source.
:SourceStrand: ``BOT|TOP|PLUS|MINUs``
:SourceSeq:
:TopGenomicSeq:
:BeadSetID:
:Exp_Clusters:
:Intensity_Only:
:RefStrand: Reference strand ``(\+|\-)``

https://support.illumina.com/bulletins/2017/06/how-to-interpret-dna-strand-and-allele-information-for-infinium-.html

.. C

contam.out (verifyIDintensity)
------------------------------

.. warning::
  This output has no official documentation.

Sample level contamination score.

.. rubric:: Fields

:ID: sample_ID
:%MIX: Percent mixture with another sample (i.e. Alpha from Jun et al. 2012)
:LLN: The minimized log-likelihood. *not really sure what this is*
:LLN0: The log-likelihood at alpha = 0 (i.e. no contamination).

csv (Illumina Sample Sheet)
---------------------------

.. warning::
  This output has no official documentation.

I am assuming this is the sample sheet used to run ``GenomeStudio``.

.. rubric:: Fields

:Sample_ID: ``SC\d+_PB\d+_[A-H][0-1][1-9]``
:SentrixBarcode_A: ``\d+``
:SentrixPosition_A: ``R\d{2}C\d{2}``
:Sample_Plate: ``WG\d+-DNA``
:Sample_Well: ``[A-H][0-1][1-9]``
:Sample_Group: ``.*``
:Identifier_Sex: ``(M|F|U)``
:Sample_Name: ``\w{2} \d{4} \d{4}``
:Replicate: ``.*``
:Parent1: ``SC\d+_PB\d+_[A-H][0-1][1-9]``
:Parent2: ``SC\d+_PB\d+_[A-H][0-1][1-9]``
:SR_Subject_ID: ``SI\d+``
:LIMSSample_ID: ``SC\d+``
:Sample_Type: ``[\w\s]+``
:LIMS_Specimen_ID:
:Project: ``GP\d+-IN\d{2}``
:Project-Sample ID: ``PS-(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})-\d{6}``
:Array: ``GSAMD-24v1-0``
:LIMS_Individual_ID: ``I-\d+``
:PI_Subject_ID: ``\w-\d+-\d``
:PI_Study_ID:
:Age:
:Expected_Sex: ``(M|F)``
:Ancestry_S1:
:Ancestry_S2:
:Ancestry_S3:
:POPGROUP:
:Case/Control_Status: ``(Control|Case)``

.. D

.. E

eigenstratgeno (``eigensoft``)
------------------------------

The genotype data for each individual at each SNP. Each column represent an individual sample (same order as ind file). Each row represents a SNP (same order as snp file). Where values are 0-2 or 9.

.. rubric:: Encoding

:0: Zero copies of the reference allele
:1: One copies of the reference allele
:2: Two copies of the reference allele
:9: Missing data

.. F

fam (``plink --make-bed``)
--------------------------

The ``fam`` sample information and accompanies a ``bed`` file.

.. epigraph::
    A text file with no header line, and one line per sample ...

    If there are any numeric phenotype values other than {-9, 0, 1, 2}, the phenotype is interpreted as a quantitative trait instead of case/control status.
    In this case, -9 normally still designates a missing phenotype; use --missing-phenotype if this is problematic.

    If your case/control phenotype is encoded as '0' = control and '1' = case, you'll need to specify --1 to load it properly.

    -- https://www.cog-genomics.org/plink2/formats#fam

.. rubric:: Fields

:Family ID: ('FID')
:Within-family ID: ('IID'; cannot be '0')
:Within-family ID of father: ('0' if father isn't in dataset)
:Within-family ID of mother: ('0' if mother isn't in dataset)
:Sex code: ('1' = male, '2' = female, '0' = unknown)
:Phenotype value: ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

frq (``plink --freq``)
----------------------

The allele frequency report.

.. rubric:: Fields

:CHR: Chromosome
:SNP: Variant identifier
:A1: Allele 1 (usually minor)
:A2: Allele 2 (usually major)
:MAF: Allele 1 frequency
:NCHROBS: Number of allele observations

https://www.cog-genomics.org/plink/1.9/formats#frq

.. G

genome (``plink --genome``)
---------------------------

The IBS/IBD report.

.. rubric:: Fields

:FID1: Family ID for first sample
:IID1: Individual ID for first sample
:FID2: Family ID for second sample
:IID2: Individual ID for second sample
:RT: Relationship type inferred from .fam/.ped file
:EZ: IBD sharing expected value, based on just .fam/.ped relationship
:Z0: P(IBD=0)
:Z1: P(IBD=1)
:Z2: P(IBD=2)
:PI_HAT: Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
:PHE: Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)
:DST: IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)
:PPC: IBS binomial test
:RATIO: HETHET : IBS0 SNP ratio (expected value 2)

.. rubric:: Optional Fields (given ``--genome full``)

:IBS0: Number of nonmissing variants where both allele are different
:IBS1: Number of nonmissing variants where 1 allele is the same
:IBS2: Number of nonmissing variants where both alleles are the same
:HOMHOM: Number of IBS 0 SNP pairs used in PPC test
:HETHET: Number of IBS 2 het/het SNP pairs used in PPC test

gtc (Illumina Genotype Calls)
-----------------------------

The Illumina Infinium genotype (GTC) file format. This format is output from Illumina's genotype calling software (either Autocall or Autoconvert). This is a complicated data format containing large amounts of metadata describing the run, data describing each probe on the array, as well as sample level genotype calls for probe.

Each sample has a GTC file associated with it.

:download:`Format Spec <https://github.com/Illumina/BeadArrayFiles/blob/develop/docs/GTC_File_Format_v5.pdf>`

.. H

hwe (``plink --hardy``)
-----------------------

Exact test results for Hardy-Weinberg equilibrium.

.. rubric:: Fields

:CHR: Chromosome
:SNP: Variant identifier
:TEST: Type of test: one of {'ALL', 'AFF', 'UNAFF', 'ALL(QT)', 'ALL(NP)'}
:A1: Allele 1 (usually minor)
:A2: Allele 2 (usually major)
:GENO: '/'-separated genotype counts (A1 hom, het, A2 hom)
:O(HET): Observed heterozygote frequency
:E(HET): Expected heterozygote frequency
:P: Hardy-Weinberg equilibrium exact test p-value

When the samples are case/control, three separate sets of Hardy-Weinberg equilibrium statistics are computed: one considering both cases and controls, one considering only cases, and one considering only controls. These are distinguished by 'ALL', 'AFF', and 'UNAFF' in the TEST column, respectively. If the phenotype is quantitative or nonexistent instead, there is just one line per variant, labeled 'ALL(QT)' or 'ALL(NP)' respectively.

https://www.cog-genomics.org/plink/1.9/formats#hwe

.. I

idat (Illumina)
---------------

A binary format of intensities. This file includes various types of metadata (i.e., array information, software versions, the type of BeadChip). Data fields are made up of 4 values: ID of each probe on the array, mean intensity, intensity standard deviation, and the number of beads with each probe.


imiss (``plink --missing``)
---------------------------

Sample-based missing data report.

.. rubric:: Fields

:FID: Family ID
:IID: Within-family ID
:MISS_PHENO: [Y/N] Phenotype is missing
:N_MISS: The number of missing genotype calls not including obligatory missing or heterozygous haploids.
:N_GENO: Number of potentially valid calls
:F_MISS: Missing call rate

https://www.cog-genomics.org/plink/1.9/formats#imiss

ind (``eigensoft``)
-------------------

The individual sample description file. Each row represents a sample.

.. rubric:: Fields

:Sample ID: The sample identifier
:gender: The gender of the sample {M, F, U}
:phenotype: The sample phenotype or population information.



.. J

.. K

.. L

lmiss (``plink --missing``)
---------------------------

Variant-based missing data report.

.. rubric:: Fields

:CHR: Chromosome
:SNP: Variant ID
:N_MISS: The number of missing genotype calls no including obligatory missing or heterozygous haploids.
:N_GENO: Number of potentially valid calls
:F_MISS: Missing call rate
:Optional fields: if run with ``--within/--family``

  :CLST: Cluster identifier
  :N_CLST: Cluster size (does not include non-males on chrY)

https://www.cog-genomics.org/plink/1.9/formats#.miss

.. M

map (``plink``)
---------------

The ``map`` file describes the location of variants.

.. epigraph::
    Variant information file accompanying a .ped.

    A text file with no header file, and one line per variant with the following 3-4 fields.

    All lines must have the same number of columns (so either no lines contain the morgans/centimorgans column, or all of them do).

    -- https://www.cog-genomics.org/plink2/formats#map

.. rubric:: Fields

:Chromosome code: PLINK 1.9 also permits contig names here, but most older programs do not.
:Variant identifier:
:Position in morgans or centimorgans: (optional; also safe to use dummy value of '0')
:Base-pair coordinate:

.. N

.. O

.. P

ped (``plink``)
---------------

The ``ped`` format is a text pedigree and genotype table used by PLINK/MERLIN/Haploview.

.. epigraph::
    Original standard text format for sample pedigree information and genotype calls. Normally must be accompanied by a .map file ...

    Contains no header line, and one line per sample with **2V+6 fields** where V is the number of variants.
    The first six fields are the same as *fam* files.
    The seventh and eighth fields are allele calls for the first variant in the .map file ('0' = no call); the 9th and 10th are allele calls for the second variant; and so on.

    -- https://www.cog-genomics.org/plink2/formats#ped

.. rubric:: Fields

:Family ID: ('FID')
:Within-family ID: ('IID'; cannot be '0')
:Within-family ID of father: ('0' if father isn't in dataset)
:Within-family ID of mother: ('0' if mother isn't in dataset)
:Sex code: ('1' = male, '2' = female, '0' = unknown)
:Phenotype value: ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
:Allele A marker_{1..n}:
:Allele B marker_{1..n}:

prune.in (``plink --indep-pairwise``)
-------------------------------------

A list of marker IDs to include after running LD pruning.

prune.out (``plink --indep-pairwise``)
--------------------------------------

A list of marker IDs to exclude after running LD pruning.



.. Q

.. R

.. S

sexcheck (``plink --check-sex``)
--------------------------------

X chromosome based sex sanity checks.

.. rubric:: Fields

:FID: Family ID
:IID: Within-family ID
:PEDSEX: Sex code in input file
:SNPSEX: Imputed sex code (1 = male, 2 = female, 0 = unknown)
:STATUS: [OK/PROBLEM] if PEDSEX and SNPSEX match then OK

.. rubric:: Optional Fields (given ``--ycount`` or ``-y-only``)

:F: Inbreeding coefficient based on X chromosome
:YCOUNT: Number of non-missing genotypes calls on Y chromosome

https://www.cog-genomics.org/plink/1.9/formats#sexcheck

snp (``eigensoft``)
-------------------

A file describing each SNP (one SNP per row).

.. rubric:: Fields

:SNP ID:
:Chromosome: The chromosome number (X: 23, Y: 24, mtDNA: 90, XY: 91)
:Genetic Position: Given in morgans or 0.0 if unknown.
:Physical Position: Given in bases.

.. rubric:: Optional Fields

:Reference Allele:
:Variant Allele: For monomorphic SNPs can be encoded as X (unknown).

snpwt.* (``SNPweights``)
------------------------

.. rubric:: Header Rows

:first: shrinkage for the predicted PCs (1 and 2)
:second: ancestral populations
:third: number of ancestral samples for each population
:fourth: average PCs for each ancestral population
:fifth: parameter for linear transformation of PCs to % ancestry

.. rubric:: Data Fields (remaining rows)

:SNP rs number:
:reference allele:
:variant allele:
:reference allele frequency:
:SNP weight for PC1:
:SNP weight for PC2:


snpweights (``SNPweights``)
---------------------------

Per sample SNP weights based on an external reference panel.

.. rubric:: Fields

:Sample ID:
:Population Label:
:Number of SNPs: The number of SNPs used for inference of SNP weight.
:Predicted PC1:
:Predicted PC2:
:Percent YRI Ancestry:
:Percent CEU Ancestry:
:Percent ASI (CHB + CHD) Ancestry:


.. T

.. U

.. V

.. W

.. X

.. Y

.. Z
