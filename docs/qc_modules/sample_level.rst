Sample QC
=========

Genotype Call Rate Filters
--------------------------

This module uses ``plink`` to look at Sample and SNP level genotype call rates. Two call rate thresholds are defined in the config (80%/95%). Low call rate is an indicator of poor sample quality or sample contamination resulting in a large amount of missing data. Quality control metrics (MAF, pct missing, HWE, autosomal homozygosity, X-linked homozygosity) are calculated for (starting data, filtered data at 80% thresholds, and filtered data at 95% threshold) using ``plink``.

.. mermaid::

    graph TD
        D1{Genotype Calls}
        D2{Illumina Manifest}
        S1[gtc_to_ped]
        S2[merge_sample_peds]
        S3[plink_sample_qc_stats::start]
        S4[filter_call_rate_1]
        S5[plink_sample_qc_stats::filter1]
        S6[filter_call_rate_2]
        S7[plink_sample_qc_stats::filter2]
        D1 -->|"{sample_id}.gtc"| S1
        D2 -->|bpm| S1
        S1 -->|"{sample_id}.ped + {sample_id}.map"| S2
        S2 -->|bed + bim + fam| S3 & S4
        S4 -->|bed + bim + fam| S5 & S6
        S6 -->|bed + bim + fam| S7

Sample Contamination Detection
------------------------------

Mixing of 2+ samples or swapping sample labels is a common problem for high-throughput genomics. Contamination is detected by modeling allele intensities accounting for population level allele frequencies. Here we use the tool ``verifyIDintensity`` to perform this analysis on individual
samples.

Sample swaps between male and female samples can be detected by looking at X-linked SNP heterozygosity.

.. mermaid::

    graph TD
        D1{Genotype Calls}
        D2{Illumina Manifest}
        D3{1000 Genome<br>VCF}

        A[gtc_to_adpc]
        D1 -->|<sample_id>.gtc| A
        D2 -->|bmp| A

        B[make_abf_txt]
        D2 -->|bmp| B
        D3 -->|vcf.gz| B

        C[contam_test_GSA_1000g]
        A  -->|<sample_id>.adpc.bin| C
        B -->|abf.txt| C

        D[combine_contam_one_samp_b_1000g]
        C -->|<sample_id>.contam.out| D

        %% I don't know where these are from yet
        intens(all_sample_idat_intensity/idat_intensity.csv)
        imiss3(plink_filter_call_rate_2/samples_filter2.imiss)
        intens --> D
        imiss3 --> D
