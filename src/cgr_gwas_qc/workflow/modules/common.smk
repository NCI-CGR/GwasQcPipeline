from cgr_gwas_qc.parsers.illumina import BeadPoolManifest


numSNPs = BeadPoolManifest(cfg.config.reference_files.illumina_manifest_file).num_loci
