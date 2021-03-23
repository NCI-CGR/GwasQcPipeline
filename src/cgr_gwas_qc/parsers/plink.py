from pathlib import Path

import pandas as pd


def read_het(filename: Path) -> pd.DataFrame:
    """Parse PLINK's het file format.

    Returns:
        pd.DataFrame:
            A (n x 5) table with the following columns

            .. csv-table::
                :header: name, dtype, description

                **ID** (*index*), object, Sample or Subject ID
                O(HOM), int, Observed number of homozygotes
                E(HOM), int, Expected number of homozygotes
                N(NM), int, Number of non-missing autosomal genotypes
                F, float, Method-of-moments F coefficient estimate

    References:
        - https://www.cog-genomics.org/plink/1.9/basic_stats#ibc
        - https://www.cog-genomics.org/plink/1.9/formats#het
    """
    return (
        pd.read_csv(filename, delim_whitespace=True)
        .drop("FID", axis=1)
        .rename({"IID": "ID"}, axis=1)
        .set_index("ID")
    )


def read_imiss(filename: Path) -> pd.DataFrame:
    """Parse PLINK's imiss file format.

    This reports sample level missing data.

    Returns:
        pd.DataFrame:
            A (n x 5) table with the following columns

            .. csv-table::
                :header: name, dtype, description

                **ID** (*index*), object, Sample or Subject ID
                MISS_PHENO, str, [Y/N] if the phenotype is missing
                N_MISS, int, The number of missing genotype calls not including obligatory missing or heterozygous haploids.
                N_GENO, int, Number of potentially valid calls
                F_MISS, float, Missing call rate

    References:
        - https://www.cog-genomics.org/plink/1.9/basic_stats#missing
        - https://www.cog-genomics.org/plink/1.9/formats#imiss
    """
    return (
        pd.read_csv(filename, delim_whitespace=True)
        .drop("FID", axis=1)
        .rename({"IID": "ID"}, axis=1)
        .set_index("ID")
    )


def read_lmiss(filename: Path) -> pd.DataFrame:
    """Parse PLINK's lmiss file format.

    This reports snp level missing data.

    Returns:
        pd.DataFrame:
            A (n x 5) table with the following columns

            .. csv-table::
                :header: name, dtype, description

                CHR, str, Chromosome code
                SNP, str, Variant identifier
                N_MISS, int, The number of missing genotype calls not including obligatory missing or heterozygous haploids.
                N_GENO, int, Number of potentially valid calls
                F_MISS, float, Missing call rate

    References:
        - https://www.cog-genomics.org/plink/1.9/basic_stats#missing
        - https://www.cog-genomics.org/plink/1.9/formats#lmiss
    """
    return pd.read_csv(filename, delim_whitespace=True)


def read_sexcheck(filename: Path) -> pd.DataFrame:
    """Parse PLINK's sexcheck file format.

    Returns:
        pd.DataFrame:
            A (n x 5) table with the following columns

            .. csv-table::
                :header: name, dtype, description

                **ID** (*index*), object, Sample or Subject ID
                PEDSEX, int, Sex code in input file
                SNPSEX, int, Imputed sex code (1 = male, 2 = female, 0 = unknown)
                STATUS, str, "OK" if PEDSEX and SNPSEX match and are nonzero, "PROBLEM" otherwise
                F, int, Inbreeding coefficient, considering only X chromosome.

    References:
        - https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex
        - https://www.cog-genomics.org/plink/1.9/formats#sexcheck
    """
    return (
        pd.read_csv(filename, delim_whitespace=True)
        .drop("FID", axis=1)
        .rename({"IID": "ID"}, axis=1)
        .set_index("ID")
    )
