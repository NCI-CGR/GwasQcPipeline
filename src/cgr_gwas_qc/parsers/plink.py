import numpy as np
import pandas as pd

from cgr_gwas_qc.typing import PathLike


def read_het(filename: PathLike) -> pd.DataFrame:
    """Parse PLINK's het file format.

    Returns:
        pd.DataFrame:
            A (n x 5) table with the following columns

            .. csv-table::
                :header: name, dtype, description

                **ID** (*index*), object, Sample or Subject ID
                O_HOM, int, Observed number of homozygotes
                E_HOM, int, Expected number of homozygotes
                N_NM, int, Number of non-missing autosomal genotypes
                F, float, Method-of-moments F coefficient estimate

    References:
        - https://www.cog-genomics.org/plink/1.9/basic_stats#ibc
        - https://www.cog-genomics.org/plink/1.9/formats#het
    """
    return (
        pd.read_csv(filename, delim_whitespace=True, dtype={"FID": "string", "IID": "string"})
        .drop("FID", axis=1)
        .rename({"IID": "ID", "O(HOM)": "O_HOM", "E(HOM)": "E_HOM", "N(NM)": "N_NM"}, axis=1)
        .set_index("ID")
    )


def read_hwe(filename: PathLike) -> pd.DataFrame:
    """Parse PLINK's hwe file format.

    Returns:
        pd.DataFrame:
            A (# SNP x 9) table with the following columns

            .. csv-table::
                :header: name, dtype, description

                **SNP** (*index*), string, SNP ID
                CHR, string, Chromosome code
                TEST, string, Type of test; one of {ALL, AFF, UNAFF, ALL(QT), ALL(NP)}
                A1, string, Allele 1 (usually minor)
                A2, string, Allele 2 (usually major)
                GENO, string, '/'- separated genotype counts (A1 hom, het, A2 hom)
                O_HET, float, Observed heterozygote frequency
                E_HET, float, Expected heterozygote frequency
                P, float, Hardy-Weinberg equilibrium exact test p-value

    References:
        - https://www.cog-genomics.org/plink/1.9/basic_stats#hardy
        - https://www.cog-genomics.org/plink/1.9/formats#hwe
    """
    dtypes = {
        "CHR": "string",
        "SNP": "string",
        "TEST": "string",
        "A1": "string",
        "A2": "string",
        "GENO": "string",
        "O(HET)": "float",
        "E(HET)": "float",
        "P": "float",
    }
    return (
        pd.read_csv(filename, delim_whitespace=True, dtype=dtypes)
        .rename({"O(HET)": "O_HET", "E(HET)": "E_HET"}, axis=1)
        .set_index("SNP")
    )


def read_genome(filename: PathLike, required_cols=None, chunksize=1000000) -> pd.DataFrame:
    """Parse PLINK's genome file format.

    Each row of the genome file is a pairwise combinations of
    samples/subjects.

    Returns:
        pd.DataFrame:
            A (n x 17) table with the following columns

            .. csv-table::
                :header: name, dtype, description

                ID1, string, First Sample or Subject ID (alphanumerically)
                ID2, string, Second Sample or Subject ID (alphanumerically)
                RT, category, Relationship type inferred from .fam/.ped file {FS: Full Sib, HS, Half Sib, PO: Parent-Offspring, OT; Other}
                EZ, object, IBD sharing expected value, based on just .fam/.ped relationship
                Z0, float, P(IBD=0)
                Z1, float, P(IBD=1)
                Z2, float, P(IBD=2)
                PI_HAT, float, Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)
                PHE, category, Pairwise phenotypic code (1, 0, -1 = case-case, case-ctrl, and ctrl-ctrl pairs, respectively)
                DST, float, IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)
                PPC, float, IBS binomial test
                RATIO, float, HETHET: IBS0 SNP ratio (expected value 2)
                IBS0, int, Number of IBS 0 nonmissing variants
                IBS1, int, Number of IBS 1 nonmissing variants
                IBS2, int, Number of IBS 2 nonmissing variants
                HOMHOM, int, Number of IBS 0 SNP pairs used in PPC test
                HETHET, int, Number of IBS 2 het/het SNP pairs used in PPC test

    References:
        - https://www.cog-genomics.org/plink/1.9/ibd
        - https://www.cog-genomics.org/plink/1.9/formats#genome
    """

    DTYPES = {
        "IID1": "string",
        "IID2": "string",
        "RT": "category",
        "EZ": "category",
        "Z0": np.float32,
        "Z1": np.float32,
        "Z2": np.float32,
        "PI_HAT": np.float32,
        "PHE": "category",
        "DST": np.float32,
        "PPC": np.float32,
        "RATIO": np.float32,
        "IBS0": np.uint32,
        "IBS1": np.uint32,
        "IBS2": np.uint32,
        "HOMHOM": np.uint32,
        "HETHET": np.uint32,
    }

    if required_cols is None:

        def required_cols(col_name):
            return lambda col_name: col_name not in ["FID1", "FID2"]

    return pd.read_csv(
        filename, delim_whitespace=True, dtype=DTYPES, usecols=required_cols, chunksize=chunksize
    )


def read_imiss(filename: PathLike) -> pd.DataFrame:
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
        pd.read_csv(filename, delim_whitespace=True, dtype={"FID": "string", "IID": "string"})
        .drop("FID", axis=1)
        .rename({"IID": "ID"}, axis=1)
        .set_index("ID")
    )


def read_lmiss(filename: PathLike) -> pd.DataFrame:
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
    dtypes = {
        "CHR": "string",
        "SNP": "string",
        "N_MISS": "UInt32",
        "N_GENO": "UInt32",
        "F_MISS": "float",
    }
    return pd.read_csv(filename, delim_whitespace=True, dtype=dtypes)


def read_sexcheck(filename: PathLike) -> pd.DataFrame:
    """Parse PLINK's sexcheck file format.

    Returns:
        pd.DataFrame:
            A (n x 5) table with the following columns

            .. csv-table::
                :header: name, dtype, description

                **ID** (*index*), object, Sample or Subject ID
                PEDSEX, int, Sex code in input file
                SNPSEX, int, Imputed sex code (1 = male; 2 = female; 0 = unknown)
                STATUS, str, OK if PEDSEX and SNPSEX match and are nonzero PROBLEM otherwise
                F, int, Inbreeding coefficient considering only X chromosome.

    References:
        - https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex
        - https://www.cog-genomics.org/plink/1.9/formats#sexcheck
    """
    return (
        pd.read_csv(filename, delim_whitespace=True, dtype={"FID": "string", "IID": "string"})
        .drop("FID", axis=1)
        .rename({"IID": "ID"}, axis=1)
        .set_index("ID")
    )
