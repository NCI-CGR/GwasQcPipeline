"""Test parsing of Illumina BPM files."""
import pytest

from cgr_gwas_qc.parsers.illumina import BeadPoolManifest


@pytest.fixture(scope="module")
def bpm(bpm_file) -> BeadPoolManifest:
    return BeadPoolManifest(bpm_file)


bpm_list_attributes = [
    "names",
    "snps",
    "chroms",
    "map_infos",
    "addresses",
    "normalization_ids",
    "assay_types",
    "normalization_lookups",
    "ref_strands",
    "source_strands",
]


@pytest.mark.parametrize("attribute", bpm_list_attributes)
def test_bpm_attribute_lengths(bpm, attribute):
    """Check that the list element attributes all have the same length."""
    assert len(bpm.__dict__[attribute]) == bpm.num_loci
