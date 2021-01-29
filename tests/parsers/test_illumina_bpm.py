"""Test parsing of Illumina BPM files."""
import pytest

################################################################################
# Most BPM attrributes have an entry for each loci. Check that these entries
# are equal to the number of loci.
################################################################################
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
    # GIVEN: attributes of the Illumina BeadPoolManifest
    # WHEN: these attributes are a list (i.e. a column of values from the file)
    # THEN: the number of values is equal to the number of SNPs in the file (i.e., number of rows)
    assert len(bpm.__dict__[attribute]) == bpm.num_loci
