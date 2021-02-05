from typing import List

import pytest

from cgr_gwas_qc.parsers import bim, vcf
from cgr_gwas_qc.workflow.scripts.update_snps_to_1kg_rsID import extract_rsID, update_record_id

records_for_update = [
    pytest.param(
        bim.BimRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs1kg", "1", 100, "A", ("G",))],
        "rs1kg",
        id="Update-Exact_Match",
    ),
    pytest.param(
        bim.BimRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs1kg", "1", 100, "G", ("A",))],
        "rs1kg",
        id="Update-Alleles_Flipped",
    ),
    pytest.param(
        bim.BimRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs1kg", "1", 100, "T", ("C",))],
        "rs1kg",
        id="Update-Complement",
    ),
    pytest.param(
        bim.BimRecord("GSA-rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs1kg", "1", 999, "A", ("G",))],
        "rs123",
        id="Update-Extracted_RSID",  # NOTE: update is based on the bim record not the vcf which doesn't match
    ),
    pytest.param(
        bim.BimRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs1kg", "1", 999, "A", ("G",))],
        "rs123",
        id="NoUpdate-Different_Position",
    ),
    pytest.param(
        bim.BimRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs1kg", "1", 100, "A", ("G", "C"))],
        "rs123",
        id="NoUpdate-Multi_Allelic",
    ),
    pytest.param(
        bim.BimRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("rs1kg", "1", 100, "A", ("GTT",))],
        "rs123",
        id="NoUpdate-Not_Snp",
    ),
    pytest.param(
        bim.BimRecord("rs123", "1", 100, "A", "G"),
        [vcf.VcfRecord("kg123", "1", 100, "A", ("G",))],
        "rs123",
        id="NoUpdate-Not_RsID",
    ),
]


@pytest.mark.parametrize("bim_record,vcf_records,updated_id", records_for_update)
def test_update_record_id(
    bim_record: bim.BimRecord, vcf_records: List[vcf.VcfRecord], updated_id: str, vcf_mock
):
    update_record_id(bim_record, vcf_mock(vcf_records))
    assert updated_id == bim_record.id


@pytest.mark.parametrize(
    "chip_id,rsID",
    [
        pytest.param("rs123", "rs123", id="A simple rsID"),
        pytest.param("GSA-rs123", "rs123", id="GSA- prefix"),
        pytest.param("GSArs123", "rs123", id="GSA prefix"),
        pytest.param("rs123GSA", "rs123", id="GSA suffix"),
        pytest.param("GSArs123GSA", "rs123", id="GSA prefix and suffix"),
        pytest.param("ID123", "ID123", id="No rsID gives ID untouched"),
    ],
)
def test_extract_rsID(chip_id, rsID):
    assert rsID == extract_rsID(chip_id)
