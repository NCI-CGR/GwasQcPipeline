from cgr_gwas_qc.workflow.scripts import trim_ped_map_ids


def test_trim_ped_ids(fake_data_cache):
    filename = fake_data_cache / "plink/samples.ped"
    trimmed = trim_ped_map_ids.trim_ped_ids(filename, 2)
    row = next(trimmed)
    columns = row.split(" ")
    assert 2 == len(columns[0])
    assert 2 == len(columns[1])


def test_trim_map_ids(fake_data_cache):
    filename = fake_data_cache / "plink/samples.map"
    trimmed = trim_ped_map_ids.trim_map_ids(filename, 2)
    row = next(trimmed)
    columns = row.split("\t")
    assert 2 == len(columns[1])
