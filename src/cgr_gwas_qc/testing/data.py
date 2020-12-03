"""Cache real data for testing.

We have a set of real data for use in regression testing. Currently, we are
keeping these data private and will only be used for testing on NIH machines.
"""
import os
import shutil
import subprocess
from collections import defaultdict
from logging import getLogger
from pathlib import Path
from typing import Optional, Union

from cgr_gwas_qc.parsers.sample_sheet import SampleSheet
from cgr_gwas_qc.testing import make_test_config

logger = getLogger(__name__)


class TestData:

    test_data_path = (Path(__file__).parents[3] / "tests/data").absolute()
    gtc_pattern = "{Sample_ID}.gtc"
    idat_red_pattern = "{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat"
    idat_green_pattern = "{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat"

    def __init__(self):
        self.ss = SampleSheet(self.test_data_path / "example_sample_sheet.csv")

        self._config = defaultdict(dict)
        self._config["project_name"] = "Test Data"

    def copy(self, source: str, destination: Path):
        if (self.test_data_path / source).is_dir():
            shutil.copytree(self.test_data_path / source, destination.absolute())
        elif (self.test_data_path / source).is_file():
            shutil.copyfile(self.test_data_path / source, destination.absolute())
        else:
            raise FileNotFoundError(f"{(self.test_data_path / source).as_posix()} does not exist.")

    def copy_sample_sheet(self, working_dir: Path) -> "TestData":
        self._config["sample_sheet"] = "sample_sheet.csv"
        self.copy(
            "example_sample_sheet.csv", working_dir / self._config["sample_sheet"],
        )

        return self

    def copy_reference_files(
        self, working_dir: Path, sample_sheet: bool = True, manifest: bool = True, vcf: bool = True
    ) -> "TestData":
        if sample_sheet:
            self._config["sample_sheet"] = "sample_sheet.csv"
            self.copy(
                "example_sample_sheet.csv", working_dir / self._config["sample_sheet"],
            )

        if manifest:
            self._config["reference_files"]["illumina_manifest_file"] = "manifest.bpm"
            self.copy(
                "illumina/bpm/small_manifest.bpm",
                working_dir / self._config["reference_files"]["illumina_manifest_file"],
            )

        if vcf:
            self._config["reference_files"]["thousand_genome_vcf"] = "thousG.vcf.gz"
            self.copy(
                "1KG/small_1KG.vcf.gz",
                working_dir / self._config["reference_files"]["thousand_genome_vcf"],
            )

            self._config["reference_files"]["thousand_genome_tbi"] = "thousG.vcf.gz.tbi"
            self.copy(
                "1KG/small_1KG.vcf.gz.tbi",
                working_dir / self._config["reference_files"]["thousand_genome_tbi"],
            )

        return self

    def copy_user_files(self, working_dir: Path, entry_point: str = "bed") -> "TestData":
        if entry_point == "gtc":
            self._create_gtcs(working_dir)
            self._create_idats(working_dir)

        elif entry_point == "ped":
            self._config["user_files"]["ped"] = "samples.ped"
            self.copy(
                "plink/samples.ped", working_dir / self._config["user_files"]["ped"],
            )

            self._config["user_files"]["map"] = "samples.map"
            self.copy(
                "plink/samples.map", working_dir / self._config["user_files"]["map"],
            )

        elif entry_point == "bed":
            self._config["user_files"]["bed"] = "samples.bed"
            self.copy(
                "plink/samples.bed", working_dir / self._config["user_files"]["bed"],
            )

            self._config["user_files"]["bim"] = "samples.bim"
            self.copy(
                "plink/samples.bim", working_dir / self._config["user_files"]["bim"],
            )

            self._config["user_files"]["fam"] = "samples.fam"
            self.copy(
                "plink/samples.fam", working_dir / self._config["user_files"]["fam"],
            )

        return self

    def make_config(self, working_dir, **kwargs) -> "TestData":
        self._config.update(kwargs)
        make_test_config(working_dir, **self._config)
        return self

    def _create_gtcs(self, working_dir: Path):
        self._config["user_files"]["gtc_pattern"] = self.gtc_pattern

        for _, r in self.ss.data.iterrows():
            self.copy(
                "illumina/gtc/small_genotype.gtc", working_dir / self.gtc_pattern.format(**dict(r)),
            )

    def _create_idats(self, working_dir: Path):
        self._config["user_files"]["idat_pattern"] = {}
        self._config["user_files"]["idat_pattern"]["red"] = self.idat_red_pattern
        self._config["user_files"]["idat_pattern"]["green"] = self.idat_green_pattern

        for _, r in self.ss.data.iterrows():
            self.copy(
                "illumina/idat/small_intensity.idat",
                working_dir / self.idat_red_pattern.format(**dict(r)),
            )

            self.copy(
                "illumina/idat/small_intensity.idat",
                working_dir / self.idat_green_pattern.format(**dict(r)),
            )


class RealDataCache:

    user = os.environ.get("TEST_DATA_USER", None) or os.environ.get("USER")
    server = os.environ.get("TEST_DATA_SERVER", "cgemsiii.nci.nih.gov")
    test_data_path = Path(
        os.environ.get(
            "TEST_DATA_PATH", "/DCEG/CGF/Bioinformatics/Production/fearjm/gwas_test_data"
        )
    )

    def __init__(self, cache_path: Optional[Union[str, Path]] = None):
        """Cache real test data.

        If test data path exists then create a symbolic link to test data. If
        the test data path does not exist then try to rsync data to cache.

        Args:
            cache_path: Path to use for cache. Defaults to
              ``GwasQcPipeline/.cache``.

        Example:
            >>> cache = ReadDataCache()  # automatically symlink or rsync data

            >>> cache / existing_file.txt  # If file exists then returns path to the file
            Path("../../../.cache/cgr_gwas_qc/test_data/existing_file.txt")

            >>> cache / non_existing_file.txt  # If file does not exist then raises exception
            FileNotFoundError: "../../../.cache/cgr_gwas_qc/test_data/non_existing_file.txt"
        """
        self._cache_path = self._make_cache(cache_path)
        self._sync_data()

    def _make_cache(self, user_path: Optional[Union[str, Path]]) -> Path:
        """Create cache folder."""
        suffix = "cgr_gwas_qc/test_data"

        if user_path:
            _cache_path = Path(user_path) / suffix
        elif os.environ.get("XDG_CACHE_HOME"):
            _cache_path = Path(os.environ["XDG_CACHE_HOME"]) / suffix
        else:
            _cache_path = Path(__file__).absolute().parents[3] / ".cache" / suffix

        _cache_path.mkdir(parents=True, exist_ok=True)
        return _cache_path

    def _sync_data(self):
        """Symlink or Rsync data to cache folder."""
        if self.test_data_path.exists():
            self._cache_path.symlink_to(self.test_data_path)
            return

        # Rsync from remote server
        cmd = " ".join(
            [
                "rsync",
                "-a",
                f"{self.user}@{self.server}:{self.test_data_path.as_posix()}/",
                f"{self._cache_path.as_posix()}/",
            ]
        )
        try:
            subprocess.run(cmd, shell=True, check=True, stderr=subprocess.STDOUT)
        except subprocess.SubprocessError:
            logger.warn(f"Could not connect to {self.server}.")

    def copy(self, source: str, destination: Path):
        if (self / source).is_dir():
            shutil.copytree(self / source, destination.absolute())
        elif (self / source).is_file():
            shutil.copyfile(self / source, destination.absolute())
        else:
            raise FileNotFoundError(f"{(self / source).as_posix()} does not exist.")

    def __truediv__(self, value):
        """Override division operator to build paths."""
        pth = self._cache_path / value
        if pth.exists():
            return pth
        raise FileNotFoundError(pth.as_posix())
