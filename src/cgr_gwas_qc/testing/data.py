"""Testing data access.

The goal of this module is to make test data easily accessible in fixtures
and tests. Here we are hiding some of the complexity and provide a simple API
for copying test data into a working directory and creating a config.yml.

Currently, we have two sources for test data ``TestData`` and ``RealData``.

The ``TestData`` contains a set of very small test datasets that I collected
from the internet. These data are ideal for testing specific functionality
like converting file types and running upstream parts of the workflow.
However, the data are synthetic so they will not work with filtering steps
and various parts of the workflow.

The ``RealData`` contains ~200 real samples. These data are potentially
re-identifiable so they cannot be distributed outside of the network.
However, because they are real data they will behave well for testing
filtering steps and downstream parts of the workflow. This data also includes
all outputs for these samples from the Production workflow. This will allow
regression testing of most parts of the workflow.
"""
import os
import shutil
import subprocess
from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path
from typing import MutableMapping, Optional, TypeVar, Union
from warnings import warn

from cgr_gwas_qc.parsers.sample_sheet import SampleSheet
from cgr_gwas_qc.testing import make_snakefile, make_test_config

DEFAULT_TEST_DATA_SERVER = "cgemsiii.nci.nih.gov"
DEFAULT_TEST_DATA_PATH = "/DCEG/CGF/Bioinformatics/Production/fearjm/gwas_test_data"

# This is a trick to add a type annotation that says that the method is
# returning the class itself. This will magically update to the subclass too.
T = TypeVar("T", bound="DataRepo")
U = Union["DataRepo", T]  # I had to use the Union instead of `T` to get mypy to pass


class DataRepo(ABC):
    """Abstract Base Class for Test Data.

    My goal is to create a set of classes to easily copy test datasets to the
    test's working directory. The general functionality is very similar for
    different sources of tests so I opted to use class inheritance instead of
    duplicating code/logic.

    Similar to a regular class, an Abstract Base Class (ABC), allows you to
    create attributes/properties that are inheritable by subclasses. However,
    the ABC also allows you to define the attributes/methods that your
    subclasses need but not implement them.
    """

    # Required Class Attributes, make sure to define in your subclass.
    _data_type: str  # This will be the `project_name` in the config
    _data_path: Path  # Path to the data repository
    _sample_sheet: str  # name of the sample sheet relative to the data repository
    _illumina_manifest_file: str  # name of the BPM file
    _thousand_genome_vcf: str  # name of the VCF file
    _thousand_genome_tbi: str  # name of the VCF.TBI file
    _gtc_pattern: str  # GTC file name pattern
    _idat_red_pattern: str  # Red IDAT file name pattern
    _idat_green_pattern: str  # Green IDAT file name pattern
    _ped: str  # Name of the aggregate samples PED
    _map: str  # Name of the aggregate samples MAP
    _bed: str  # Name of the aggregate samples BED
    _bim: str  # Name of the aggregate samples BIM
    _fam: str  # Name of the aggregate samples FAM

    def __init__(self, working_dir: Optional[Path] = None):
        self.working_dir = working_dir
        self.ss = SampleSheet(self / self._sample_sheet)
        self._config: MutableMapping = defaultdict(dict)
        self._config["project_name"] = self._data_type

    @abstractmethod
    def _add_gtcs(self):
        """Each subclass must define this method.

        Handling of copying GTCs is very specific to the source of the test
        data. So I require that each subclass defines this method.
        """
        raise NotImplementedError

    @abstractmethod
    def _add_idats(self):
        """Each subclass must define this method.

        Handling of copying IDATs is very specific to the source of the test
        data. So I require that each subclass defines this method.
        """
        raise NotImplementedError

    def add_sample_sheet(self, copy: bool = True) -> U:
        """Add sample sheet.

        The sample sheet can either be copied into the working directory or
        have its full path referenced in the config. If a working directory
        is given (``self.working_dir``) then it will be copied otherwise it
        will be referenced in the config.
        """
        if self.working_dir is None or copy is False:
            # Don't copy. Add full path to sample sheet in the data repository.
            self._config["sample_sheet"] = (self / self._sample_sheet).absolute()
        else:
            self._config["sample_sheet"] = "sample_sheet.csv"
            self.copy(self._sample_sheet, self._config["sample_sheet"])

        return self

    def add_reference_files(self, manifest: bool = True, vcf: bool = True, copy: bool = True) -> U:
        """Add the various reference files.

        The reference files can either be copied into the working directory or
        have their full path referenced in the config. If a working directory is
        given (``self.working_dir``) then they will be copied otherwise they
        will be referenced in the config.

        Args:
            manifest: Include the Illumina manifest file (BPM). Defaults to True.
            vcf: Include the Thousand Genomes VCF/TBI files. Defaults to True.
        """
        if self.working_dir is None or copy is False:
            # Don't copy. Add full path to reference files in the data repository.
            self._config["reference_files"]["illumina_manifest_file"] = (
                self / self._illumina_manifest_file
            ).absolute()
            self._config["reference_files"]["thousand_genome_vcf"] = (
                self / self._thousand_genome_vcf
            ).absolute()
            self._config["reference_files"]["thousand_genome_tbi"] = (
                self / self._thousand_genome_tbi
            ).absolute()
        else:
            if manifest:
                self._config["reference_files"]["illumina_manifest_file"] = "manifest.bpm"
                self.copy(
                    self._illumina_manifest_file,
                    self._config["reference_files"]["illumina_manifest_file"],
                )

            if vcf:
                self._config["reference_files"]["thousand_genome_vcf"] = "thousG.vcf.gz"
                self.copy(
                    self._thousand_genome_vcf,
                    self._config["reference_files"]["thousand_genome_vcf"],
                )

                self._config["reference_files"]["thousand_genome_tbi"] = "thousG.vcf.gz.tbi"
                self.copy(
                    self._thousand_genome_tbi,
                    self._config["reference_files"]["thousand_genome_tbi"],
                )
        return self

    def add_user_files(self, entry_point: str = "bed", copy: bool = True) -> U:
        """Add user provided files.

        The user files can either be copied into the working directory or
        have their full path referenced in the config. If a working directory is
        given (``self.working_dir``) then they will be copied otherwise they
        will be referenced in the config.

        Args:
            entry_point: Which entry point to use ("gtc", "ped", or "bed").
              Defaults to "bed".

        Returns:
            T: [description]
        """
        if self.working_dir is None or copy is False:
            # Don't copy. Add full path to user files in the data repository.
            if entry_point == "gtc":
                self._config["user_files"]["gtc_pattern"] = (
                    self._data_path.absolute() / self._gtc_pattern
                ).as_posix()
                self._config["user_files"]["idat_pattern"] = {}
                self._config["user_files"]["idat_pattern"]["red"] = (
                    self._data_path.absolute() / self._idat_red_pattern
                ).as_posix()
                self._config["user_files"]["idat_pattern"]["green"] = (
                    self._data_path.absolute() / self._idat_green_pattern
                ).as_posix()

            elif entry_point == "ped":
                self._config["user_files"]["ped"] = (self / self._ped).absolute()
                self._config["user_files"]["map"] = (self / self._map).absolute()

            elif entry_point == "bed":
                self._config["user_files"]["bed"] = (self / self._bed).absolute()
                self._config["user_files"]["bim"] = (self / self._bim).absolute()
                self._config["user_files"]["fam"] = (self / self._fam).absolute()
        else:
            if entry_point == "gtc":
                self._config["user_files"]["gtc_pattern"] = self._gtc_pattern
                self._add_gtcs()

                self._config["user_files"]["idat_pattern"] = {}
                self._config["user_files"]["idat_pattern"]["red"] = self._idat_red_pattern
                self._config["user_files"]["idat_pattern"]["green"] = self._idat_green_pattern
                self._add_idats()

            elif entry_point == "ped":
                self._config["user_files"]["ped"] = "samples.ped"
                self.copy(self._ped, self._config["user_files"]["ped"])

                self._config["user_files"]["map"] = "samples.map"
                self.copy(self._map, self._config["user_files"]["map"])

            elif entry_point == "bed":
                self._config["user_files"]["bed"] = "samples.bed"
                self.copy(self._bed, self._config["user_files"]["bed"])

                self._config["user_files"]["bim"] = "samples.bim"
                self.copy(self._bim, self._config["user_files"]["bim"])

                self._config["user_files"]["fam"] = "samples.fam"
                self.copy(self._fam, self._config["user_files"]["fam"])

        return self

    def make_config(self, **kwargs) -> U:
        """Create a config.yml in ``self.working_dir``.

        The config.yml will be based on which objects were added using the
        ``self.add_*`` methods.

        Args:
            kwargs: Config options that you would like to set when creating
              the config.yml
        """
        self._config.update(kwargs)
        if self.working_dir:
            make_test_config(self.working_dir, **self._config)
        else:
            raise ValueError("You need to have set ``self.working_dir``.")
        return self

    def make_snakefile(self, contents: str) -> U:
        """Convenience wrapper of ``cgr_gwas_qc.testing.make_snakefile``."""
        if self.working_dir:
            make_snakefile(self.working_dir, contents)
        else:
            raise ValueError("You need to have set ``self.working_dir``.")
        return self

    def copy(self, source: str, destination: str) -> None:
        """Copy folder or file to a new directory

        Args:
            source: The name of a folder or file in the repository
              ``self._data_path`` that you want to copy to
              ``self.working_dir``.
            destination: The name where you want to copy the file or
              folder in ``self.working_dir. Default None.

        Raises:
            ValueError: If the ``source`` does not exist in the repository.
        """
        if self.working_dir:
            source_path = (self._data_path / source).absolute()
            destination_path = (self.working_dir / destination).absolute()
        else:
            raise ValueError("You need to have set ``self.working_dir``.")

        if source_path.is_dir():
            shutil.copytree(source_path, destination_path)
        elif source_path.is_file():
            shutil.copyfile(source_path, destination_path)
        else:
            raise ValueError(f"{source_path.as_posix()} is not a file or directory.")

    def __truediv__(self, value) -> Path:
        """Override division operator to build paths."""
        pth = self._data_path / value
        if pth.exists():
            return pth
        raise FileNotFoundError(pth.as_posix())


class FakeData(DataRepo):
    """A repository of synthetic data stored in ``tests/data``.

    Example:
        >>> cache = FakeData(working_dir)

        >>> cache / existing_file.txt  # If file exists then returns path to the file
        Path(".../existing_file.txt")

        >>> cache / non_existing_file.txt  # If file does not exist then raises exception
        FileNotFoundError: ".../non_existing_file.txt"

        >>> cache.add_sample_sheet()  # copy sample sheet to working dir
        >>> cache.add_reference_files()  # add reference files to working dir
        >>> cache.add_user_files()  # add BED/BIM/FAM to working dir
        >>> cache.make_config()  # make the config for files added using cache
    """

    _data_type = "Fake Data"
    _data_path = (Path(__file__).parents[3] / "tests/data").absolute()

    _sample_sheet = "example_sample_sheet.csv"

    _illumina_manifest_file = "illumina/bpm/small_manifest.bpm"
    _thousand_genome_vcf = "1KG/small_1KG.vcf.gz"
    _thousand_genome_tbi = "1KG/small_1KG.vcf.gz.tbi"

    _gtc_pattern = "{Sample_ID}.gtc"
    _idat_red_pattern = "{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat"
    _idat_green_pattern = "{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat"

    _ped = "plink/samples.ped"
    _map = "plink/samples.map"

    _bed = "plink/samples.bed"
    _bim = "plink/samples.bim"
    _fam = "plink/samples.fam"

    _test_gtc = "illumina/gtc/small_genotype.gtc"
    _test_idat = "illumina/idat/small_intensity.idat"

    def add_user_files(self, entry_point: str = "bed", copy: bool = True) -> U:
        if (self.working_dir is None or not copy) and entry_point == "gtc":
            # TestData only has 1 GTC file that I copy multiple times to
            # simulate samples. So this entry point must be copied into the
            # working directory.
            raise ValueError("Test GTC files must be copied to the working directory.")

        return super().add_user_files(entry_point, copy)

    def _add_gtcs(self):
        for _, r in self.ss.data.iterrows():
            self.copy(self._test_gtc, self._gtc_pattern.format(**dict(r)))

    def _add_idats(self):
        for _, r in self.ss.data.iterrows():
            self.copy(self._test_idat, self._idat_red_pattern.format(**dict(r)))
            self.copy(self._test_idat, self._idat_green_pattern.format(**dict(r)))


def _make_real_data_cache() -> Path:
    """Create cache folder."""
    suffix = "cgr_gwas_qc/test_data"

    if os.environ.get("XDG_CACHE_HOME"):
        _cache_path = Path(os.environ["XDG_CACHE_HOME"]) / suffix
    else:
        _cache_path = Path(__file__).absolute().parents[3] / ".cache" / suffix

    _cache_path.mkdir(parents=True, exist_ok=True)
    return _cache_path


class RealData(DataRepo):
    """A repository of real data stored on cgems.

    Example:
        >>> cache = ReadData(working_dir, sync=True)  # automatically symlink (on cgems) or rsync data locally

        >>> cache / existing_file.txt  # If file exists then returns path to the file
        Path("../../../.cache/cgr_gwas_qc/test_data/existing_file.txt")

        >>> cache / non_existing_file.txt  # If file does not exist then raises exception
        FileNotFoundError: "../../../.cache/cgr_gwas_qc/test_data/non_existing_file.txt"

        >>> cache.add_sample_sheet()  # copy sample sheet to working dir
        >>> cache.add_reference_files(copy=False)  # add full path of reference files to config
        >>> cache.add_user_files(copy=False)  # add full path of user files to config
        >>> cache.make_config()  # make the config for files added using cache
    """

    _data_type = "Real Data"
    _data_path = _make_real_data_cache()

    _sample_sheet = "original_data/manifest_short.csv"

    _illumina_manifest_file = "reference_data/GSAMD-24v1-0_20011747_A1.bpm"
    _thousand_genome_vcf = (
        "reference_data/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
    )
    _thousand_genome_tbi = (
        "reference_data/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi"
    )

    _gtc_pattern = "original_data/{SentrixBarcode_A}_{SentrixPosition_A}.gtc"
    _idat_red_pattern = "original_data/{SentrixBarcode_A}_{SentrixPosition_A}_Red.idat"
    _idat_green_pattern = "original_data/{SentrixBarcode_A}_{SentrixPosition_A}_Grn.idat"

    _ped = "production_outputs/plink_start/samples.ped"
    _map = "production_outputs/plink_start/samples.map"

    _bed = "production_outputs/plink_start/samples.bed"
    _bim = "production_outputs/plink_start/samples.bim"
    _fam = "production_outputs/plink_start/samples.fam"

    def __init__(self, working_dir: Optional[Path] = None, sync: bool = False):
        """A real test data repository of ~200 samples.

        These data are stored on cgems. If you are running this on a computer
        with access to cgems then providing ``sync=True`` will rsync the data
        into a local cache. If you are running this on cgems then
        ``sync=True`` will create a symlink into the cache folder.

        Args:
            sync: If True we will try to cache data locally. Defaults to False.

        """
        if sync:
            self._sync_data()

        super().__init__(working_dir)

    def add_sample_sheet(self, copy: bool = True, full_sample_sheet: bool = True) -> U:
        """Add sample sheet.

         We have two different sample sheets for the Real Data. The first
         `manifest_short.csv` that contains only the 2 samples that I have
         GTC/IDAT files for. The second `manifest_full.csv` contains the ~200
         samples. Essentially, you will only want the short sample sheet when
         using the `gtc` entry_point.

        Args:
            working_dir: If the path to the working directory is given then
              the sample sheet is copied into the working directory. Defaults to None.
            full_sample_sheet: Should you use the full or short sample sheet.
              Defaults to True.

        Returns:
            T: [description]
        """
        if full_sample_sheet:
            self._sample_sheet = "original_data/manifest_full.csv"

        return super().add_sample_sheet(copy=copy)

    def _add_gtcs(self):
        target_folder = self.working_dir / "original_data"
        target_folder.mkdir(exist_ok=True)
        for file_name in (self / "original_data").glob("*.gtc"):
            shutil.copy(file_name, target_folder)

    def _add_idats(self):
        target_folder = self.working_dir / "original_data"
        target_folder.mkdir(exist_ok=True)
        for file_name in (self / "original_data").glob("*.idat"):
            shutil.copy(file_name, target_folder)

    def _sync_data(self):
        """Symlink or Rsync data to cache folder."""
        user = os.environ.get("TEST_DATA_USER", None) or os.environ.get("USER")
        server = os.environ.get("TEST_DATA_SERVER", DEFAULT_TEST_DATA_SERVER)
        remote_path = Path(os.environ.get("TEST_DATA_PATH", DEFAULT_TEST_DATA_PATH))

        if remote_path.exists():
            self._data_path.symlink_to(remote_path)
            return

        # Rsync from remote server
        cmd = " ".join(
            [
                "rsync",
                "-a",
                f"{user}@{server}:{remote_path.as_posix()}/",
                f"{self._data_path.as_posix()}/",
            ]
        )
        try:
            subprocess.run(cmd, shell=True, check=True, stderr=subprocess.STDOUT)
        except subprocess.SubprocessError:
            warn(
                f"Could not connect to {user}@{server}:{remote_path}. To "
                "change set the environmental variable `TEST_DATA_USER`, "
                "`TEST_DATA_SERVER`, and `TEST_DATA_PATH."
            )
