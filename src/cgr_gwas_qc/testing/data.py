"""Cache real data for testing.

We have a set of real data for use in regression testing. Currently, we are
keeping these data private and will only be used for testing on NIH machines.
"""
import os
import shutil
import subprocess
from logging import getLogger
from pathlib import Path
from typing import Optional, Union

USER = os.environ.get("TEST_DATA_USER", None) or os.environ.get("USER")
SERVER = os.environ.get("TEST_DATA_SERVER", "cgemsiii.nci.nih.gov")
TEST_DATA_PATH = os.environ.get(
    "TEST_DATA_PATH", "/DCEG/CGF/Bioinformatics/Production/fearjm/gwas_test_data"
)

logger = getLogger(__name__)


class RealDataCache:

    test_data_path = Path(TEST_DATA_PATH)

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
                f"{USER}@{SERVER}:{self.test_data_path.as_posix()}/",
                f"{self._cache_path.as_posix()}/",
            ]
        )
        try:
            subprocess.run(cmd, shell=True, check=True, stderr=subprocess.STDOUT)
        except subprocess.SubprocessError:
            logger.warn(f"Could not connect to {SERVER}.")

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
