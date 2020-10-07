import struct
from pathlib import Path
from typing import Union

from numpy import float32, uint16


class AdpcRecord:
    def __init__(
        self,
        x_raw: int,
        y_raw: int,
        x_norm: float,
        y_norm: float,
        genotype_score: float,
        genotype: int,
    ):
        """A record from an Adpc file.

        Args:
            x_raw: A allele intensity
            y_raw: B allele intensity
            x_norm: A allele normalized intensity
            y_norm: B allele normalize intensity
            genotype_score: Genotype call score (a.k.a. clustering confidence)
            genotype: Genotype call (0: AA, 1: AB, 2: BB, 3: unknown/missing)
        """
        # Coerce to specific types b/c I am dealing with reading/writing binary data.
        self.x_raw = uint16(x_raw)
        self.y_raw = uint16(y_raw)
        self.x_norm = float32(x_norm)
        self.y_norm = float32(y_norm)
        self.genotype_score = float32(genotype_score)
        self.genotype = uint16(genotype)

        self.validate()

    def validate(self):
        """Simple validation, not sure if needed."""
        if self.x_raw < 0:
            raise ValueError(f"Raw X intensity must be positive integer: {self.x_raw}")

        if self.y_raw < 0:
            raise ValueError(f"Raw Y intensity must be positive integer: {self.y_raw}")

        if self.x_norm < 0:
            raise ValueError(f"Norm X intensity must be positive number: {self.x_norm}")

        if self.y_norm < 0:
            raise ValueError(f"Norm Y intensity must be positive number: {self.y_norm}")

        if self.genotype_score < 0:
            raise ValueError(f"Genotype score must be positive number: {self.genotype_score}")

        if self.genotype not in [0, 1, 2, 3]:
            raise ValueError(f"Genotype must be [0, 1, 2, 3]: {self.genotype}")

    def __iter__(self):
        """This allows tuple unpacking"""
        return iter(
            (self.x_raw, self.y_raw, self.x_norm, self.y_norm, self.genotype_score, self.genotype,)
        )

    def __str__(self):
        return "\t".join(
            (self.x_raw, self.y_raw, self.x_norm, self.y_norm, self.genotype_score, self.genotype,)
        )


class AdpcBase:
    """Base class for Adpc files.

    See https://javadoc.io/static/org.broadinstitute/gatk/4.1.4.1/picard/arrays/illumina/IlluminaAdpcFileWriter.html
    """

    padding_struct = struct.Struct("<8H")  # There is a 16 byte padding
    record_struct = struct.Struct("<HHfffH")  # Followed by 18 byte records


class AdpcReader(AdpcBase):
    def __init__(self, file_name: Union[str, Path]):
        """Adpc reader context manager.

        Args:
            file_name: Path to an adpc file.

        Example:
            >>> file_name = "../../../tests/data/small.adpc.bin"
            >>> with AdpcReader(file_name) as fh:
            ...     for record in fh:
            ...         print(record)
        """
        self.file_name = Path(file_name)

    def __enter__(self):
        self.file_handler = self.file_name.open("rb")
        self.file_handler.read(16)  # skip padding
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        self.file_handler.close()

    def __iter__(self):
        return self

    def __next__(self):
        data_bytes = self.file_handler.read(18)
        if data_bytes:
            return AdpcRecord(*self.record_struct.unpack(data_bytes))
        else:
            raise StopIteration


class AdpcWriter(AdpcBase):
    def __init__(self, file_name: Union[str, Path]):
        """Adpc writer context manager.

        Args:
            file_name: Path to save an adpc file.

        Example:
            >>> file_name = "/tmp/small.adpc.bin"
            >>> with AdpcReader(file_name) as fh:
            ...    fh.write(record)
        """
        self.file_name = Path(file_name)

    def __enter__(self):
        self.file_handler = self.file_name.open("wb")
        self.file_handler.write(self.padding_struct.pack(*[0] * 8))  # add padding
        return self

    def __exit__(self, *args):
        self.file_handler.close()

    def close(self):
        self.file_handler.close()

    def write(self, record: AdpcRecord):
        self.file_handler.write(self.record_struct.pack(*record))
