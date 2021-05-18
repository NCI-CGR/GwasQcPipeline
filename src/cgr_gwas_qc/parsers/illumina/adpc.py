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
        """A record from Illumina's adpc.bin file.

        Args:
            x_raw: A allele intensity
            y_raw: B allele intensity
            x_norm: A allele normalized intensity
            y_norm: B allele normalize intensity
            genotype_score: Genotype call score (a.k.a. clustering confidence)
            genotype: Genotype call (0: AA, 1: AB, 2: BB, 3: unknown/missing)
        """
        # Coerce to specific C-types b/c I am dealing with reading/writing binary data.
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
            (
                self.x_raw,
                self.y_raw,
                self.x_norm,
                self.y_norm,
                self.genotype_score,
                self.genotype,
            )
        )

    def __str__(self):
        return "\t".join(
            (
                str(self.x_raw),
                str(self.y_raw),
                str(self.x_norm),
                str(self.y_norm),
                str(self.genotype_score),
                str(self.genotype),
            )
        )


class AdpcBase:
    """Base class for Illumina's adpc.bin files.

    Based on information found on Picard's website:

    https://javadoc.io/static/org.broadinstitute/gatk/4.1.4.1/picard/arrays/illumina/IlluminaAdpcFileWriter.html
    """

    padding_struct = struct.Struct("<8H")  # There is a 16 byte padding
    record_struct = struct.Struct("<HHfffH")  # Followed by 18 byte records


class AdpcReader(AdpcBase):
    def __init__(self, file_name: Union[str, Path]):
        """A context manager for reading adpc.bin formatted files.

        Args:
            file_name: Path to an adpc.bin file.

        Example:
            >>> with AdpcReader(file_name) as fh:
            ...     for record in fh:
            ...         print(record)
        """
        self.file_name = Path(file_name)

    def open(self):
        self.file_handler = self.file_name.open("rb")

    def close(self):
        self.file_handler.close()

    def read(self, bytes):
        return self.file_handler.read(bytes)

    def __enter__(self):
        self.open()
        self.read(16)  # skip padding
        return self

    def __exit__(self, *exc):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        data_bytes = self.read(18)
        if data_bytes:
            return AdpcRecord(*self.record_struct.unpack(data_bytes))
        else:
            raise StopIteration


class AdpcWriter(AdpcBase):
    def __init__(self, file_name: Union[str, Path]):
        """A context manager for writing adpc.bin formatted files.

        Args:
            file_name: Path to save an adpc file.

        Example:
            >>> with AdpcReader(file_name) as fh:
            ...    fh.write(record)
        """
        self.file_name = Path(file_name)

    def open(self):
        self.file_handler = self.file_name.open("wb")

    def close(self):
        self.file_handler.close()

    def write(self, record: AdpcRecord):
        self.file_handler.write(self.record_struct.pack(*record))

    def _add_padding(self):
        """The adpc.bin format has 16 bytes of padding at the beginning of the file.

        The Picard group did not know what this is for, but they said it is
        not needed so I will just add 8 0s at 2 bytes each.
        """
        self.file_handler.write(self.padding_struct.pack(*[0] * 8))

    def __enter__(self):
        self.open()
        self._add_padding()
        return self

    def __exit__(self, *exc):
        self.close()
