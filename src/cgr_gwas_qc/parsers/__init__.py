import builtins
from abc import ABC, abstractmethod
from typing import Optional


class CgrFile(ABC):
    endchar: Optional[str] = None

    def __init__(self, filename, mode="r"):

        if hasattr(filename, "readlines"):
            self.fileobj = filename
        else:
            self.fileobj = builtins.open(filename, mode)

    def close(self):
        self.fileobj.close()

    def write(self, record, endchar=None):
        payload = str(record)
        _endchar = endchar or self.endchar

        if _endchar is not None:
            payload += _endchar

        self.fileobj.write(payload)

    @abstractmethod
    def __iter__(self):
        pass
