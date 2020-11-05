from pydantic import BaseModel


class EnvModules(BaseModel):
    """Environmental modules found on CGEMS or other HPC Systems."""

    plink2: str
    eigensoft: str
    r: str
