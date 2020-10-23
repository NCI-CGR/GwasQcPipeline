from pydantic import BaseModel


class EnvModules(BaseModel):
    """Environmental modules found on CGEMS or other HPC Systems."""

    plink2: str = "plink2/1.90b5"
    eigensoft: str = "eigensoft/7.2.1"
    r: str = "R/3.4.0"
