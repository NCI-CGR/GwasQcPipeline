from pydantic import BaseModel


class EnvModules(BaseModel):
    """Environmental modules found on CGEMS or other HPC Systems."""

    eigensoft: str
    graf: str
    plink2: str
    r: str
