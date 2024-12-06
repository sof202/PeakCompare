import pandas as pd
from typing import Optional, TypeVar
pandasDataFrame = TypeVar('pandas.core.frame.DataFrame')


def read_bedgraph(file_path: str) -> Optional[pandasDataFrame]:
    """Reads a bedgraph file into a pandas DataFrame.

    Args:
        file_path (str): The path to the bedgraph file.

    Returns:
        Optional[pandasDataFrame]: The DataFrame containing the bedgraph data,
        or None if an error occurred.
    """
    try:
        bedgraph = pd.read_table(file_path, sep="\t")
        if bedgraph.shape[1] != 4:
            raise ValueError(f"Bedgraph file at {file_path}"
                             "does not have exactly 4 columns.")
        bedgraph.columns = ["CHR", "START", "END", "SCORE"]
        return bedgraph
    except (FileNotFoundError, IOError):
        print(f"{file_path} does not exist or could not be read.")
        return None
    except IsADirectoryError:
        print(f"{file_path} is a directory.")
        return None
    except PermissionError:
        print(f"Permission denied for {file_path}")
        return None
    except OSError as e:
        print(f"OS error occurred: {e}")
        return None


def read_bedbase(file_path: str) -> Optional[pd.DataFrame]:
    """Reads a bedbase file into a pandas DataFrame. Bedbase files are like
    bed files but each base is shown instead of genomic windows.

    Args:
        file_path (str): The path to the bedbase file.

    Returns:
        Optional[pandasDataFrame]: The DataFrame containing the bedbase data,
        or None if an error occurred.
    """
    try:
        bedbase = pd.read_table(file_path, sep="\t")
        if bedbase.shape[1] != 3:
            raise ValueError(f"Bedbase file at {file_path}"
                             "does not have exactly 3 columns.")
        bedbase.columns = ["CHR", "BASE", "SCORE"]
        return bedbase
    except (FileNotFoundError, IOError):
        print(f"{file_path} does not exist or could not be read.")
        return None
    except IsADirectoryError:
        print(f"{file_path} is a directory.")
        return None
    except PermissionError:
        print(f"Permission denied for {file_path}")
        return None
    except OSError as e:
        print(f"OS error occurred: {e}")
        return None


def read_narrow_peak_file(file_path: str) -> Optional[pd.DataFrame]:
    """Reads a BED3+7 file from MACS into a pandas DataFrame.

    Args:
        file_path (str): The path to the bedbase file.

    Returns:
        Optional[pd.DataFrame]: The DataFrame containing the bed data,
        or None if an error occurred.
    """
    try:
        bed = pd.read_table(file_path, sep="\t", skiprows=1)
        if bed.shape[1] != 10:
            raise ValueError(f"Bed file at {file_path}"
                             "does not have exactly 10 columns.")
        bed.columns = ["CHR", "START", "END", "NAME", "SCORE", "STRAND",
                       "SIGNAL_VALUE", "PVALUE", "QVALUE", "SUMMIT"]
        return bed
    except (FileNotFoundError, IOError):
        print(f"{file_path} does not exist or could not be read.")
        return None
    except IsADirectoryError:
        print(f"{file_path} is a directory.")
        return None
    except PermissionError:
        print(f"Permission denied for {file_path}")
        return None
    except OSError as e:
        print(f"OS error occurred: {e}")
        return None


def write_file(data: pd.DataFrame, file_path: str) -> None:
    """Writes a pandas DataFrame to a file in tab-separated format.

    Args:
        data (pandasDataFrame): The DataFrame to write to the file.
        file_path (str): The path to the output file.
    """
    try:
        with open(file_path, 'w') as file:
            data.to_csv(file, sep="\t", header=False, index=False)
    except (FileNotFoundError, IOError):
        print(f"{data} could not be written to {file_path}.")
    except IsADirectoryError:
        print(f"{file_path} is a directory.")
    except PermissionError:
        print(f"Permission denied for {file_path}")
    except OSError as e:
        print(f"OS error occurred: {e}")