from pathlib import Path
import urllib.request 


def download(url: str, filename: str) -> Path:
    """
    Download a file to the exact path given by `filename`.
    Ensures that the parent directory exists.
    """
    dest_path = Path(filename)
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    if not dest_path.exists():
        print(f"Downloading: {dest_path.name}")
        urllib.request.urlretrieve(url, dest_path)
    else:
        print(f"Found existing file: {dest_path}")

    return dest_path