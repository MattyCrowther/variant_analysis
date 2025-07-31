import subprocess
from pathlib import Path

def decompress_gzip(filepath: Path) -> Path:
    """Decompress a .gz file in-place and return the decompressed file path."""
    filepath = Path(filepath)

    if filepath.suffix != ".gz":
        raise ValueError(f"Expected .gz file, got: {filepath}")

    decompressed_path = filepath.with_suffix("")

    if decompressed_path.exists():
        print(f"Already decompressed: {decompressed_path}")
        return decompressed_path

    print(f"Decompressing: {filepath}")
    subprocess.run(["gunzip", str(filepath)], check=True)

    return decompressed_path
