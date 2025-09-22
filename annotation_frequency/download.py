
import subprocess
import os

def download(url,output):
    # Download and index the reference VCF
    if os.path.isfile(output):
        return
    os.makedirs(output.parent, exist_ok=True)
    subprocess.run(["wget", url, "-O", str(output)], check=True)

    if output.suffix.lower() == "ycf":
        subprocess.run(["tabix", "-p", "vcf", str(output)], check=True)
    