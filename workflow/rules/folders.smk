RESULTS = Path("results")

PRE = RESULTS / "preprocess"
PRE_READS = PRE / "reads"
PRE_HOSTS = PRE / "hosts"
PRE_FASTP = PRE / "fastp"
PRE_KRAKEN2 = PRE / "kraken2"
PRE_RIBODETECTOR = PRE / "ribodetector/"
PRE_INDEX = PRE / "index"
PRE_STAR = PRE / "star"
PRE_CLEAN = PRE / "clean"

QUANT = RESULTS / "quantify"
QUANT_MAGS = QUANT / "mags"
QUANT_INDEX = QUANT / "index"
QUANT_BOWTIE2 = QUANT / "bowtie2"
QUANT_SUBREAD = QUANT / "subread"
