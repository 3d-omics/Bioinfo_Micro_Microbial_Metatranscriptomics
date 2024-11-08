RESULTS = Path("results")


PRE = RESULTS / "preprocess"
READS = PRE / "reads"
HOSTS = PRE / "hosts"
FASTP = PRE / "fastp"
KRAKEN2 = PRE / "kraken2"
RIBODETECTOR = PRE / "ribodetector/"
STAR_INDEX = PRE / "index"
STAR = PRE / "star"
CLEAN = PRE / "clean"

QUANT = Path("results/quantify")
MAGS = QUANT / "mags"
BOWTIE2_INDEX = QUANT / "index"
BOWTIE2 = QUANT / "bowtie2"
COVERM = QUANT / "coverm"
BEDTOOLS = QUANT / "bedtools"
HTSEQ = QUANT / "htseq"
SUBREAD = QUANT / "subread"

REPORT_STEP = Path("reports/by_step/")
REPORT_LIBRARY = Path("reports/by_library/")
