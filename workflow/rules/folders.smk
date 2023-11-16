READS = Path("results/reads/")
REFERENCE = Path("results/reference/")

PRE = Path("results/preprocessing")
FASTP = PRE / "fastp"
KRAKEN2 = PRE / "kraken2"
RIBODETECTOR = PRE / "ribodetector/"
INDEX = PRE / "index"
STAR = PRE / "star"

QUANT = Path("results/quantification")
BOWTIE2 = QUANT / "bowtie2"
COVERM = QUANT / "coverm"

REPORT_STEP = Path("reports/by_step/")
REPORT_LIBRARY = Path("reports/by_library/")
