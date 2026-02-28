# Test Fixtures

## sacCer3.fa.gz

S. cerevisiae reference genome (SacCer3), 12,157,105 bp across 17 chromosomes.
Used by `bench_guide_design.py`, `recall_test.py`, and `test/saccer3.test.ts`.

To decompress:
```bash
gunzip -k test/fixtures/sacCer3.fa.gz
```

If missing, download from UCSC:
```bash
curl -L "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz" \
  -o test/fixtures/sacCer3.fa.gz
```

## test_genome.fa

Small synthetic 6-line FASTA with many ATGATG patterns. Used by unit tests.

## phiX174.fa

PhiX174 bacteriophage genome (~5.4 kb). Used by basic integration tests.
