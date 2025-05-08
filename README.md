# csci4444_final

Note: the BAM file was too large to include in the repository.

## Usage
```
usage: microdna.py [-h] <-f> FILE [-t] THRESHOLD [-c] CUTOFF

required arguments:
  -f, --file            file to read
optional arguments:
  -h, --help            show this help message and exit
  -t, --threshold       minimum observations of a read to continue
  -c, --cutoff          score reporting cutoff
```

## Examples
```
python3 src/microdna.py -f "data/SRR413984.sorted.NC_000001.10.bam"
```

```
1. 121484503 -- 121484749: Score 3.8178

2. 121484614 -- 121484955: Score 3.7541

3. 121484200 -- 121484513: Score 3.7040

4. 121484599 -- 121484955: Score 3.7033

5. 121484602 -- 121484955: Score 3.6371

6. 121485043 -- 121485412: Score 3.5008

7. 121484060 -- 121484419: Score 3.4162

8. 121484781 -- 121485152: Score 3.4114

9. 36780729 -- 36780995: Score 3.3953

10. 121484601 -- 121484955: Score 3.3346
```
