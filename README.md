# fqplot: Fastq Quality plot

Quality, gc, error rate stats and plot from fastq files. Both gzipped and plain text files are supported for read1 and read2 individually.

## Getting started

Help message

```shell
fqplot 0.1.0
slyo <sean.lyo@outlook.com>
Fast paired fastq quality, base percent, error rate plot.

USAGE:
    fqplot --prefix <FILE> --read1 <FILE> --read2 <FILE>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -p, --prefix <FILE>    Output prefix.
    -1, --read1 <FILE>     Fastq of read1.
    -2, --read2 <FILE>     Fastq of read2.
```

To make it work, font `wqy-zenhei` need installing first. See directory `data`.

Test examples.

```shell
fqplot -p sample -1 tests/sample_R1.fq.gz -2 tests/sample_R2.fq.gz
```

## Todo

- [ ] Make `--read2` optional, when read2 is absent, only plot read1.
- [ ] Less time.

## Benchmark

For gziped fastq,

```
~ 3.5min/Gb
```
