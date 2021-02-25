# fqplot: Fastq Quality plot

Quality, gc, error rate stats and plot form fastq file.

## Getting started

Help message

```shell
fqplot
slyo <sean.lyo@outlook.com>
Fastq(s) quality, base percent, error rate plotting.

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

## Todo

- [ ] Make `--read2` optional, when read2 is absent, only plot read1.
- [ ] Less time.

## Benchmark

For gziped fastq,

```
~ 3.5min/Gb
```
