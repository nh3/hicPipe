# hicPipe

**Map and process pair-end C-seq reads**
* _Split read mapping for longer reads_
* _Iterative trimming and mapping for shorter reads_
<p></p>

    Usage:
        hicPipe.py (split|iter|short) [options] <genome> <read1> <read2> [<output>]

    Options:
        <output>             output prefix, append "." if not empty [default: ]
        --dryRun             print commands without execution
        --resume             start from where it failed last time
        --start <INT>        index of the first step to run
        --end <INT>          index of the last step to run
        --nThread <INT>      number of threads [default: 1]
        --exChrom <INT>      0-based index of excluded chromosomes in <genome>, index of MT by default [default: auto]
        --minD <INT>         minimum mapping distance (bp) for informative read pairs [default: 600]
        --minMQ <INT>        minimum MAPQ [default: 30]
        --maxNM <INT>        maximum allowed number of mismatches per read [default: 2]

    split options:
        --maxGap <INT>       maximum allowed number of bases per read not covered by any alignment [default: 19]
        --maxOverlap <INT>   maximum allowed number of bases per read covered by multiple alignments [default: 5]
        --multiSplit         allow more than two split parts per read [default: False]

    iter options:
        --readLen <INT>      typical read length [default: 100]
        --trimStep <INT>     trim step [default: 4]
        --minReadLen <INT>   minimum read length after trimming [default: 20]
