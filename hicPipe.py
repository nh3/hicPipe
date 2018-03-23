#!/usr/bin/env python
'''
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
    --exReg <FILE>       excluded region, in BED format (3 columns minimum)
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
'''

from __future__ import print_function
import sys
import os.path
import os
import subprocess as sbp
import re
from os import getcwd
from datetime import datetime
import signal
import logging
from subprocess import CalledProcessError
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
logging.basicConfig(
        level=logging.WARN,
        format='%(asctime)s; %(levelname)s; %(funcName)s(); %(message)s',
        datefmt='%y-%m-%d %H:%M:%S')

def currTime():
    return datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')


class Step(object):
    def __init__(self, stepName, cmd, defaultArgs=None):
        self.name = stepName
        self.cmdTemplate = cmd
        self.defaultArgs = defaultArgs
        self.cmd = None
        self.args = None
        self.execute = False
        self.retcode = None

    def assembleCmd(self, **args):
        if self.cmd is None or len(args) > 0:
            for k in self.defaultArgs:
                if k not in args:
                    args[k] = self.defaultArgs[k]
            assert 'inputs' in args, '{} lack <inputs>'.format(self.name)
            assert 'output' in args, '{} lack <output>'.format(self.name)
            self.cmd = self.cmdTemplate.format(**args)
            self.args = args

    def run(self, force=False, dryRun=True, stdout=None, stderr=None, **args):
        self.assembleCmd(**args)
        if not dryRun:
            self.execute = True
            self.exeTime = currTime()
            inputs = self.args['inputs']
            inputModifyTime = self._checkInput(inputs)
            if inputModifyTime is None:
                self.execute = False
                logging.error('{} (inputs) does not exist or is empty'.format(inputs))
                self.retcode = 1
            output = self.args['output']
            outputModifyTime = self._checkOutput(output)
            if outputModifyTime is not None and outputModifyTime > inputModifyTime:
                if not force:
                    self.execute = False
                    logging.warn('{} (output) exists and is newer than {} (inputs), skip {} without forcing'.format(output, inputs, self.name))
                    self.retcode = 1
                else:
                    logging.warn('{} (output) exists and is newer than {} (inputs), run {} by force'.format(output, inputs, self.name))
            if self.execute:
                self.execute = False
                try:
                    self.retcode = sbp.call(self.cmd, stdout=stdout, stderr=stderr, shell=True, executable='/bin/bash')
                except KeyboardInterrupt:
                    self.retcode = 1
                    self.retTime = currTime()
                    raise KeyboardInterrupt
            self.retTime = currTime()
        else:
            return self.cmd

    def _checkInput(self, inputs):
        if type(inputs) is str:
            inputs = [inputs]
        elif type(inputs) is dict:
            inputs = inputs.values()
        mtime = None
        if all([os.path.exists(f) and os.path.getsize(f)>0 for f in inputs]):
            mtime = max([os.path.getmtime(f) for f in inputs])
        return mtime

    def _checkOutput(self, output):
        if type(output) is str:
            output = [output]
        elif type(output) is dict:
            output = output.values()
        mtime = None
        if any([os.path.exists(f) for f in output]):
            mtime = min([os.path.getmtime(f) for f in output if os.path.exists(f)])
        return mtime


class HiCProject(object):
    def __init__(self, genome, read1, read2, output,
                 split=False, iter=False, short=False,
                 nThread=1, exChrom=None, exReg=None, minMQ=30, maxNM=2, minD=600,
                 maxGap=19, maxOverlap=5, multiSplit=False, readLen=None, trimStep=4, minReadLen=20):
        self.genome = genome
        self.read1 = read1
        self.read2 = read2
        self.prefix = output
        self.outdir = os.path.dirname(output)
        self.nThread = nThread
        if exChrom is None or exChrom == 'auto':
            self.exChrom = self.findMtInGenome()
        else:
            self.exChrom = int(exChrom)
        self.exReg = exReg
        self.minMQ = minMQ
        self.maxNM = maxNM
        self.minD = minD
        if split:
            self.mode = 'split'
            self.maxGap = int(maxGap)
            self.maxOverlap = int(maxOverlap)
            self.multiSplit = multiSplit
        elif iter:
            self.mode = 'iter'
            self.readLen = int(readLen)
            self.trimStep = int(trimStep)
            self.minReadLen = int(minReadLen)
        else:
            self.mode = 'short'
        self.preparePrefix()
        self.prepareFileNames()
        self.pipeline = None

    def preparePrefix(self):
        path = self.prefix
        if '/' in path:
            if path.endswith('/'):
                self.prefix = path
            else:
                self.prefix = path + '.'
        elif len(path) > 0:
            self.prefix = '{}/{}.'.format(path,path) 
        else:
            self.prefix = ''

    def prepareFileNames(self):
        self.chromSize = os.path.join(self.outdir, os.path.basename(self.genome).replace('.fa','.chromSize.txt'))
        self.fastq1    = '{}read1.fastq.gz'.format(self.prefix)
        self.fastq2    = '{}read2.fastq.gz'.format(self.prefix)
        self.aligned1  = '{}read1.bam'.format(self.prefix)
        self.aligned2  = '{}read2.bam'.format(self.prefix)
        self.merged    = '{}merged.sortn.bam'.format(self.prefix)
        self.paired    = '{}exMT.mapq.fixmate.paired.sortp.bam'.format(self.prefix)
        self.markdup   = '{}exMT.mapq.fixmate.paired.sortp.markdup.bam'.format(self.prefix)
        self.valid     = '{}valid.bam'.format(self.prefix)
        self.validBW   = '{}valid.bw'.format(self.prefix)
        self.transInfo = '{}info.trans.bam'.format(self.prefix)
        self.cisInfo   = '{}info.cis.bam'.format(self.prefix)
        self.cisInfoBW = '{}info.cis.bw'.format(self.prefix)
        self.seqStats  = '{}stats.txt'.format(self.prefix)
        self.log       = '{}log'.format(self.prefix)
        self.statusLog = '{}status'.format(self.prefix)

    def findMtInGenome(self):
        fai = self.genome + '.fai'
        assert os.path.exists(fai), 'genome index ({}) not found'.format(fai)
        mtPattern = re.compile(r'(?:chr|CHROMOSOME_)M(?:[Tt]|[Tt]_DNA)?')
        with open(fai) as f:
            chroms = [(line.rstrip().split('\t',1))[0]for line in f]
            idxMT = [i for i,c in enumerate(chroms) if mtPattern.match(c)]
        if len(idxMT) == 1:
            idxMT = idxMT[0]
        elif len(idxMT) > 1:
            idxMT = idxMT[-1]
            logging.warn('Multiple MT detected in {}, the last occurrence ({}) is taken'.format(self.genome, chroms[idxMT]))
        else:
            idxMT = -1
            logging.warn('No MT detected in {}'.format(self.genome))
        return idxMT

    def preparePipeline(self):
        READ1 = '0x40'
        READ2 = '0x80'
        steps = []
        steps.append(self.warmUp('warmUp', [self.read1, self.read2], [self.fastq1, self.fastq2]))
        if self.mode == 'split':
            steps.append(self.bwaMem('bwaMem1', inputs=self.fastq1, output=self.aligned1, flag=READ1))
            steps.append(self.bwaMem('bwaMem2', inputs=self.fastq2, output=self.aligned2, flag=READ2))
            steps.append(self.mergeReads('mergeReads', inputs=[self.aligned1,self.aligned2], output=self.merged))
            steps.append(self.pairSplitReads('pairSplit', inputs=self.merged, output=self.paired, exReg=self.exReg))
        elif self.mode == 'iter':
            bams = []
            i = 0
            for rlen in range(self.readLen, self.minReadLen-1, -self.trimStep):
                if rlen < self.readLen:
                    if i == 1:
                        prevBam1 = '{}read1.bam'.format(self.prefix)
                        prevBam2 = '{}read2.bam'.format(self.prefix)
                    else:
                        prevBam1 = '{}read1.{}bp.bam'.format(self.prefix, rlen+self.trimStep)
                        prevBam2 = '{}read2.{}bp.bam'.format(self.prefix, rlen+self.trimStep)
                    fq1 = '{}read1.{}bp.fastq.gz'.format(self.prefix, rlen)
                    fq2 = '{}read2.{}bp.fastq.gz'.format(self.prefix, rlen)
                    bam1 = '{}read1.{}bp.bam'.format(self.prefix, rlen)
                    bam2 = '{}read2.{}bp.bam'.format(self.prefix, rlen)
                    steps.append(self.trimUnmapped('trim1-{}'.format(rlen), inputs=prevBam1, output=fq1, length=rlen))
                    steps.append(self.trimUnmapped('trim2-{}'.format(rlen), inputs=prevBam2, output=fq2, length=rlen))
                else:
                    fq1 = self.fastq1
                    fq2 = self.fastq2
                    bam1 = self.aligned1
                    bam2 = self.aligned2
                steps.append(self.bwaAln('bwaAln1-{}'.format(rlen), inputs=fq1, output=bam1, flag=READ1))
                steps.append(self.bwaAln('bwaAln2-{}'.format(rlen), inputs=fq2, output=bam2, flag=READ2))
                bams.append(bam1)
                bams.append(bam2)
                i += 1
            steps.append(self.mergeReads('mergeReads', inputs=bams, output=self.merged))
            steps.append(self.pairReads('pairReads', inputs=self.merged, output=self.paired))
        else:
            steps.append(self.bwaAln('bwaAln1', inputs=self.fastq1, output=self.aligned1, flag=READ1))
            steps.append(self.bwaAln('bwaAln2', inputs=self.fastq2, output=self.aligned2, flag=READ2))
            steps.append(self.mergeReads('mergeReads', inputs=[self.aligned1,self.aligned2], output=self.merged))
            steps.append(self.pairReads('pairReads', inputs=self.merged, output=self.paired))
        steps.append(self.markDup('markDup', inputs=self.paired, output=self.markdup))
        steps.append(self.removeDup('removeDup', inputs=self.markdup, output=self.valid))
        steps.append(self.extractTransInfo('getTrans', inputs=self.valid, output=self.transInfo))
        steps.append(self.extractCisInfo('getCis', inputs=self.valid, output=self.cisInfo))
        steps.append(self.gatherStats('gatherStats', inputs=[self.fastq1, self.fastq2,
                                                             self.merged,self.paired,self.valid,
                                                             self.cisInfo,self.transInfo], output=self.seqStats))
        steps.append(self.makeBigwig('makeValidBW', inputs=self.valid, output=self.validBW))
        steps.append(self.makeBigwig('makeCisInfoBW', inputs=self.cisInfo, output=self.cisInfoBW))
        self.pipeline = steps

    def checkStatus(self):
        assert self.pipeline is not None and len(self.pipeline) > 0, 'pipeline is not defined, run preparePipeline() first'
        if os.path.exists(self.statusLog):
            self.status = {}
            with open(self.statusLog) as f:
                for line in f:
                    time,i,st = line.rstrip().split('\t')
                    self.status[int(i)] = int(st)
        else:
            self.status = {i:None for i in xrange(len(self.pipeline))}

    def printSteps(self):
        assert self.pipeline is not None and len(self.pipeline) > 0, 'pipeline is not defined, run preparePipeline() first'
        print('# {}'.format(currTime()))
        print('# CWD: {}\n'.format(getcwd()))
        for i,step in enumerate(self.pipeline):
            cmd = step.run(dryRun=True)
            print('# {}: {}\n{}\n'.format(i, step.name, cmd))

    def run(self, resume=False, start=None, end=None):
        assert self.pipeline is not None and len(self.pipeline) > 0, 'pipeline is not defined, run preparePipeline() first'
        with open(self.statusLog, 'a') as fstatus, open(self.log, 'w') as flog:
            if start is not None:
                for i,step in enumerate(self.pipeline):
                    if i >= start and (end is None or i<= end):
                        try:
                            step.run(dryRun=False, force=True, stdout=flog, stderr=flog)
                        except KeyboardInterrupt:
                            raise KeyboardInterrupt
                        finally:
                            print('{}\t{}\t{}'.format(step.retTime, i, step.retcode), file=fstatus)
                        assert step.retcode == 0
            elif resume:
                self.checkStatus()
                x = [i for i,st in self.status.iteritems() if st != 0]
                if len(x) > 0:
                    start = min(x)
                    for i,step in enumerate(self.pipeline):
                        if i >= start:
                            try:
                                step.run(dryRun=False, stdout=flog, stderr=flog)
                            except KeyboardInterrupt:
                                raise KeyboardInterrupt
                            finally:
                                print('{}\t{}\t{}'.format(step.retTime, i, step.retcode), file=fstatus)
                            assert step.retcode == 0
            else:
                for i,step in enumerate(self.pipeline):
                    try:
                        step.run(dryRun=False, stdout=flog, stderr=flog)
                    except KeyboardInterrupt:
                        raise KeyboardInterrupt
                    finally:
                        print('{}\t{}\t{}'.format(step.retTime, i, step.retcode), file=fstatus)
                    assert step.retcode == 0

    def writeConfig(self, filename=None):
        entries = []
        for k,v in vars(self).iteritems():
            if type(v) is not dict:
                entries.append('{} = {}'.format(k,v))
            else:
                for dk,dv in v.iteritems():
                    entries.append('{} = {}'.format(dk,dv))
        if filename is None:
            filename = '{}config'.format(self.prefix)
        entries = sorted(entries)
        with open(filename, 'w') as f:
            print('# {}'.format(currTime()), file=f)
            for e in entries:
                print(e, file=f)

    def readConfig(self, filename):
        pass

    def assignDefaultArgs(self, args):
        if 'self' in args:
            del args['self']
        if 'stepName' in args:
            del args['stepName']
        for k in args:
            if hasattr(self, k):
                args[k] = getattr(self, k)
        return args

    def warmUp(self, stepName, inputs, output, executables=[]):
        args = self.assignDefaultArgs(locals())
        dir = os.path.dirname(self.prefix)
        cmds = []
        for prog in executables:
            cmds.append('which {}'.format(prog))
        if not os.path.exists(dir):
            cmds.append('mkdir {}'.format(dir))
        for read,fq in zip(inputs, output):
            if not os.path.lexists(fq) and fq != read:
                cmds.append('ln -s {} {}'.format(os.path.abspath(read), fq))
        cmds.append('cut -f1,2 {} > {}'.format(self.genome+'.fai', self.chromSize))
        cmd = ' && '.join(cmds)
        return Step(stepName, cmd, defaultArgs=args)

    def bwaMem(self, stepName, inputs, output, flag, nThread=1, genome=None, otherBwaOpt=''):
        args     = self.assignDefaultArgs(locals())
        align    = "bwa mem -t {nThread} -Y {otherBwaOpt} {genome} {inputs}"
        setFlag  = "awk -v OFS='\\t' '{{if ($0!~/^@/) {{$2=or($2, {flag})}}; print}}'"
        sam2Bam  = "samtools view -@ {nThread} -b -"
        sortName = "sambamba sort -n -o {output} /dev/stdin"
        cmd      = ' | '.join([align, setFlag, sam2Bam, sortName])
        return Step(stepName, cmd, defaultArgs=args)

    def bwaAln(self, stepName, inputs, output, flag, nThread=1, genome=None, otherBwaOpt=''):
        args     = self.assignDefaultArgs(locals())
        align    = "bwa aln -t {nThread} {otherBwaOpt} {genome} {inputs} | bwa samse {genome} - {inputs}"
        setFlag  = "awk -v OFS='\\t' '{{if ($0!~/^@/) {{$2=or($2, {flag})}}; print}}'"
        sam2Bam  = "samtools view -@ {nThread} -b -"
        sortName = "sambamba sort -n -o {output} /dev/stdin"
        cmd      = ' | '.join([align, setFlag, sam2Bam, sortName])
        return Step(stepName, cmd, defaultArgs=args)

    def trimUnmapped(self, stepName, inputs, output, length, nThread=1):
        args     = self.assignDefaultArgs(locals())
        filters  = "samtools view -u -f 0x4 {inputs}"
        bam2Fq   = "samtools bam2fq -n -s /dev/stdout -"
        trim     = "trimFq -5 {length}"
        compress = "pigz -p {nThread} -c > {output}"
        cmd      = ' | '.join([filters, bam2Fq, trim, compress])
        return Step(stepName, cmd, defaultArgs=args)

    def mergeReads(self, stepName, inputs, output, nThread=1):
        args     = self.assignDefaultArgs(locals())
        inputs   = ' '.join(['{{inputs[{}]}}'.format(i) for i in xrange(len(inputs))])
        merge    = "samtools merge -@ {nThread} - " + inputs
        sort     = "sambamba sort -t {nThread} -F 'not unmapped' -n -o {output} /dev/stdin"
        cmd      = ' | '.join([merge, sort])
        return Step(stepName, cmd, defaultArgs=args)

    def pairSplitReads(self, stepName, inputs, output, nThread=1, exChrom=-1, exReg=None, minMQ=30, maxNM=2, maxGap=19, maxOverlap=5, minD=600):
        args     = self.assignDefaultArgs(locals())
        if exReg is None:
            filters = "sambamba view -t {nThread} -F 'ref_id!={exChrom}' -h {inputs}"
        else:
            filters = "bedtools intersect -v -a {inputs} -b {exReg} | sambamba view -F 'ref_id!={exChrom}' -h /dev/stdin"
        pair     = "pairSplitRead -m {minMQ} -n {maxNM} -g {maxGap} -v {maxOverlap} -d {minD} -i <(%s)" % filters
        fixmate  = "samtools fixmate -r -p -O bam - -"
        sort     = "sambamba sort -t {nThread} -o {output} /dev/stdin"
        cmd      = ' | '.join([pair, fixmate, sort])
        return Step(stepName, cmd, defaultArgs=args)

    def pairReads(self, stepName, inputs, output, nThread=1, exChrom=-1, minMQ=30, maxNM=2):
        args     = self.assignDefaultArgs(locals())
        filters  = "sambamba view -t {nThread} -F 'ref_id!={exChrom} and mapping_quality>={minMQ} and [NM]<={maxNM}' -f bam {inputs}"
        fixmate  = "samtools fixmate -r -p -O bam - -"
        pair     = "sambamba sort -t {nThread} -F 'paired' -o {output} /dev/stdin"
        cmd      = ' | '.join([filters, fixmate, pair])
        return Step(stepName, cmd, defaultArgs=args)

    def markDup(self, stepName, inputs, output, nThread=1):
        args     = self.assignDefaultArgs(locals())
        cmd      = "sambamba markdup -t {nThread} --hash-table-size=1000000 --overflow-list-size=1000000 {inputs} {output}"
        return Step(stepName, cmd, defaultArgs=args)

    def removeDup(self, stepName, inputs, output, nThread=1):
        args     = self.assignDefaultArgs(locals())
        cmd      = "samtools view -@ {nThread} -F 0x400 -b -o {output} {inputs}"
        return Step(stepName, cmd, defaultArgs=args)

    def extractTransInfo(self, stepName, inputs, output, nThread=1):
        args     = self.assignDefaultArgs(locals())
        cmd      = "sambamba view -t {nThread} -F 'chimeric' -f bam -o {output} {inputs}"
        return Step(stepName, cmd, defaultArgs=args)

    def extractCisInfo(self, stepName, inputs, output, nThread=1, minD=600):
        args     = self.assignDefaultArgs(locals())
        cmd      = "sambamba view -t {nThread} -F 'template_length>{minD} or template_length<-{minD}' -f bam -o {output} {inputs}"
        return Step(stepName, cmd, defaultArgs=args)

    def makeBigwig(self, stepName, inputs, output, nThread=1, chromSize=None):
        args     = self.assignDefaultArgs(locals())
        cmd      = "bam2bw -c {chromSize} -o {output} {inputs}"
        return Step(stepName, cmd, defaultArgs=args)

    def gatherStats(self, stepName, inputs, output, nThread=1, exChrom=-1):
        args     = self.assignDefaultArgs(locals())
        inputs   = ' '.join(['{{inputs[{}]}}'.format(i) for i in xrange(len(inputs))])
        cmd      = "gatherCSeqStats " + inputs + " {exChrom} {output} {nThread}"
        return Step(stepName, cmd, defaultArgs=args)

    def makeContactMatrix(self, stepName, inputs, output, nThread=1, chromSize=None):
        args     = self.assignDefaultArgs(locals())
        sam2Mat  = 'samtools view {inputs} {{}} | sam2contact -g {chromSize} -w {resolution} -s {step} -r {{}} -o {prefix}'
        cmd      = "cut -f1 {chromSize} | parallel -j {nThread} '%s'" % sam2Mat
        return Step(stepName, cmd, defaultArgs=args)

    # TODO
    def plotContactMap(self, stepName, inputs, output):
        args     = self.assignDefaultArgs(locals())
        plotCMap = 'plotCMap -m {inputs[0]} -b {inputs[1]} -o {output}'
        return Step(stepName, cmd, defaultArgs=args)

    # TODO
    def callSignificantInteractions(self):
        pass

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def check_executables():
    required_exec = ['bwa', 'samtools', 'sambamba', 'pigz', 'trimFq', 'pairSplitRead', 'gatherCSeqStats', 'bam2bw']
    missing_exec = False
    for exe in required_exec:
        exe_path = which(exe)
        if exe_path is None:
            missing_exec = True
            logging.error('*{} not found*'.format(exe))
        else:
            logging.info('{} found at {}'.format(exe, exe_path))
    if missing_exec:
        return 1
    else:
        return 0

def main(args):
    logging.debug(args)
    dryRun = args['dryRun']
    resume = args['resume']
    start  = args['start']
    end    = args['end']
    del args['dryRun']
    del args['resume']
    del args['start']
    del args['end']

    if check_executables() != 0:
        return 1
    hic = HiCProject(**args)
    hic.writeConfig()
    hic.preparePipeline()
    if dryRun:
        hic.printSteps()
    else:
        if start is not None:
            start = int(start)
        if end is not None:
            end = int(end)
        hic.run(resume=resume, start=start, end=end)
    logging.info('all done')
    return 0


if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)
    args = {k.lstrip('-<').rstrip('>'):args[k] for k in args}
    try:
        main(args)
    except KeyboardInterrupt:
        logging.warning('Interrupted')
        sys.exit(1)
