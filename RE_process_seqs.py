import os
import sys
import time
import argparse
import re
import logging
import subprocess
import fnmatch

sys.path.insert(0, "/data/boehm/group/pipelines/ruffus")
from ruffus import *
from ruffus.proxy_logger import *
from ruffus.combinatorics import *
from ruffus.drmaa_wrapper import run_job, run_job_using_drmaa, error_drmaa_job

from PIL import Image
import string

parser = argparse.ArgumentParser(prog="REpipe", version="0.1.0", description="processes RE amplicon sequences", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-ri", "--readIn", dest="readdir", action="store", default=False, help="input read folder")
parser.add_argument("-as", "--aSet", dest="aSet", action="store", default=False, help="amplicon set information")
parser.add_argument("-wd", "--wdir", dest="wdir", action="store", default=False, help="output folder")
parser.add_argument("-bs", "--batchSize", dest="bsize", action="store",metavar="INT",type=int,default=10, help="number of samples to process in parallel")
parser.add_argument("-nt", "--numThr", dest="nthreads", action="store",metavar="INT",type=int,default=8, help="number of threads to use per sample")
parser.add_argument("-ec", "--errorCorrect", dest="errorCorrect", action="store_true", help="run error correction of sequences (experimental)")
parser.add_argument("--touchOnly", dest="touchOnly", action="store_true", help="only touch files")
parser.add_argument("--target_tasks", dest="target_tasks", action="store",default=[], help="target tasks")
parser.add_argument("--forcedtorun_tasks", dest="forcedtorun_tasks", action="store",default=[], help="forced to run tasks")

args = parser.parse_args()

#setup central working directory
wdir=args.wdir
if not os.path.exists(wdir):
    os.makedirs(wdir)
os.chdir(wdir)

#setup logging
logger = logging.getLogger(__name__)
fhandler = logging.FileHandler(filename=os.path.join(wdir,'pipeline.log'), mode='a')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
logger.setLevel(logging.DEBUG)

logger.debug(subprocess.check_output('echo $DRMAA_LIBRARY_PATH',shell=True))

import drmaa

#initiate 1 drmaa session for the whole pipeline
mySession=drmaa.Session()
mySession.initialize()

#identify pipeline input files
readdir=args.readdir

libset2 = []
# Walk through directory
for dName, sdName, fList in os.walk(readdir):
    for fileName in fList:
        if fnmatch.fnmatch(fileName, "*fastq.gz"): # Match search string
            libset2.append(os.path.join(dName, fileName))

libset2_R1=filter(lambda x:'_R1.fastq.gz' in x, libset2)
libset2_R1.sort()
libset2_R2=filter(lambda x:'_R2.fastq.gz' in x, libset2)
libset2_R2.sort()
read_root=[ re.sub('_R1.fastq.gz','',R1f) for R1f in libset2_R1 ]
INfiles=list(zip(libset2_R1,libset2_R2))
    
logger.debug(INfiles)  

##################PATHS TO EXECUTABLES###############################################################
FQCpath='/package/FastQC-0.11.3'
cutpath='/package/cutadapt-1.9.1/bin'
prinpath='/data/boehm/sikora/tools/prinseq-lite-0.20.4' #/data/boehm/group/pipelines/tools/COMPLETE_PATH
flashpath='/data/boehm/sikora/tools/FLASH-1.2.11'
blastpath='/data/boehm/sikora/tools/ncbi-blast-2.2.30+/bin'
FXpath='/package/fastx_toolkit_0.0.13'
Rpath='/package/R-3.3.1/bin'

###############configure amplicon and primer sequences##############################################
if args.aSet=='171214RE':
    comb_db='/home/sikora/works/boehm/ines/RE_sequences/171214RE/RE8.RE13.comb.db'
    RE8e6g='/home/sikora/works/boehm/ines/RE_sequences/171214RE/RE8.e6.cutadapt.g.fasta'
    RE8e6a='/home/sikora/works/boehm/ines/RE_sequences/171214RE/RE8.e6.cutadapt.a.fasta'
    RE8e7g='/home/sikora/works/boehm/ines/RE_sequences/171214RE/RE8.e7.cutadapt.g.fasta'
    RE8e7a='/home/sikora/works/boehm/ines/RE_sequences/171214RE/RE8.e7.cutadapt.a.fasta'
    RE13e2g='/home/sikora/works/boehm/ines/RE_sequences/171214RE/RE13.e2.cutadapt.g.fasta'
    RE13e2a='/home/sikora/works/boehm/ines/RE_sequences/171214RE/RE13.e2.cutadapt.a.fasta'
    RE13e3g='/home/sikora/works/boehm/ines/RE_sequences/171214RE/RE13.e3.cutadapt.g.fasta'
    RE13e3a='/home/sikora/works/boehm/ines/RE_sequences/171214RE/RE13.e3.cutadapt.a.fasta'
    gList=[RE8e6g,RE8e7g,RE13e2g,RE13e3g]
    aList=[RE8e6a,RE8e7a,RE13e2a,RE13e3a]
    RoutList=['.RE8.e6.hits.fasta','.RE8.e7.hits.fasta','.RE13.e2.hits.fasta','.RE13.e3.hits.fasta']

elif args.aSet=='170905RE':
    comb_db='/home/sikora/works/boehm/ines/RE_sequences/170905RE/RE8.RE13.comb.db'
    RE8e6g='/home/sikora/works/boehm/ines/RE_sequences/170905RE/RE8.e6.cutadapt.g.fasta'
    RE8e6a='/home/sikora/works/boehm/ines/RE_sequences/170905RE/RE8.e6.cutadapt.a.fasta'
    RE13e2g='/home/sikora/works/boehm/ines/RE_sequences/170905RE/RE13.e2.cutadapt.g.fasta'
    RE13e2a='/home/sikora/works/boehm/ines/RE_sequences/170905RE/RE13.e2.cutadapt.a.fasta'
    gList=[RE8e6g,RE13e2g]
    aList=[RE8e6a,RE13e2a]
    RoutList=['.RE8.e6.hits.fasta','.RE13.e2.hits.fasta']


##########define folder structure ##############################
cutout=os.path.join(wdir,'reads_cut')
fqcout=os.path.join(wdir,'fastqc_cut')
#blastout=os.path.join(wdir,"megablast")
Rout=os.path.join(wdir,"workspaceR")
Fout=os.path.join(wdir,"final_sequences")

########adapter-trim and BQ-trim sequencing reads########
@mkdir(cutout,os.path.join(cutout,'logs'))
@transform(INfiles,suffix('_R1.fastq.gz'),['_R1.fastq.gz','_R2.fastq.gz'],output_dir=cutout)
def cut_reads(infiles, outfiles):
    ii1=infiles[0]
    ii2=infiles[1]
    oo1=outfiles[0]
    oo2=outfiles[1]
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(cutpath,'cutadapt') + ' -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 30  -n 5  -o ' + oo1 + ' -p ' + oo2 + ' ' + ii1 + ' ' + ii2
    logger.info(bshcmd)       
    with open(os.path.join(cutout,"logs","%s.cut_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.cut_reads.err" % read_root),'w+') as stderrF:
        try:
        
            stdout_res, stderr_res  = run_job(cmd_str   = bshcmd,
                                      job_name          = 'cut_reads',
                                      logger            = logger,
                                      drmaa_session     = mySession,
                                      run_locally       = False,
                                      working_directory = os.getcwd(),
                                      job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logger.error("Cut_reads error: %s" % err)
            raise
        else:
           logger.info('Adapter trimming complete')

@transform(cut_reads,suffix('_R1.fastq.gz'),['_prin_1.fastq.gz','_prin_2.fastq.gz'],output_dir=cutout)
def trim_BQ(infiles,outfiles):
    ii1=infiles[0]
    ii2=infiles[1]
    oo1=outfiles[0]
    oo2=outfiles[1]
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(ii1))
    uzcmd1='zcat -v '+ ii1 + ' > ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1)))
    uzcmd2='zcat -v '+ ii2 + ' > ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2)))
    bshcmd='perl '+ os.path.join(prinpath,'prinseq-lite.pl') + ' -fastq ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1))) + ' -fastq2 ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2))) + ' -out_good ' + os.path.join(cutout,re.sub('_1.fastq.gz','',os.path.basename(oo1))) +' -trim_qual_right 24 -trim_qual_type mean -trim_qual_window 6 -trim_qual_step 3 -min_len 50 -ns_max_p 10 -out_bad null'
    zcmd1='gzip -c '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(oo1))) + ' > ' + oo1
    zcmd2='gzip -c '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(oo2))) + ' > ' + oo2
    clcmd='rm -v '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(oo1))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(oo2)))
    cmd_all=';'.join([uzcmd1,uzcmd2,bshcmd,zcmd1,zcmd2,clcmd])
    logger.info(cmd_all)           
    with open(os.path.join(cutout,"logs","%s.BQtrim_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.BQtrim_reads.err" % read_root),'w+') as stderrF:
        try:
        
            stdout_res, stderr_res  = run_job(cmd_str   = cmd_all,
                                      job_name          = 'BQtrim_reads',
                                      logger            = logger,
                                      drmaa_session     = mySession,
                                      run_locally       = False,
                                      working_directory = os.getcwd(),
                                      job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logger.error("BQtrim_reads error: %s" % err)
            raise
        else:
           logger.info('Base quality trimming complete')
           
@mkdir(fqcout,os.path.join(fqcout,'logs'))            
@transform(trim_BQ,suffix('_1.fastq.gz'),output=['_1.zip','_1.html','_2.zip','_2.html'],output_dir=fqcout)
def postTrim_fqc(input_files,output_files):
    ii1 = input_files[0]
    ii2 = input_files[1]
    read_root=re.sub('_prin_1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(FQCpath,'fastqc ')+' --outdir ' + fqcout + ' -t '+ str(args.nthreads) + ' ' + ii1 + ' ' + ii2
    logger.info(bshcmd)       
    with open(os.path.join(fqcout,"logs","%s.post_fqc.out" % read_root),'w+') as stdoutF, open(os.path.join(fqcout,"logs","%s.post_fqc.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'post_fqc',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(args.nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Post_trim_fastqc error: %s" % err)
            raise
        else:
            logger.info('Post trim fastqc complete')

#######merge forward and reverse mates########
@transform(trim_BQ,suffix('_1.fastq.gz'),output='_flash.extendedFrags.fastq.gz',output_dir=cutout)
def merge_mates(input_files,output_file):
    ii1 = input_files[0]
    ii2 = input_files[1]
    oo = re.sub('.extendedFrags.fastq.gz','',os.path.basename(output_file))
    read_root=re.sub('_prin_1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(flashpath,'flash')+ ' -z -M 300 -t ' + str(args.nthreads) + ' -o '+ oo + ' -d ' + cutout + ' ' + ii1 + ' ' + ii2 
    logger.info(bshcmd)
    with open(os.path.join(cutout,"logs","%s.flash.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.flash.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'flash',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(args.nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Flash error: %s" % err)
            raise
        else:
            logger.info('Merging mates complete')


####modify read names to replace spaces " " with underscores "_"
@transform(merge_mates,suffix('_prin_flash.extendedFrags.fastq.gz'),output='_prin_flash.extendedFrags.sed.fasta',output_dir=cutout)#
def mod_Rnames(input_file,output_file):
    ii = input_file
    oo = output_file
    read_root=re.sub('_prin_flash.extendedFrags.fastq.gz','',os.path.basename(ii))
    bshcmd='zcat ' + ii + ' | sed \'s/\ /_/g\' - | awk \'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}\' -  > ' + oo 
    logger.info(bshcmd)
    with open(os.path.join(cutout,"logs","%s.sed.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.sed.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'sed',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Sed error: %s" % err)
            raise
        else:
            logger.info('Renaming reads complete')

########megablast against the amplicon sequences###########
#@mkdir(blastout,os.path.join(blastout,'logs'))            
@transform(mod_Rnames,suffix('_prin_flash.extendedFrags.sed.fasta'),output='.RE8.RE13.bla',output_dir=cutout)#
def mega_blast(input_file,output_file):
    ii = input_file
    oo = output_file
    read_root=re.sub('_prin_flash.extendedFrags.fasta','',os.path.basename(ii))
    bshcmd=os.path.join(blastpath,'blastn ')+ '-task megablast -num_threads ' + str(args.nthreads) + ' -db ' + comb_db + ' -query ' + ii + ' -outfmt "6 qseqid qlen sseqid slen pident qstart qend sstart send length mismatch gaps evalue bitscore " -out ' + oo
    logger.info(bshcmd)
    with open(os.path.join(cutout,"logs","%s.megablast.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.megablast.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'megablast',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(args.nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Megablast error: %s" % err)
            raise
        else:
            logger.info('Blasting reads complete')

########process blast results in R; split by amplicon##############
#if not os.path.exists(Rout):
#    os.makedirs(Rout)
#os.chdir(Rout)
#@mkdir(Rout,os.path.join(Rout,'logs'))            
@transform(mega_blast,suffix('.RE8.RE13.bla'),add_inputs(r"\1_prin_flash.extendedFrags.sed.fasta"), RoutList)#,output_dir=Rout #for x in RoutList 
def proc_blast(input_files,output_files):
    ii1 = input_files[0]
    ii2 = input_files[1]
    oo = output_files[0]
    read_root=re.sub('.RE8.RE13.bla','',os.path.basename(ii1))
    Rcmd=os.path.join(Rpath,'Rscript')+' --no-save --no-restore /data/boehm/group/pipelines/RE_amplicon_seq/v0.1.0/process_blast.R ' + cutout + ' ' + ii1 + ' ' + ii2 + ';sleep 300'
    logger.info(Rcmd)
    with open(os.path.join(cutout,"logs","%s.megablast.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.megablast.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = Rcmd,
                                          job_name          = 'blast_proc',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Postprocessing error: %s" % err)
            raise
        else:
            logger.info('Processing blast results complete')

########cut primers##############
#@mkdir(Rout,os.path.join(Rout,'logs'))   
@transform(proc_blast,suffix('.RE8.e6.hits.fasta'),output=[re.sub('.hits.fasta','.hits.primcut.fasta', x ) for x in RoutList])#,output_dir=cutout
def cut_primers(input_files,output_files):
    ii1 = input_files[0]
    oo1 = output_files[0]
    read_root=re.sub('.hits.fasta','',os.path.basename(ii1))
    cmd_all=[os.path.join(cutpath,'cutadapt') + ' -g file:' + REgFasta + ' -a file:' + REaFasta + ' --minimum-length 100  -n 2 -e 0.2 -o ' + oo + ' ' + ii for REgFasta,REaFasta,oo,ii in zip(gList,aList,output_files,input_files)]
    bshcmd=';'.join(cmd_all)
    logger.info(bshcmd)
    with open(os.path.join(cutout,"logs","%s.primcut.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.primcut.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'primcut',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Primer cutting error: %s" % err)
            raise
        else:
            logger.info('Primer cutting complete')

########collapse identical sequences##############
@mkdir(Fout,os.path.join(Fout,'logs'))   
@transform(cut_primers,suffix('.RE8.e6.hits.primcut.fasta'),output=[re.sub('.hits.fasta','.hits.dedup.fasta', x ) for x in RoutList],output_dir=Fout)#
def dedup_seqs(input_files,output_files):
    ii1 = input_files[0]
    oo1 = output_files[0]
    read_root=re.sub('.hits.primcut.fasta','',os.path.basename(ii1))
    cmd_all=[os.path.join(FXpath,'fastx_collapser') + ' -i ' + ii + ' -o ' + oo for ii,oo in zip(input_files,output_files)]
    bshcmd=';'.join(cmd_all)  
    logger.info(bshcmd)
    with open(os.path.join(Fout,"logs","%s.fastx.out" % read_root),'w+') as stdoutF, open(os.path.join(Fout,"logs","%s.fastx.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'dedup',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Deduplication error: %s" % err)
            raise
        else:
            logger.info('Deduplication complete')

########modify cluster names to include sample name##############
@mkdir(Fout,os.path.join(Fout,'logs'))   
@transform(dedup_seqs,suffix('.RE8.e6.hits.dedup.fasta'),output=[re.sub('.hits.fasta','.hits.dedup.reN.fasta', x ) for x in RoutList],output_dir=Fout)#
def reN_clusters(input_files,output_files):
    ii1 = input_files[0]
    oo1 = output_files[0]
    read_root=re.sub('.hits.dedup.fasta','',os.path.basename(ii1))
    cmd_all=['sed \'s/>/>\'"'+ read_root +'"\'./g\' '+ ii + ' > ' + oo for ii,oo in zip(input_files,output_files)]
    bshcmd=';'.join(cmd_all)  
    logger.info(bshcmd)
    with open(os.path.join(Fout,"logs","%s.reN.out" % read_root),'w+') as stdoutF, open(os.path.join(Fout,"logs","%s.reN.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'reN',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Cluster renaming error: %s" % err)
            raise
        else:
            logger.info('Cluster renaming complete')

#####main

if __name__ == '__main__':
    with open(os.path.join(wdir,"pipelineGraph.png"),'w') as pipeGraph:
        pipeline_printout_graph(stream=pipeGraph,output_format='png',pipeline_name='WGBS',target_tasks=args.target_tasks)
    with open (os.path.join(wdir,"pipelinePrint.txt"),'w') as pipePrint:
        pipeline_printout(verbose_abbreviated_path=0,output_stream=pipePrint,target_tasks=args.target_tasks)    

    pipeline_run(touch_files_only=args.touchOnly,multiprocess=args.bsize,target_tasks=args.target_tasks,forcedtorun_tasks=args.forcedtorun_tasks,logger=logger)
    mySession.exit() 

