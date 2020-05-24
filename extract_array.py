#!/usr/bin/env python3
import os
import errno
import glob
from argparse import ArgumentParser

import numpy as np
import h5py


def get_info(f5path, rna_name, column='length', strand='+'):
    tombo_key = 'RawGenomeCorrected_000' # actual tombo default
    #tombo_key = "RawGenomeCorrected_params_default" # SGB key
    mapped_start = None
    event_bases = None
    data = None
    f = h5py.File(f5path, "r")
    try:
        mapped_chrom = f['/Analyses/{}/BaseCalled_template/Alignment'.format(tombo_key)].attrs['mapped_chrom']
        #print(mapped_chrom)
        assert mapped_chrom == rna_name
        mapped_strand = f['/Analyses/{}/BaseCalled_template/Alignment'.format(tombo_key)].attrs['mapped_strand']
        assert mapped_strand == strand
        mapped_start = f['/Analyses/{}/BaseCalled_template/Alignment'.format(tombo_key)].attrs['mapped_start']
        events = f['/Analyses/{}/BaseCalled_template/Events'.format(tombo_key)]
        event_bases = events['base']
        event_bases = ''.join(np.array([nuc.decode('utf-8').replace('T','U') for nuc in event_bases]).tolist())
        data = events[column]
    except (KeyError, AssertionError):
        pass
    return mapped_start, event_bases, data


def iter_f5_folder(f5path, rna_name, column=None, quiet=False, sort_reads=True, max_files=None):
    fcount = 0
    mapped_count = 0  
    if sort_reads:
        filenames = []
        fcount = 0
        break_early = False
        for folder in glob.iglob('{}/*/'.format(f5path)):
            #print(folder)
            for filename in glob.iglob('{}/*.fast5'.format(folder)):
                #print(filename)
                f = h5py.File(filename, "r")
                reads = f['/Raw/Reads']
                read_number = reads[list(reads.keys())[0]].attrs['read_number']
                filenames.append((filename, read_number))
                fcount += 1
                if not quiet and fcount % 100 == 0:
                    print("located {} fast5 files so far".format(fcount))
                if max_files is not None and fcount >= max_files:
                    break_early = True
                    break
            if break_early:
                break
        if not quiet:
            print("located {} total fast5 files".format(len(filenames)))
        filenames = sorted(filenames, key=lambda x: x[1])
        filenames = [tup[0] for tup in filenames]
        if not quiet:
            print("sorted fast5 filenames by read number in run")
    else:
        filenames = []
        break_early = False
        for folder in glob.iglob('{}/*/'.format(f5path)):
            for filename in glob.iglob('{}/*.fast5'.format(folder)):
                filenames.append(filename)
                if max_files is not None and len(filenames) == max_files:
                    break_early = True
                    break
            if break_early:
                break

    fcount = 0
    mapped_count = 0
    for filename in filenames:
        mapped_start, event_bases, data = get_info(filename, rna_name, column=column)
        fcount += 1
        if data is not None:
            mapped_count += 1
        if not quiet:
            if fcount % 100 == 0:
                print("    parsed {} fast5 files".format(fcount))
                print("    {} / {} reads mapped".format(mapped_count, fcount))
        if data is not None:
            yield filename, mapped_start, event_bases, data


def extract_data(inpath, outpath, target_seq, rna_name, data_type="dwell", sort_reads=True, check_seq_aln=True, max_reads=None, max_files=None):
    slen = len(target_seq)
    print("sequence length: {} nucleotides".format(slen))

    if data_type=='dwell':
        column='length'
    elif data_type=='current':
        column='norm_mean'
    elif data_type=='stdev':
        column='norm_stdev'

    # first do one pass to determine required array size (total read count with mapped bases)
    print("determining array size . . .")
    read_count = 0
    for _ in iter_f5_folder(inpath, rna_name, column=column, sort_reads=False, max_files=max_files):        
        read_count += 1
        if max_reads is not None and read_count > max_reads:
            read_count = max_reads
            break
    print(". . . {} total reads with mapped basecalls".format(read_count))

    # create memory-mapped array file
    try:
        os.makedirs(os.path.split(outpath)[0])
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
    arr = np.memmap(outpath, dtype=float, mode='w+', shape=(read_count,slen))
    arr[:,:] = np.nan
    print("created empty array with shape {}.".format(arr.shape))
    print("filling array with data . . .")

    row_index = 0
    for fname, mapped_start, event_bases, data in iter_f5_folder(inpath, rna_name, column=column, sort_reads=sort_reads, quiet=True, max_files=max_files):
        left = mapped_start
        right = mapped_start+len(data)
        
        if check_seq_aln:
            subseq = target_seq[left:right]
            try:
                assert(subseq == event_bases)
            except AssertionError:
                print(fname)
                print(mapped_start)
                print(target_seq)
                print(event_bases)
                print(subseq)
                print("exiting early, sequence mismatch")
                exit()
        try:
            arr[row_index, left:right] = data
            #print(dwell_times[:10])
            #print(arr[row_index, left:right][:10])
        except ValueError:
            print(fname)
            print(mapped_start)
            print(target_seq)
            print(event_bases)
            print(subseq)
            print("slen = {}".format(slen))
            print("len(subseq) = {}".format(len(subseq)))
            print("len(event_bases) = {}".format(len(event_bases)))
            print("left:right = {}:{}".format(left,right))
            print("exiting early, tried to access array out of bounds")
            exit()
        
        row_index += 1
        if max_reads is not None and row_index >= max_reads:
            break
        if row_index % 100 == 0:
            print("    processed {}/{} mapped reads".format(row_index, read_count))
    print(". . . array filled.")

ap = ArgumentParser()
ap.add_argument('--rna-name', type=str, required=True)
ap.add_argument('--fa', type=str, required=True)
ap.add_argument('--f5', type=str, required=True)
ap.add_argument('--out', type=str, required=True)
ap.add_argument('--max-reads', type=int, default=None)
ap.add_argument('--max-files', type=int, default=None)
ap.add_argument('--sort-reads', action='store_true')
ap.add_argument('--current', action='store_true')
ap.add_argument('--dwell', action='store_true')
ap.add_argument('--stdev', action='store_true') ##
pa = ap.parse_args()

data_type=None
if pa.current:
    data_type='current'
if pa.dwell:
    data_type='dwell'
if pa.stdev: ##
    data_type='stdev' ##
    
###### pick up here    
    
if pa.current and pa.dwell:
    raise RuntimeError("Can't provide --current and --dwell at the same time")
if pa.current and pa.stdev:
    raise RuntimeError("Can't provide --current and --stdev at the same time")
if data_type is None:
    raise RuntimeError("Use --current, --stdev, or --dwell to specify type of data to extract")

target_seq = ''.join([line.rstrip() for line in open(pa.fa,'r')][1:]).replace('T','U')
extract_data(pa.f5, pa.out, target_seq, pa.rna_name, data_type=data_type, sort_reads=pa.sort_reads, max_reads=pa.max_reads, max_files=pa.max_files)
print("Successfully extracted data for reads in {}".format(pa.f5))