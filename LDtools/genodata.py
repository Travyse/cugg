# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/00_Genodata.ipynb (unless otherwise specified).

__all__ = ['read_bgen', 'read_bim', 'bgen2dask', 'pybgen_region', 'extract_bed', 'Genodata', 'write_plink', 'write_fam',
           'write_bim', 'write_bed']

# Cell
import numpy as np
import pandas as pd
import dask.array as da
from bgen_reader import open_bgen
from pandas_plink import read_plink
from pandas_plink.bed_reader import lib, ffi
try:
    from pybgen.parallel import ParallelPyBGEN as PyBGEN
except:
    print('Can not import ParallelPyBGEN. import PyBGEN instead')
    from pybgen import PyBGEN

# Cell
from math import floor
from pathlib import Path
from typing import Optional, Union
from tqdm import tqdm
from numpy import ascontiguousarray, empty, float32, float64, nan_to_num, uint8, uint64, arange, full
from xarray import DataArray
from pandas import DataFrame, array

# Cell
def read_bgen(file, sample_file=None,pybgen=True):
    '''the function to read genotype data'''
    if pybgen:
        bg = PyBGEN(file,probs_only=True)
        bim = []
        for i,t in enumerate(bg.iter_variant_info()):
            bim.append([int(t.chrom),t.name,0.0,t.pos,t.a1,t.a2,i])
        bim = pd.DataFrame(bim,columns=['chrom','snp','cm','pos','a0','a1','i'])
        bim.snp = 'chr'+bim[['chrom','pos','a0','a1']].astype(str).agg(':'.join, axis=1)
    else:
        bg = open_bgen(file,verbose=False)
        snp,aa0,aa1 = [],[],[]
        for c,p,alleles in zip(bg.chromosomes,bg.positions,bg.allele_ids):
            a0,a1 = alleles.split(',')
            aa0.append(a0)
            aa1.append(a1)
            snp.append(':'.join(['chr'+str(int(c)),str(p),a0,a1]))  # '05' first change to int, then change to str
        bim = pd.DataFrame({'chrom':bg.chromosomes.astype(int),'snp':snp,'pos':bg.positions,'a0':aa0,'a1':aa1})
    if sample_file is None:
        fam = None
    else:
        fam = pd.read_csv(sample_file, header=0, delim_whitespace=True, quotechar='"',skiprows=1)
        fam.columns = ['fid','iid','missing','sex'] #Fix me
        fam = fam
    return bim,fam,bg

# Cell
def read_bim(fn):
    header = ["chrom", "snp", "cm","pos","a0", "a1"]
    df = pd.read_csv(fn,delim_whitespace=True,header=None,names=header,compression=None,engine="c",iterator=False)
    df["i"] = range(df.shape[0])
    return df

# Cell
def bgen2dask(bgen,index,step=500):
    '''The function to covert bgen to dask array'''
    genos = []
    n = len(index)
    for i in range(0,n,step):
        onecode_geno = bgen.read(index[i:min(n,i+step)])  #samples x variants
        geno = onecode_geno.argmax(axis=2).astype(np.int8)
        genos.append(da.from_array(geno))
    return(da.concatenate(genos,axis=1).T)

# Cell
def pybgen_region(bgen,region,step=100):
    genos,geno=[],[]
    i = 1
    for _,v in bgen.iter_variants_in_region('0'+str(region[0]) if region[0]<10 else str(region[0]),region[1],region[2]):
        if i % step == 0:
            genos.append(da.from_array(geno))
            geno = []
        geno.append(v.argmax(axis=1).astype(np.int8))
        i += 1
    genos.append(da.from_array(geno))
    return(da.concatenate(genos,axis=0))

# Cell
def extract_bed(geno,idx,row=True,step=500,region=None):  #row = True by variants, row = False by samples
    if isinstance(geno,da.core.Array):
        if row:
            geno = geno[idx,:]
        else:
            geno = geno[:,idx]
    elif isinstance(geno,PyBGEN):
        geno = pybgen_region(geno,region,step)
    else:
        if row:
            #must be numric index
            if type(list(idx)[0]) is bool:
                pd_idx = pd.Series(idx)
                idx = list(pd_idx[pd_idx].index)
            geno = bgen2dask(geno,idx,step)
        else:
            geno = geno.read() # read all variants
            geno = geno[:,idx]
    return geno

# Cell
class Genodata:
    def __init__(self,geno_path,sample_path=None):
        self.bim,self.fam,self.bed = self.read_geno(geno_path,sample_path)

    def __repr__(self):
        return "bim:% s \n fam:% s \n bed:%s" % (self.bim, self.fam, self.bed)

    def read_geno(self,geno_file,sample_file):
        if geno_file.endswith('.bed'):
            bim,fam,bed =  read_plink(geno_file[:-4], verbose=False)
            bim.snp = 'chr'+bim[['chrom','pos','a0','a1']].astype(str).agg(':'.join, axis=1)
        elif geno_file.endswith('.bgen'):
            if sample_file is None:
                sample_file = geno_file.replace('.bgen', '.sample')
            bim,fam,bed = read_bgen(geno_file,sample_file)
        else:
            raise ValueError('Plesae provide the genotype files with PLINK binary format or BGEN format')
        bim.chrom = bim.chrom.astype(int)
        bim.pos = bim.pos.astype(int)
        return bim,fam,bed


    def geno_in_stat(self,stat,notin=False):
        '''The function to find an overlap region between geno data with sumstat'''
        variants = stat.SNP
        self.extractbyvariants(variants,notin)


    def geno_in_unr(self,unr):
        '''The function to find an overlap samples between geno data with unr'''
        samples = unr.IID
        self.extractbysamples(samples)

    def extractbyregion(self,region):
        bim = self.bim
        idx = (bim.chrom == region[0]) & (bim.pos >= region[1]) & (bim.pos <= region[2])
        print('this region',region,'has',sum(idx),'SNPs in Genodata')
        if sum(idx) == 0:
            raise ValueError('The extraction is empty')
        #update bim,bed
        self.extractbyidx(idx,row=True,region=region)

    def extractbyvariants(self,variants,notin=False):  #variants is list or pd.Series
        idx = self.bim.snp.isin(variants)
        if notin:
            idx = idx == False
        if sum(idx) == 0:
            raise ValueError('The extraction is empty')
        #update bim,bed
        self.extractbyidx(idx,row=True)

    def extractbysamples(self,samples,notin=False): #samples is list or pd.Series
        samples = pd.Series(samples,dtype=str)
        idx = self.fam.iid.astype(str).isin(samples)
        if notin:
            idx = idx == False
        if sum(idx) == 0:
            raise ValueError('The extraction is empty')
        #update fam,bed
        self.extractbyidx(idx,row=False)

    def extractbyidx(self,idx,row=True,region=None):
        '''get subset of genodata by index
        if index is numbers, the order of genodata will be sorted by the order of index.
        if row = True, extract by variants. Otherwise, extract by samples.'''
        idx = list(idx)
        self.idx = idx
        if row:
            #update bim
            if type(idx[0]) is bool:
                self.bim = self.bim[idx]
            else:
                self.bim = self.bim.iloc[idx]
        else:
            #update fam
            if type(idx[0]) is bool:
                self.fam = self.fam[idx]
            else:
                self.fam = self.fam.iloc[idx]
        self.bed = extract_bed(self.bed,idx,row,region=region)

    def export_plink(self, bed: Union[str, Path], bim: Optional[Union[str, Path]] = None, fam: Optional[Union[str, Path]] = None,row: str = "variant",verbose: bool = True):
        bed = Path(bed)
        if bim is None:
            bim = bed.with_suffix(".bim")
        if fam is None:
            fam = bed.with_suffix(".fam")
        bim = Path(bim)
        fam = Path(fam)

        write_bed(bed, self.bed, row, verbose)

        _echo("Writing FAM... ", end="", disable=not verbose)
        write_fam(fam, self.fam)
        _echo("done.", disable=not verbose)

        _echo("Writing BIM... ", end="", disable=not verbose)
        write_bim(bim, self.bim)
        _echo("done.", disable=not verbose)



# Cell
def write_plink(
    G,
    bed: Union[str, Path],
    bim: Optional[Union[str, Path]] = None,
    fam: Optional[Union[str, Path]] = None,
    row: str = "variant",
    verbose: bool = True,
):
    """
    Write PLINK 1 binary files into a data array.

    A PLINK 1 binary file set consists of three files:

    - BED: containing the genotype.
    - BIM: containing variant information.
    - FAM: containing sample information.

    The user must provide the genotype (dosage) via a :class:`xarray.DataArray` matrix
    with data type :const:`numpy.float32` or :const:`numpy.float64`. That matrix must
    have two named dimensions: **sample** and **variant**. The only allowed values for
    the genotype are: :const:`0`, :const:`1`, :const:`2`, and :data:`math.nan`.

    Parameters
    ----------
    G
        Genotype with bim, bed, and fam.
    bed
        Path to a BED file.
    bim
        Path to a BIM file.It defaults to :const:`None`, in which case it will try to be
        inferred.
    fam
        Path to a FAM file. It defaults to :const:`None`, in which case it will try to
        be inferred.
    major
        It can be either :const:`"sample"` or :const:`"variant"` (recommended and
        default). Specify the matrix layout on the BED file.
    verbose
        :const:`True` for progress information; :const:`False` otherwise.
    """
    if G.bed.ndim != 2:
        raise ValueError("G has to be bidimensional")

    bed = Path(bed)
    if bim is None:
        bim = bed.with_suffix(".bim")
    if fam is None:
        fam = bed.with_suffix(".fam")
    bim = Path(bim)
    fam = Path(fam)

    write_bed(bed, G.bed, row, verbose)

    _echo("Writing FAM... ", end="", disable=not verbose)
    write_fam(fam, G.fam)
    _echo("done.", disable=not verbose)

    _echo("Writing BIM... ", end="", disable=not verbose)
    write_bim(bim, G.bim)
    _echo("done.", disable=not verbose)


def _echo(msg: str, end: str = "\n", disable: bool = False):
    if not disable:
        print(msg, end=end, flush=True)


def write_fam(filepath: Path, df):
    cols = ["fid", "iid", "father","mother","gender","trait"]
    df = df[cols]
    df.to_csv(
        filepath,
        index=False,
        sep="\t",
        header=False,
        encoding="ascii",
        line_terminator="\n",
    )


def write_bim(filepath: Path, df):
    cols = ["chrom","snp","cm","pos","a0","a1"]
    df = df[cols]
    df.to_csv(
        filepath,
        index=False,
        sep="\t",
        header=False,
        encoding="ascii",
        line_terminator="\n",
    )

# Cell
def write_bed(filepath: Path, G, row='variant', verbose=True):
    """
    Write BED file.
    It assumes that ``X`` is a variant-by-sample matrix.
    """
    if not isinstance(G,da.core.Array):
        G = da.asanyarray(G)

    if row != "variant":
        G = G.T

    row_code = 1 if row == "variant" else 0
    e = lib.write_bed_header(str(filepath).encode(), row_code)
    if e != 0:
        raise RuntimeError(f"Failure while writing BED file {filepath}.")

    nrows = G.shape[0]
    ncols = G.shape[1]

    row_chunk = max(1, floor((1024 * 1024 * 256) / ncols))
    row_chunk = min(row_chunk, nrows)

    G = G.rechunk((row_chunk, ncols))

    row_start = 0
    for chunk in tqdm(G.chunks[0], "Writing BED", disable=not verbose):
        data = G[row_start : row_start + chunk, :].compute()
        if data.dtype not in [float32, float64]:
            msg = "Unsupported data type. "
            msg += "Please, provide a dosage matrix in either "
            msg += "float32 or float64 format."
            raise ValueError(msg)

        _write_bed_chunk(filepath, data)
        row_start += chunk


def _write_bed_chunk(filepath: Path, X):

    base_type = uint8
    base_size = base_type().nbytes
    base_repr = "uint8_t"

    nan_to_num(X, False, 3.0)
    G = ascontiguousarray(X, base_type)
    assert G.flags.aligned

    strides = empty(2, uint64)
    strides[:] = G.strides
    strides //= base_size

    e = lib.write_bed_chunk(
        str(filepath).encode(),
        G.shape[1],
        G.shape[0],
        ffi.cast(f"{base_repr} *", G.ctypes.data),
        ffi.cast("uint64_t *", strides.ctypes.data),
    )
    if e != 0:
        raise RuntimeError(f"Failure while writing BED file {filepath}.")