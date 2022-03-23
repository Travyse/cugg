# In-house Python utilities from Gao Wang's group



### 4 modules
- Genodata
- Sumstats
- Liftover
- LDmatrix

## Install

`pip install cugg`

## How to use

```
lf = Liftover('hg38','hg19')
```

```
vcf ='/home/yh3455/Github/SEQLinkage/MWE/small_sample_ii_coding.vcf.gz'
```

```
lf.vcf_liftover(vcf)
```

```
!which python
```

    /home/yh3455/miniconda3/bin/python


```
region = [5,272741,1213528-900000]
geno_path = 'MWE_region_extraction/ukb23156_c5.merged.filtered.5_272741_1213528.bed'
sumstats_path = 'MWE_region_extraction/090321_UKBB_Hearing_aid_f3393_expandedwhite_6436cases_96601ctrl_PC1_2_f3393.regenie.snp_stats'
pheno_path = None
unr_path = 'MWE_region_extraction/UKB_genotypedatadownloaded083019.090221_sample_variant_qc_final_callrate90.filtered.extracted.white_europeans.filtered.092821_ldprun_unrelated.filtered.prune.txt'
imp_geno_path = 'MWE_region_extraction/ukb_imp_chr5_v3_05_272856_1213643.bgen'
imp_sumstats_path = 'MWE_region_extraction/100521_UKBB_Hearing_aid_f3393_expandedwhite_15601cases_237318ctrl_500k_PC1_PC2_f3393.regenie.snp_stats.gz'
imp_ref = 'hg19'
bgen_sample_path = 'MWE_region_extraction/ukb_imp_chr5_v3_05_272856_1213643.sample'
output_sumstats = 'test.snp_stats.gz'
output_LD = 'test_corr.csv.gz'

#main(region,geno_path,sumstats_path,pheno_path,unr_path,imp_geno_path,imp_sumstats_path,imp_ref,output_sumstats,output_LD)
```

```
    def main(region,geno_path,sumstats_path,pheno_path,unr_path,imp_geno_path,imp_sumstats_path,imp_ref,output_sumstats,output_LD,bgen_sample_path):

        print('1. Preprocess sumstats (regenie format) and extract it from a region')
        if pheno_path is not None:
            # Load phenotype file
            pheno = pd.read_csv(pheno_path, header=0, delim_whitespace=True, quotechar='"')
        if unr_path is not None:
            # Load unrelated sample file
            unr = pd.read_csv(unr_path, header=0, delim_whitespace=True, quotechar='"')  
        # Load the file of summary statistics and standardize it.
        exome_sumstats = Sumstat(sumstats_path)
        exome_geno = Genodata(geno_path,bgen_sample_path)

        print('1.1. Region extraction')
        exome_sumstats.extractbyregion(region)
        exome_geno.extractbyregion(region)
        exome_sumstats.match_ss(exome_geno.bim)
        exome_geno.geno_in_stat(exome_sumstats.ss)

        if imp_geno_path is not None:
            #two genotype data
            imput_sumstats = Sumstat(imp_sumstats_path)
            imput_geno = Genodata(imp_geno_path,bgen_sample_path)

            if imp_ref is None:
                imput_sumstats.extractbyregion(region)
                imput_geno.extractbyregion(region)
                imput_sumstats.match_ss(imput_geno.bim)
                imput_geno.geno_in_stat(imput_sumstats.ss)
            else:
                print('1.2. LiftOver the region')
                hg38toimpref = Liftover('hg38',imp_ref)
                imp_region = hg38toimpref.region_liftover(region)
                imput_sumstats.extractbyregion(imp_region)
                imput_geno.extractbyregion(imp_region)
                imput_sumstats.match_ss(imput_geno.bim)
                imput_geno.geno_in_stat(imput_sumstats.ss)

                print('1.3. Regional SNPs Liftover')
                impreftohg38 = Liftover(imp_ref,'hg38') #oppsite with hg38toimpref
                imput_geno.bim = impreftohg38.bim_liftover(imput_geno.bim)
                imput_sumstats.ss.POS = list(imput_geno.bim.pos)
                imput_sumstats.ss.SNP = list(imput_geno.bim.snp)

            print('1.1.1 Get exome unique sumstats and geno and Combine sumstats')
            snp_match = compare_snps(exome_sumstats.ss,imput_sumstats.ss)
            exome_sumstats.ss = exome_sumstats.ss.loc[snp_match.qidx[snp_match.exact==False].drop_duplicates()] #remove by exact match. can be improved.
            exome_sumstats.extractbyvariants(imput_sumstats.ss.SNP,notin=True)
            exome_geno.geno_in_stat(exome_sumstats.ss)
            sumstats = pd.concat([exome_sumstats.ss,imput_sumstats.ss])
        else:
            #one genotype data
            sumstats = exome_sumstats

        print('2. Remove relative samples')
        if unr_path is not None:
            exome_geno.geno_in_unr(unr)
            if imp_geno_path is not None:
                imput_geno.geno_in_unr(unr)
        else:
            print('Warning:There is no file of relative sample. All sample are included in computing LD matrix')

        if pheno_path is not None:
            print('Warning: This function has been implementd yet.')
            pass #sld and pld

        print('3. Calculate LD matrix')
        if imp_geno_path is None:
            cor_da = geno_corr(exome_geno.bed.T)
        else:
            xx = geno_corr(exome_geno.bed.T)
            yy = geno_corr(imput_geno.bed.T,step=500)

            imput_fam = imput_geno.fam.copy()
            imput_fam.index = list(imput_fam.iid.astype(str))
            imput_fam['i'] = list(range(imput_fam.shape[0]))
            imput_fam_comm = imput_fam.loc[list(exome_geno.fam.iid.astype(str))]
            imput_geno.extractbyidx(list(imput_fam_comm.i),row=False)
            xy = geno_corr(exome_geno.bed.T,imput_geno.bed.T,step=500)
            cor_da = da.concatenate([da.concatenate([xx,xy],axis=1),da.concatenate([xy.T,yy],axis=1)],axis=0)

        print('4. Output sumstats and LD matrix')
        index = list(sumstats.SNP.apply(shorten_id))
        sumstats.SNP = index
        sumstats.index = list(range(sumstats.shape[0]))
        sumstats.to_csv(output_sumstats, sep = "\t", header = True, index = False,compression='gzip')

        corr = cor_da.compute()
        np.fill_diagonal(corr, 1)
        corr = pd.DataFrame(corr, columns=index)
        corr.to_csv(output_LD, sep = "\t", header = True, index = False,compression='gzip')


```

```
main(region,geno_path,sumstats_path,pheno_path,unr_path,imp_geno_path,imp_sumstats_path,imp_ref,output_sumstats,output_LD,bgen_sample_path)
```

    1. Preprocess sumstats (regenie format) and extract it from a region
    1.1. Region extraction
    this region [5, 272741, 313528] has 498 SNPs
    Total SNPs 119 . Flip SNPs 118


    /home/yh3455/miniconda3/lib/python3.8/site-packages/pandas/core/generic.py:5516: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      self[name] = value


    1.2. LiftOver the region
    this region (5, 272856, 313643) has 1736 SNPs
    Total SNPs 385 . Flip SNPs 0
    1.3. Regional SNPs Liftover
    1.1.1 Get exome unique sumstats and geno and Combine sumstats
    keep   exact  flip   reverse  both   complement
    False  False  False  False    False  False         106
    True   False  True   False    False  False          12
           True   False  False    False  False           1
    dtype: int64
    2. Remove relative samples
    3. Calculate LD matrix


    /mnt/mfs/statgen/yin/Github/LDtools/LDtools/ldmatrix.py:29: RuntimeWarning: invalid value encountered in true_divide
      geno_i = (geno_i - np.nanmean(geno_i,axis=0)[None,:])/np.nanstd(geno_i,axis=0)[None,:]
    /mnt/mfs/statgen/yin/Github/LDtools/LDtools/genodata.py:74: PerformanceWarning: Slicing with an out-of-order index is generating 22098 times more chunks
      geno = geno[:,idx]
    /mnt/mfs/statgen/yin/Github/LDtools/LDtools/ldmatrix.py:67: RuntimeWarning: invalid value encountered in true_divide
      geno_i = (geno_i - np.nanmean(geno_i,axis=0)[None,:])/np.nanstd(geno_i,axis=0)[None,:]


    4. Output sumstats and LD matrix

