
# Local ancestry

The formatting of mixed data as well as the reference data for the inference of local alcestrality are detailed in the following steps.

```bash
# 01 subsets de referencias e ISA GRCh37

plink1.9
  --bfile merge_isa1000g
  --keep afr95.iid
  --make-bed
  --out afr/afr

plink1.9
  --bfile merge_isa1000g
  --keep amr95.iid
  --make-bed
  --out amr/amr

plink1.9
  --bfile merge_isa1000g
  --keep eas95.iid
  --make-bed
  --out eas/eas

plink1.9
  --bfile merge_isa1000g
  --keep eur95.iid
  --make-bed
  --out eur/eur

plink1.9
  --bfile merge_isa1000g
  --keep isa.grm_0.125.id
  --make-bed
  --out isa/isa


# 02 split binaries

for k in $(seq 1 22)
do 

  plink1.9 \
  --bfile afr/afr \
  --chr $k \
  --make-bed \
  --out afr/afr.chr${k}

  plink1.9 \
  --bfile amr/amr \
  --chr $k \
  --make-bed \
  --out amr/amr.chr${k}

  plink1.9 \
  --bfile eas/eas \
  --chr $k \
  --make-bed \
  --out eas/eas.chr${k}

  plink1.9 
  --bfile isa/isa \
  --chr $k \
  --make-bed \
  --out isa/isa.chr${k} \

done


# 03 preparing phased inputs for o rfmix

for k in $(seq 22)
do 
  shapeit -B afr/afr.chr${k} --force -M 1000GP_Phase3/genetic_map_chr${k}_combined_b37.txt   -O afr.phased/afr.phased.chr${k} -T 15

  shapeit -B amr/amr.chr${k} --force -M 1000GP_Phase3/genetic_map_chr${k}_combined_b37.txt -O amr.phased/amr.phased.chr${k} -T 15

  shapeit -B eas/eas.chr${k} --force -M 1000GP_Phase3/genetic_map_chr${k}_combined_b37.txt -O eas.phased/eas.phased.chr${k} -T 15

  shapeit -B eur/eur.chr${k} --force -M 1000GP_Phase3/genetic_map_chr${k}_combined_b37.txt -O eur.phased/eur.phased.chr${k} -T 15

  shapeit -B isa/isa.chr${k} --force -M 1000GP_Phase3/genetic_map_chr${k}_combined_b37.txt -O isa.phased/isa.phased.chr${k} -T 15 \

done


# 04 Local Ancestry inference

conda activate py2

for k in $(seq 1 22)
do

  python rfmix/ancestry_pipeline-master/shapeit2rfmix.py \
    --shapeit_hap_ref afr.phased/afr.phased.chr${k}.haps,amr.phased/amr.phased.chr${k}.haps,eas.phased/eas.phased.chr${k}.haps,eur.phased/eur.phased.chr${k}.haps \
    --shapeit_hap_admixed  isa.phased/isa.phased.chr${k}.haps \
    --shapeit_sample_ref afr.phased/afr.phased.chr${k}.sample,amr.phased/amr.phased.chr${k}.sample,eas.phased/eas.phased.chr${k}.sample,eur.phased/eur.phased.chr${k}.sample \
    --shapeit_sample_admixed isa.phased/isa.phased.chr${k}.sample \
    --ref_keep ref/all.ref \
    --admixed_keep isa/isa.iid \
    --chr ${k} \
    --genetic_map 1000GP_Phase3/genetic_map_chr${k}_combined_b37.txt \
    --out out.shapeit2rfmix/isa.rfmix

  python rfmix/RFMix_v1.5.4/RunRFMix.py \
  -e 1 \
  -w 0.2 \
  -n 5 \
  --num-threads 15 \
  --forward-backward  PopPhased \
  out.shapeit2rfmix/isa.rfmix_chr${k}.alleles \
  out.shapeit2rfmix/isa.rfmix.classes \
  out.shapeit2rfmix/isa.rfmix_chr${k}.snp_locations \
  -o outs.rfmix/isa.rfmix_chr${k}

done
```