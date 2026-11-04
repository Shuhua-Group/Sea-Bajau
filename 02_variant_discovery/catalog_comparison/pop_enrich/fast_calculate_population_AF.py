#!/usr/bin/env python3
"""Process variant calls and summarize variant catalogs."""



import argparse
import pandas as pd
import numpy as np
import pysam



parser=argparse.ArgumentParser()


parser.add_argument(
    "--vcf",
    required=True
)

parser.add_argument(
    "--meta",
    required=True
)

parser.add_argument(
    "--out",
    required=True
)

parser.add_argument(
    "--missing",
    type=float,
    default=0.8
)


args=parser.parse_args()



##################################################
# metadata
##################################################

meta=pd.read_csv(
    args.meta,
    sep="\t"
)



##################################################
# open VCF
##################################################

vcf=pysam.VariantFile(
    args.vcf
)


samples=list(
    vcf.header.samples
)



##################################################
# build population index
##################################################

sample_index={
    s:i
    for i,s in enumerate(samples)
}



pop_index={}


for _,row in meta.iterrows():

    pop=row["Location"]
    iid=row["IID"]

    if iid not in sample_index:
        raise ValueError(
            f"{iid} not found in VCF"
        )

    pop_index.setdefault(
        pop,
        []
    ).append(
        sample_index[iid]
    )



populations=list(
    pop_index.keys()
)



print(
    f"{len(samples)} samples"
)

print(
    f"{len(populations)} populations"
)


##################################################
# genotype parser
##################################################

def gt_to_ac_an(gt):

    """
    Calculate ALT allele count (AC)
    and called allele number (AN)

    GT examples:
    (0,0) -> AC=0 AN=2
    (0,1) -> AC=1 AN=2
    (1,1) -> AC=2 AN=2
    (None,1) -> AC=1 AN=1
    (None,None) -> AC=0 AN=0

    pysam returns GT as tuple:
    (0,1)
    (1,1)
    (None,1)
    """

    AC = 0
    AN = 0


    for allele in gt:

        # missing allele
        if allele is None:
            continue


        # called allele
        AN += 1


        # ALT allele
        if allele > 0:
            AC += 1


    return AC, AN

##################################################
# main
##################################################

records=[]



for idx,record in enumerate(vcf):


    if idx%100000==0:

        print(
            f"processed {idx} variants"
        )



    ##############################################
    # read all genotype once
    ##############################################

    dosage=np.empty(
        len(samples),
        dtype=np.int8
    )

    for i,sample in enumerate(samples):
        AC, AN = gt_to_ac_an(
            record.samples[sample]["GT"]
        )
        if AN == 0:
            dosage[i] = -1
        else:
            dosage[i] = AC

    ##############################################
    # population calculation
    ##############################################

    for pop,indices in pop_index.items():


        gt=dosage[indices]


        called=gt>=0


        called_n=called.sum()


        N=len(indices)


        missing_rate=1-called_n/N



        if called_n>0:


            AC=gt[called].sum()


            AN=called_n*2


            AF=AC/AN


        else:

            AC=0
            AN=0
            AF=np.nan



        ##########################################
        # filter
        ##########################################

        if missing_rate>args.missing:

            category="NA"
            PASS=False


        else:

            PASS=True


            if AF==0:

                category="Absent"

            elif AF<=0.05:

                category="Rare"

            elif AF<=0.2:

                category="Low-frequency"

            else:

                category="Private-common"



        records.append(
            [
                record.chrom,
                record.pos,
                record.ref,
                ",".join(record.alts),
                pop,
                N,
                int(AC),
                int(AN),
                missing_rate,
                AF,
                category,
                PASS
            ]
        )




##################################################
# save
##################################################

out=pd.DataFrame(
    records,
    columns=[
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "Population",
        "N",
        "AC",
        "AN",
        "missing_rate",
        "AF",
        "AF_category",
        "PASS"
    ]
)


out.to_csv(
    args.out,
    sep="\t",
    index=False
)
