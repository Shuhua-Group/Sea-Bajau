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
# build sample index
##################################################

sample_index={
    s:i
    for i,s in enumerate(samples)
}



##################################################
# build population index
##################################################

pop_index={}

gender_index={}



for _,row in meta.iterrows():

    pop=row["Location"]
    iid=row["IID"]
    gender=row["Gender"]


    if iid not in sample_index:

        raise ValueError(
            f"{iid} not found in VCF"
        )


    idx=sample_index[iid]


    pop_index.setdefault(
        pop,
        []
    ).append(
        idx
    )


    gender_index[idx]=gender



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

    """

    AC = 0
    AN = 0


    for allele in gt:


        if allele is None:

            continue


        AN += 1


        if allele > 0:

            AC += 1


    return AC, AN



##################################################
# chromosome ploidy function
##################################################

def expected_ploidy(chrom, gender):

    """
    Return expected allele number
    for one individual.

    Autosomes:
        Male/Female = 2

    chrX:
        Female = 2
        Male = 1

    chrY:
        Male = 1
        Female = 0

    """

    if chrom.startswith("chr"):

        c=chrom.replace(
            "chr",
            ""
        )

    else:

        c=chrom



    if c in ["X","x"]:

        if gender=="Female":

            return 2

        else:

            return 1



    elif c in ["Y","y"]:

        if gender=="Male":

            return 1

        else:

            return 0



    else:

        return 2



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


        AC,AN=gt_to_ac_an(
            record.samples[sample]["GT"]
        )


        if AN==0:

            dosage[i]=-1


        else:

            dosage[i]=AC



    ##############################################
    # population calculation
    ##############################################

    for pop,indices in pop_index.items():


        gt=dosage[indices]


        called=gt>=0


        called_n=called.sum()


        N=len(indices)



        ##########################################
        # missing rate
        ##########################################

        missing_rate=1-called_n/N



        ##########################################
        # calculate AC and AN
        ##########################################

        if called_n>0:


            AC=int(
                gt[called].sum()
            )


            AN=0


            for sample_idx in np.array(indices)[called]:


                gender=gender_index[sample_idx]


                AN += expected_ploidy(
                    record.chrom,
                    gender
                )



            if AN>0:

                AF=AC/AN

            else:

                AF=np.nan



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
