#!/usr/bin/env python3

"""
Calculate population-specific allele frequency
for APC-P1-specific variants.

Input:
    VCF with genotypes
    sample metadata

Output:
    variant x population AF table


Author:
    APC revision analysis
"""


import argparse
import pandas as pd
import numpy as np
import pysam



########################################################
# arguments
########################################################

parser = argparse.ArgumentParser()


parser.add_argument(
    "--vcf",
    required=True,
    help="Input VCF.gz"
)


parser.add_argument(
    "--meta",
    required=True,
    help="Sample information: IID Location"
)


parser.add_argument(
    "--out",
    required=True
)


parser.add_argument(
    "--missing",
    default=0.8,
    type=float,
    help="maximum missing rate"
)


args = parser.parse_args()



########################################################
# read metadata
########################################################


meta = pd.read_csv(
    args.meta,
    sep="\t"
)


required=[
    "IID",
    "Location"
]


for x in required:
    if x not in meta.columns:
        raise ValueError(
            f"Missing column {x}"
        )



sample2pop=dict(
    zip(
        meta.IID,
        meta.Location
    )
)



########################################################
# open VCF
########################################################


vcf=pysam.VariantFile(
    args.vcf
)



samples=list(
    vcf.header.samples
)



# check samples

missing_samples=set(samples)-set(sample2pop)

if len(missing_samples)>0:

    raise ValueError(
        "VCF samples missing in metadata:\n"
        +
        "\n".join(list(missing_samples)[:10])
    )



populations=sorted(
    set(meta.Location)
)



########################################################
# helper
########################################################


def gt_to_ac(gt):

    """
    genotype:
    (0,0)
    (0,1)
    (1,1)

    return ALT dosage
    """

    if gt is None:
        return np.nan


    if None in gt:
        return np.nan


    return sum(
        allele>0
        for allele in gt
    )



########################################################
# main
########################################################


results=[]


for record in vcf:


    variant_info={

        "CHROM":record.chrom,

        "POS":record.pos,

        "REF":record.ref,

        "ALT":",".join(
            record.alts
        )

    }



    # calculate genotype dosage

    dosage={}


    for sample in samples:

        dosage[sample]=gt_to_ac(
            record.samples[sample]["GT"]
        )



    dosage=pd.Series(
        dosage
    )



    ####################################################
    # population calculation
    ####################################################


    for pop in populations:


        pop_samples=[
            s for s in samples
            if sample2pop[s]==pop
        ]


        gt=dosage[
            pop_samples
        ]


        called=gt.notna().sum()


        AN=called*2


        AC=gt.sum(
            skipna=True
        )


        N=len(pop_samples)


        missing_rate=1-called/N



        if AN>0:

            AF=AC/AN

        else:

            AF=np.nan



        ############################################
        # AF category
        ############################################


        if missing_rate > args.missing:

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



        results.append(

            {

            **variant_info,

            "Population":pop,

            "N":N,

            "AC":int(AC),

            "AN":int(AN),

            "missing_rate":missing_rate,

            "AF":AF,

            "AF_category":category,

            "PASS":PASS

            }

        )



########################################################
# save
########################################################


out=pd.DataFrame(results)


out.to_csv(
    args.out,
    sep="\t",
    index=False
)
