#!/usr/bin/env bash

# Parameters used 60_48_8_36
# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/PGA_60_48_8_36"
cd "$wd" || exit 1
cut -d "," -f2,3 PGA_principal_color.csv|uniq


#1dff00,1 # Entry point
#dc0000,5 # seg2.1
#0067dd,2 # intergenic
#00aa8d,4 # seg1
#ffc500,3 # seg2.1 (PGA5)
#870098,0 # Output
#00d400,6 # The remaining fragments are short, Exclude from this step
#0000c1,7
#dcf400,8
#009fcb,9
#ff4500,10
#00a200,11

cut -d "#" -f1 PGA.ctg.csv|uniq > ./haps.list
# $ wc -l < ./haps.list
# 367

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r hap <&9 || [[ -n "$hap" ]]; do
    [[ -z "$hap" ]] && continue
    grep "^${hap}#" PGA.ctg.csv |cut -d "," -f3|uniq > temp.txt
    for node in {6..11};do
        sed -i "/${node}/d" temp.txt
        sed -i 's/\r$//' temp.txt
    done
        mapfile -t lines < temp.txt
        if [[ "${lines[0]}" == "1" && "${lines[-1]}" == "0" ]]; then
            result=("1")  #  Initialize the results list ,  Contains the start 1
            i=1
            while [[ "$i" -lt $((${#lines[@]} - 1)) ]]; do
                if [[ "${lines[$i]}" == "4" && "${lines[$((i+1))]}" == "5" && "${lines[$((i+2))]}" == "2" ]]; then
                    result+=("PGA4")
                    i=$((i + 3))
                elif [[ "${lines[$i]}" == "5" && "${lines[$((i+1))]}" == "2" ]]; then
                    result+=("PGA3")
                    i=$((i + 2))
                elif [[ "${lines[$i]}" == "4" && "${lines[$((i+1))]}" == "3" ]]; then
                    result+=("PGA5")
                    i=$((i + 2))
                else
                    i=$((i + 1))  #  Skip records that do not match
                fi
            done
            result+=("0")  #  Append at the end 0
            # Write one comma-separated row
            GT=`echo "${result[*]}" | sed 's/ /,/g'`
            echo -e "${hap}\t${GT}" >> PGA.GT.out
        else
            echo -e "${hap}\tDoes not start with 1 and end with 0 — skipping" >> PGA.GT.log
        fi
done 9< "./haps.list"



# Inspect the beginning and end
if [[ "${lines[0]}" == "1" && "${lines[-1]}" == "0" ]]; then
    result=("1")  #  Initialize the results list ,  Contains the start 1
    i=1
    while [[ "$i" -lt $((${#lines[@]} - 1)) ]]; do
        if [[ "${lines[$i]}" == "4" && "${lines[$((i+1))]}" == "5" && "${lines[$((i+2))]}" == "2" ]]; then
            result+=("PGA4")
            i=$((i + 3))
        elif [[ "${lines[$i]}" == "5" && "${lines[$((i+1))]}" == "2" ]]; then
            result+=("PGA3")
            i=$((i + 2))
        elif [[ "${lines[$i]}" == "4" && "${lines[$((i+1))]}" == "3" ]]; then
            result+=("PGA5")
            i=$((i + 2))
        else
            i=$((i + 1))  #  Skip records that do not match
        fi
    done
    result+=("0")  #  Append at the end 0

    # Write one comma-separated row
    GT=`echo "${result[*]}" | sed 's/ /,/g'`
    echo -e "${}\t${GT}" >> PGA.GT.out
else
    echo "Does not start with 1 and end with 0 — skipping"
fi

# Then

sed "s/_hap1//g" PGA.GT.out | sed "s/_hap2//g" > PGA.GT.final.out
sed "s/.hap1//g" PGA.GT.final.out |sed "s/.hap2//g" > PGA.GT.final2.out
rm PGA.GT.final.out
mv PGA.GT.final2.out PGA.GT.final.out

# Then run python Script
