#!/usr/bin/env bash

# Parameters used 60_48_8_36
# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/AMY_80_56_4_64"
cd "$wd" || exit 1
cut -d "," -f2,3 AMY_principal_color.csv|uniq

#870098,0 # Entry point
#00e200,4 # 2B
#0092dd,9 # Following 2B uninformative trailing blue arrows
#d0f800,5 # 2A.1
#ffa500,6 # 2A.2
#0000c5,8 # 2A.3
#dc0000,7 # 2A.4
#009b13,3 # 1
#001cdd,1 # Output 1
#00a3bf,2 # Output 2
### 2A.3 And 2A.4 appear separately, Compose 2Apseudo

#00a97d,10 # The following entries are uninformative
#00c700,11
#84ff00,12
#ffcd00,13
#f60000,14
#2300a6,15
#0078dd,16
#00aaa0,17
#00af00,18
#00fc00,19
#f1e700,20
#ff5100,21




cut -d "#" -f1 AMY.ctg.csv|uniq > ./haps.list
# $ wc -l < ./haps.list
# 359

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r hap <&9 || [[ -n "$hap" ]]; do
    [[ -z "$hap" ]] && continue
    grep "^${hap}#" AMY.ctg.csv |cut -d "," -f3|uniq > temp.txt
    for node in {10..21};do
        sed -i "/${node}/d" temp.txt
        sed -i 's/\r$//' temp.txt
    done
        mapfile -t lines < temp.txt
        if [[ "${lines[0]}" == "0" && "${lines[-1]}" == "2" ]]; then
            result=("1")  #  Initialize the results list ,  Contains the start 1
            i=1
            while [[ "$i" -lt $((${#lines[@]} - 1)) ]]; do
                if [[ "${lines[$i]}" == "5" && "${lines[$((i+1))]}" == "6" && "${lines[$((i+2))]}" == "8" && "${lines[$((i+3))]}" == "7" ]]; then
                    result+=("AMY2A")
                    i=$((i + 4))
                elif [[ "${lines[$i]}" == "7" && "${lines[$((i+1))]}" == "8" && "${lines[$((i+2))]}" == "6" && "${lines[$((i+3))]}" == "5" ]]; then
                    result+=("AMY2A")
                    i=$((i + 4))
                elif [[ "${lines[$i]}" == "8" && "${lines[$((i+1))]}" == "7" ]]; then
                    result+=("AMY2Ap")
                    i=$((i + 2))
                elif [[ "${lines[$i]}" == "7" && "${lines[$((i+1))]}" == "8" ]]; then
                    result+=("AMY2Ap")
                    i=$((i + 2))
                elif [[ "${lines[$i]}" == "4"  ]]; then
                    result+=("AMY2B")
                    i=$((i + 1))
                elif [[ "${lines[$i]}" == "3" ]]; then
                    result+=("AMY1")
                    i=$((i + 1))
                else
                    i=$((i + 1))  #  Skip records that do not match
                fi
            done
            result+=("0")  #  Append at the end 0
            # Write one comma-separated row
            GT=`echo "${result[*]}" | sed 's/ /,/g'`
            echo -e "${hap}\t${GT}" >> AMY.GT.out
        else
            echo -e "${hap}\tDoes not start with 1 and end with 0 — skipping" >> AMY.GT.log
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
