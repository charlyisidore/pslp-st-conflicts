#!/bin/bash
# Instance converter

indir="instances_oelschlagel"
outdir="oelschlagel"

find "${indir}" -name '*.dat' -printf '%P\0' |
    while IFS= read -r -d '' input
    do
        dir="$(dirname "${input}")"
        file="$(basename "${input}")"
        output="${dir}/${file//.dat}.json"
        in="${indir}/${input}"
        out="${outdir}/${output}"
        echo "${in} -> ${out}"
        mkdir -p "${outdir}/${dir}"
        python3 convert_oelschlagel.py -i "${in}" -o "${out}"
    done