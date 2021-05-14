#!/bin/bash
# analyze the results: bash analyze.sh [<csv_output>]

# Results directory
results_dir="results"

# CSV output file
output_csv="${1}"

header="instance,parameter,status,obj,gap,n_nodes,time"

if test -n "${output_csv}"
then echo "${header}" > "${output_csv}"
else echo "${header}"
fi

ls -1d "${results_dir}"/* |
    while IFS= read -r parameter_path
    do
        parameter="$(basename "${parameter_path}")"
        if test -n "${output_csv}"
        then echo "${parameter}"
        fi
        find "${parameter_path}" -type "f" -name "*.out.json" -print | sort |
            while IFS= read -r instance_out_json
            do
                instance_path="${instance_out_json//.out.json}"
                instance="$(realpath --relative-to="${parameter_path}" "${instance_path}")"
                jq_script=".stats | [.status, (if has(\"best_sol_orig_obj\") then .best_sol_orig_obj else .obj_value end), (if has(\"gap\") then .gap else .mip_relative_gap end), .n_nodes, .total_time] | @csv"
                row="$(jq -r "${jq_script}" "${instance_out_json}")"
                if test -n "${output_csv}"
                then echo "${instance},${parameter},${row}" >> "${output_csv}"
                else echo "${instance},${parameter},${row}"
                fi
            done
    done