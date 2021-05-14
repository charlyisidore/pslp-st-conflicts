#!/bin/bash
# run the experiments: bash run.sh <instance_list>

if test "$#" -eq 0
then
    echo "Usage: $0 <instance_list>"
    exit 0
fi

# Branch-and-price executable
bp_exe="../build/bp/pslp_bp"

# MIP executable
mip_exe="../build/mip/cplex/pslp_mip_cplex"

# Instance list
instances_list="$1"

# Instance directory
instances_dir="instances"

# Parameters directory
parameters_dir="parameters"

# Results directory
results_dir="results"

# Default settings
default_set_path="default.set"

make_batch() {
    # Load default settings
    echo "set load '${default_set_path}'"
    # Load custom settings
    echo "set load '${1}'"
    # Save non-default settings
    echo "set diffsave '${3}.set'"
    # Load instance
    echo "read '${2}'"
    # Solve
    echo "optimize"
    # Save the statistics table
    echo "write statistics '${3}.stat'"
    # Save the final LP model
    echo "write transproblem '${3}.lp'"
    # Save the final LP solution
    echo "write solution '${3}.sol'"
    # Save the JSON result
    echo "write transproblem '${3}.out.json'"
    # Quit
    echo "quit"
}

cat "${instances_list}" |
    while IFS= read -r instance_path
    do
        instance_dir="$(dirname "$(realpath --relative-to="${instances_dir}" "${instance_path}")")"
        instance_file="$(basename "${instance_path}")"
        instance="${instance_file//.json}"
        echo "====== ${instance_dir}/${instance} ======"
        find "${parameters_dir}" -type "f" -name "*.set" -print | sort |
            while IFS= read -r parameter_path
            do
                parameter_file="$(basename "${parameter_path}")"
                parameter="${parameter_file//.set}"
                result_dir="${results_dir}/${parameter}/${instance_dir}"
                result_path="${result_dir}/${instance}"
                mkdir -p "${result_dir}"
                echo "${instance_dir}/${instance} - ${parameter}"
                output_path="${result_path}.out.json"
                if test ! -f "${output_path}"
                then
                    make_batch "${parameter_path}" "${instance_path}" "${result_path}" > "${result_path}.bat"
                    "${bp_exe}" -b "${result_path}.bat" > "${result_path}.log" 2> "${result_path}.err"
                fi
            done
done

cat "${instances_list}" |
    while IFS= read -r instance_path
    do
        instance_dir="$(dirname "$(realpath --relative-to="${instances_dir}" "${instance_path}")")"
        instance_file="$(basename "${instance_path}")"
        instance="${instance_file//.json}"
        echo "====== ${instance_dir}/${instance} ======"
        for model in "3ind" "bp" "flow"
        do
            # BP model does not support non-transitive instances
            if test "${model}" = "bp"
            then
                if case ${instance_dir} in "oelschlagel/arb/"*) true;; *) false;; esac
                then
                    echo "Skip: ${model}"
                    continue
                fi
            fi
            result_dir="${results_dir}/${model}/${instance_dir}"
            result_path="${result_dir}/${instance}"
            mkdir -p "${result_dir}"
            echo "${instance_dir}/${instance} - ${model}"
            output_path="${result_path}.out.json"
            if test ! -f "${output_path}"
            then
                "${mip_exe}" -m "${model}" -i "${instance_path}" -l "${result_path}.lp" -o "${output_path}" > "${result_path}.log" 2> "${result_path}.err"
            fi
        done
done
