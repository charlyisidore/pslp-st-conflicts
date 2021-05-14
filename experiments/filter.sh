#!/bin/bash
# filter instances by number of items: bash filter.sh <min> <max>

if test "$#" -lt 2
then
    echo "Usage: $0 <min> <max>"
    exit 0
fi

find instances/ -type "f" -name "*.json" -exec sh -c "n=\`jq .n_items \$0\`; test \\( \$n -ge $1 \\) -a \\( \$n -le $2 \\)" {} \; -and -print