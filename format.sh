#!/bin/sh
find . -path ./build -prune -o -type f \( -name '*.cpp' -o -name '*.hpp' \) -exec clang-format --verbose -i '{}' \;