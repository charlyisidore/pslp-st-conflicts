#!/bin/sh
find include/ src/ -type f \( -name '*.cpp' -o -name '*.hpp' \) -exec clang-format --verbose -i '{}' \;