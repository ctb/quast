#!/bin/bash

python_interpreter=/Users/t/dev/quast/env/bin/coverage

sh ../clean.sh

for f in test_*.py; do
    $python_interpreter run "$f";
    if [ ! $? -eq 0 ]; then exit 1; fi;
    echo;
    echo;
    done
