#!/bin/bash
echo

trap "echo ' Aborting...'; exit;" SIGINT

echo "-----> find_domain.py"; parallel python find_domain.py ::: 0 1 || exit 1

echo "-----> plot_domain.py"; python plot_domain.py
