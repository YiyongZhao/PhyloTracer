#!/bin/bash

if ! command -v pip &> /dev/null; then
    echo "Error: pip is not installed."
    exit 1
fi

packages=("ete3" "pandas" "numpy" "tqdm" "matplotlib" "pypdf" "pyqt5" "phyde")
failed_packages=()
for package in "${packages[@]}"; do
    pip install $package || failed_packages+=($package)
done

if [ ${#failed_packages[@]} -eq 0 ]; then
    echo "All packages installed successfully."
else
    echo "Some packages failed to install: ${failed_packages[@]}"
fi
