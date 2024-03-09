#!/bin/bash

# 检查 pip 是否已安装
if ! command -v pip &> /dev/null; then
    echo "Error: pip is not installed."
    exit 1
fi

# 安装所需的 Python 包
pip install ete3 pandas numpy tqdm matplotlib pypdf4

echo "Installation complete."
