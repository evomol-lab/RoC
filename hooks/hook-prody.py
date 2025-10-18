# hooks/hook-prody.py
from PyInstaller.utils.hooks import collect_data_files

# Collect data for ProDy and its dependency, Biopython
datas = collect_data_files('prody')
datas += collect_data_files('Bio')