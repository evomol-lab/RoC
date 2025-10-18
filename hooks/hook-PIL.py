from PyInstaller.utils.hooks import collect_submodules

# This hook forces PyInstaller to find all of Pillow's submodules,
# including the hidden '_tkinter_finder' that is needed to display
# images in the Tkinter GUI.
hiddenimports = collect_submodules('PIL')
