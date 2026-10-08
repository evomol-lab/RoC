# -*- mode: python ; coding: utf-8 -*-
# PyInstaller build of the web interface:  pyinstaller RoCGUI.spec
# The executable starts the server and opens the browser; calculations run the same
# executable with --run-roc (see RoCGUI.main). ClustENM (OpenMM) is not bundled.
from PyInstaller.utils.hooks import collect_all

datas = [('templates', 'templates'), ('static', 'static'), ('DOCUMENTATION.md', '.')]
binaries = []
hiddenimports = ['roc']
for package in ('prody', 'MDAnalysis'):
    d, b, h = collect_all(package)
    datas += d
    binaries += b
    hiddenimports += h

a = Analysis(
    ['RoCGUI.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['tkinter', 'openmm', 'pdbfixer'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='RoomOfConformations',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
