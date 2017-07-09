# -*- mode: python -*-

block_cipher = None

mkl_dlls = [
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_vml_def.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_core.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_vml_avx.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_avx512_mic.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_rt.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_vml_avx512.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_avx.dll', ''),
('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_mc3.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_mc.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_vml_cmpt.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_vml_avx2.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_vml_mc.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_msg.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_tbb_thread.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_avx2.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_def.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_vml_mc3.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_avx512.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_sequential.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_vml_mc2.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_intel_thread.dll', ''),
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_vml_avx512_mic.dll', '')
]

a = Analysis(['startgui'],
             pathex=['..\\build\\gui'],
             binaries=mkl_dlls,
             datas=None,
             hiddenimports=['six', 'scipy.integrate'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['matplotlib', 'OpenGL', 'PyQt5.QtOpenGL'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='startgui',
          debug=False,
          strip=False,
          upx=True,
          console=True )
