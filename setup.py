from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.
buildOptions = dict(packages = [], excludes = [])

base = 'Console'

executables = [
    Executable('RSA5106A_Control_Software.py', base=base, icon="MESA_favicon.ico")
]

setup(name='RSA5106A',
      version = '0.5.0',
      description = 'RSA5106A Control Software',
      options = dict(build_exe = buildOptions),
      executables = executables)
