# https://avigailtaylor.github.io/GeneFEAST/user_guide.html
# # Create a new conda environment with Python 3.12 as GeneFeast require python 3.12
# reticulate::conda_create("py312", python_version = "3.12")
# > reticulate::use_condaenv("py312")
# > reticulate::repl_python()
# Python 3.12.11 (/Users/shenlab/Library/r-miniconda-arm64/envs/py312/bin/python)
# Reticulate 1.43.0 REPL -- A Python interpreter in R.
# Enter 'exit' or 'quit' to exit the REPL and return to R.
# >>> import subprocess
# >>> subprocess.run(["pip", "install", "--upgrade", "setuptools"])
# Requirement already satisfied: setuptools in /Users/shenlab/Library/r-miniconda-arm64/envs/py312/lib/python3.12/site-packages (80.9.0)
# CompletedProcess(args=['pip', 'install', '--upgrade', 'setuptools'], returncode=0)
# >>> subprocess.run(["pip", "install", "genefeast"])
# >>> os.chdir("genefeast_project/")

# Restart R session
# Run again
reticulate::use_condaenv("py312")
reticulate::repl_python()

from genefeast import gf

gf.gf("genefeast_03_setup.yml", "result")
# ValueError: Image size of 76040x1950 pixels is too large. It must be less than 2^16 in each direction.
# This error might be due to the too much gene for certain GO term -> Solution: do not use GOALL when retrieving annotated genes

gf.gf("mmc2_3h_setup.yml", "genefeast_example_result")
