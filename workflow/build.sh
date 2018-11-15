# convert the notebooks into actual python modules
jupyter nbconvert --to script workflow.ipynb
jupyter nbconvert --to script experiments.ipynb
jupyter nbconvert --to script fasta.ipynb
jupyter nbconvert --to script pdbfiles.ipynb
jupyter nbconvert --to script ../SeqD-HBM/SeqDHBM.ipynb

# create a package for deployment
mkdir seqdhbm
mv *.py seqdhbm
cp ../SeqD-HBM/*.py seqdhbm
cp ../SeqD-HBM/wget.exe seqdhbm
echo "__all__ = ['pdbfiles', 'SeqDHBM', 'workflow', 'fasta']" >> seqdhbm/__init__.py

echo "Package created in ./seqdhbm"
