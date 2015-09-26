#!/bin/bash
cd Dependencies

echo 'checking muscle executable...'
chmod +x muscle3.8.31

echo 'checking usearch executable...'
chmod +x usearch8.1.1756

cd cutadapt_source
python setup.py build_ext -i

cd ..
ln -s cutadapt_source/bin/cutadapt cutadapt