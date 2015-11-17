#!/bin/bash
cd Dependencies

echo 'making sure muscle is executable...'
rm muscle3.8.31
mv muscle3.8.31_i86linux64 muscle3.8.31
chmod +x muscle3.8.31

echo 'making sure usearch is executable...'
rm usearch8.1.1756
mv usearch8.1.1756_i86linux32 usearch8.1.1756
chmod +x usearch8.1.1756

echo 'compiling cutadapt...'
cd cutadapt_source
python setup.py build_ext -i

echo 'done.'