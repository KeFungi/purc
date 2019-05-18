#!/bin/bash
cd Dependencies

echo 'making sure muscle is executable...'
chmod +x muscle3.8.31

echo 'making sure usearch is executable...'
chmod +x usearch8.1.1756

echo 'compiling cutadapt...'
cd cutadapt_source
python setup.py build_ext -i

echo 'done.'