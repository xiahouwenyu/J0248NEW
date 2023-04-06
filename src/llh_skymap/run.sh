# dir=/home/lhaaso/caowy/0_Tools/sig_3ml
exe=./pixfitting_spec.py

# map=../../data/residual_all.root
# map=../../data/J0248.root
# map=../../data/gcd_new.root
map=../../data/No_diffuse.root

response=../../data/WCDA_DR_psf.root

python3 ${exe} -m ${map} -r ${response} --actBin 0 --name no_fiffuse_test