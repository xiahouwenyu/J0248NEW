#export PYTHONPATH="/home/lhaaso/gmxiang/lib/pip_lib"
# name=20210305_20220205_pinc_new3.root
# name=20210305_20220630_pinc_new.root
name=gcd_new.root

# pp=(0.3364 0.3608 1 0.6841 0.4544 0.3304 0.2555 0.218 0.1967 0.1952 0.1968)
pp=(0.42 0.32 0.25 0.22 0.18 0.15 0.30 0.27 0.22 0.20 0.17 0.15)
# pp=(0.84 0.64 0.5 0.44 0.36 0.3 0.42 )

for nn in {6..11}
do
python3 convert_root2fits.py -i ${name}  -p ${pp[${nn}]}  -n ${nn}
done
