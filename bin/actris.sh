#!/bin/bash

species_list=(
    # "absco370"
    # "acprec"
    # "lbsco450" 
    # "lsco450"
    # "od500aero"
    # "precal" 
    # "precas" 
    # "precca" 
    # "preccd" 
    # "preccl" 
    # "preccobalt" 
    # "preccr" 
    # "preccu" 
    # "precfe" 
    # "prechg"
    # "preck" 
    # "precmg"
    # "precmn" 
    # "precmsa" # No data
    # "precna"
    # "precnh4" 
    # "precni" 
    # "precno3" 
    # "precpb" 
    # "precse" # No data
    # "precso4" 
    # "precv" 
    # "preczn"
    # "pshltr" 
    # "sconcal"
    # "sconcas"
    # "sconcbap" 
    # "sconcbappm" 
    # "sconcbc" 
    # "sconcc" # Artifacts
    # "sconcc10h16" # No data
    # "sconcc2h4"
    # "sconcc2h6" 
    # "sconcc3h6" 
    # "sconcc3h8" 
    # "sconcc4h6" 
    # "sconcc4h8" 
    # "sconcc5h12" 
    # "sconcc6h14" 
    # "sconcc6h6" 
    # "sconcc7h8" # No data
    # "sconcc9h12"
    # "sconcca" 
    # "sconccd" 
    # "sconcch4" 
    # "sconccl" 
    # "sconcco" # Tower inlet height
    # "sconccobalt" 
    # "sconccr" 
    # "sconccu" 
    # "sconcdu" 
    # "sconcec" # Artifacts
    # "sconcfe" 
    # "sconchcho" 
    # "sconchcl" 
    # "sconchg" 
    # "sconchgtgm" 
    # "sconchno3" 
    # "sconcisop"
    # "sconck" 
    # "sconcmg" 
    # "sconcmn" 
    # "sconcmpxyl" # Error with API
    # "sconcmsa" # No data
    # "sconcmxyl" # No data
    # "sconcna" 
    # "sconcnh3" 
    # "sconcnh4" 
    # "sconcnh4no3"
    # "sconcni" 
    # "sconcno2" 
    # "sconcno3" 
    # "sconco3" 
    # "sconcoc" # Artifacts
    # "sconcoxyl" 
    # "sconcpb" 
    # "sconcse" 
    # "sconcso2" 
    # "sconcso4" 
    # "sconcv"
    # "sconczn" 
    # "t2" 
     "wetal" # Review station selection -> Same station with different months, but one gets selected, no plots
    # "wetas" 
    # "wetbap" 
    # "wetcd" 
    # "wetcobalt" # Review station selection -> Same station with different months, but one gets selected, no plots
    # "wetcr"
    # "wetcu" 
    # "wetfe" # Review
    # "wethg" # Review
    # "wetmn" # Review
    # "wetpb" 
    # "wetv" 
    # "wetzn"
)

cd /home/avilanov/software/providentia

# Run for each species
for species in "${species_list[@]}"; do
    echo "Running for species: $species"
    ./bin/providentia --conf=test_actris.conf --species="\"$species\"" --dl
    ./bin/providentia --conf=test_actris.conf --species="\"$species\"" --report
done
