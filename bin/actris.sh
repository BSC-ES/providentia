#!/bin/bash

species_list=(
    "absco370" "absco467" "absco470" "absco520" "absco522" "absco525" "absco528" "absco530" "absco550" "absco565"
    "absco590" "absco637" "absco652" "absco660" "absco670" "absco880" "absco950" "acprec"
    "lbsco450" "lbsco520" "lbsco525" "lbsco530" "lbsco532" "lbsco550" "lbsco635" "lbsco700" "lbsco850"
    "lsco450" "lsco520" "lsco525" "lsco530" "lsco532" "lsco550" "lsco635" "lsco700" "lsco850"
    "od1020aero" "od380aero" "od440aero" "od500aero" "od550aero" "od675aero" "od870aero"
    "precal" "precas" "precca" "preccd" "preccl" "preccobalt" "preccr" "preccu" "precfe" "prechg"
    "preck" "precmg" "precmn" "precmsa" "precna" "precnh4" "precni" "precno3" "precpb" "precse" "precso4" "precv" "preczn"
    "pshltr" "sconcal" "sconcald2" "sconcas" "sconcbap" "sconcbappm" "sconcbc" "sconcc" "sconcc10h16" "sconcc2h4"
    "sconcc2h6" "sconcc3h6" "sconcc3h8" "sconcc4h6" "sconcc4h8" "sconcc5h12" "sconcc6h14" "sconcc6h6" "sconcc7h8" "sconcc9h12"
    "sconcca" "sconccd" "sconcch4" "sconccl" "sconcco" "sconccobalt" "sconccr" "sconccu" "sconcdu" "sconcec" "sconcetoh"
    "sconcfe" "sconcglyox" "sconchcho" "sconchcl" "sconchg" "sconchggem" "sconchggom" "sconchgtgm" "sconchno3" "sconcisop"
    "sconck" "sconcmg" "sconcmn" "sconcmpxyl" "sconcmsa" "sconcmxyl" "sconcna" "sconcnh3" "sconcnh4" "sconcnh4no3"
    "sconcni" "sconcno2" "sconcno3" "sconco3" "sconcoc" "sconcoxyl" "sconcpb" "sconcse" "sconcso2" "sconcso4" "sconcv"
    "sconczn" "t2" "wetal" "wetas" "wetbap" "wetcd" "wetcobalt" "wetcr" "wetcu" "wetfe" "wethg" "wetmn" "wetpb" "wetv" 
    "wetzn"
)

cd /home/avilanov/software/Providentia

# Run for each species
for species in "${species_list[@]}"; do
    echo "Running for species: $species"
    ./bin/providentia --conf=actris.conf --species="$species --dl"
done
