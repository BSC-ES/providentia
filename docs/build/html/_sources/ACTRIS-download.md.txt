# ACTRIS

Providentia's download mode supports downloading ACTRIS in-situ measurements from the production server at NILU:

[https://prod-actris-md.nilu.no/index.html](https://prod-actris-md.nilu.no/index.html)

## Available datasets

Direct access to the files is found on [Thredds](https://thredds.nilu.no/thredds/catalog.html). Providentia only downloads the files that are contained within the `EBAS_DOI` folder.

### Available variables

These are the available species:

* `aerosol particle light absorption coefficient`
* `precipitation depth`
* `aerosol particle light hemispheric backscatter coefficient`
* `aerosol particle light scattering coefficient`
* `aerosol particle optical depth`
* `precipitation aluminium mass concentration`
* `precipitation arsenic mass concentration`
* `precipitation calcium mass concentration`
* `precipitation cadmium mass concentration`
* `precipitation chloride mass concentration`
* `precipitation cobalt mass concentration`
* `precipitation chromium mass concentration`
* `precipitation copper mass concentration`
* `precipitation iron mass concentration`
* `precipitation mercury mass concentration`
* `precipitation potassium mass concentration`
* `precipitation magnesium mass concentration`
* `precipitation manganese mass concentration`
* `precipitation methanesulfonic acid mass concentration`
* `precipitation sodium mass concentration`
* `precipitation ammonium mass concentration`
* `precipitation nickel mass concentration`
* `precipitation nitrate mass concentration`
* `precipitation lead mass concentration`
* `precipitation selenium mass concentration`
* `precipitation sulphate mass concentration`
* `precipitation vanadium mass concentration`
* `precipitation zinc mass concentration`
* `air pressure`
* `aerosol particle aluminium mass concentration`
* `acetaldehyde mass concentration`
* `aerosol particle arsenic mass concentration`
* `gas and particle benzo(a)pyrene mass concentration`
* `aerosol particle benzo(a)pyrene mass concentration`
* `aerosol particle equivalent black carbon mass concentration`
* `aerosol particle total carbon mass concentration`
* `monoterpenes amount fraction`
* `ethene amount fraction`
* `ethane amount fraction`
* `propene amount fraction`
* `propane amount fraction`
* `1,3-butadiene amount fraction`
* `1-butene amount fraction`
* `n-pentane amount fraction`
* `n-hexane amount fraction`
* `benzene amount fraction`
* `toluene`
* `1,2,4-trimethylbenzene amount fraction`
* `aerosol particle calcium mass concentration`
* `aerosol particle cadmium mass concentration`
* `methane amount fraction`
* `aerosol particle chloride mass concentration`
* `carbon monoxide amount fraction`
* `aerosol particle cobalt mass concentration`
* `aerosol particle chromium mass concentration`
* `aerosol particle copper mass concentration`
* `aerosol particle mineral dust mass concentration`
* `aerosol particle elemental carbon mass concentration`
* `ethanol mass concentration`
* `aerosol particle iron mass concentration`
* `glyoxal mass concentration`
* `formaldehyde mass concentration`
* `hydrochloric acid mass concentration`
* `aerosol particle mercury mass concentration`
* `gaseous elemental mercury mass concentration`
* `gaseous mercury oxides mass concentration`
* `gaseous elemental mercury and mercury oxides mass concentration`
* `nitric acid mass concentration`
* `isoprene amount fraction`
* `aerosol particle potassium mass concentration`
* `aerosol particle magnesium mass concentration`
* `aerosol particle manganese mass concentration`
* `m/p-xylenes amount fraction`
* `aerosol particle methanesulfonic acid mass concentration`
* `m-xylene amount fraction`
* `aerosol particle sodium mass concentration`
* `ammonia mass concentration`
* `aerosol particle ammonium mass concentration`
* `aerosol particle ammonium nitrate mass concentration`
* `aerosol particle nickel mass concentration`
* `nitrogen dioxide mass concentration`
* `aerosol particle nitrate mass concentration`
* `ozone mass concentration`
* `aerosol particle organic carbon mass concentration`
* `o-xylene amount fraction`
* `aerosol particle lead mass concentration`
* `aerosol particle selenium mass concentration`
* `sulfur dioxide mass concentration`
* `aerosol particle sulphate mass concentration`
* `aerosol particle vanadium mass concentration`
* `aerosol particle zinc mass concentration`
* `air temperature`
* `deposition aluminium mass flux`
* `deposition arsenic mass flux`
* `deposition benzo(a)pyrene mass flux`
* `deposition cadmium mass flux`
* `deposition cobalt mass flux`
* `deposition chromium mass flux`
* `deposition copper mass flux`
* `deposition iron mass flux`
* `deposition mercury mass flux`
* `deposition manganese mass flux`
* `deposition lead mass flux`
* `deposition vanadium mass flux`
* `deposition zinc mass flux`

Providentia can understand the ACTRIS vocabulary for some species and an error will be thrown for those who can be mapped to multiple GHOST species (e.g. aerosol particle optical depth). The mapping is provided in the [GHOST–ACTRIS species mapping](#ghostactris-species-mapping) section.

### Available resolutions

These are the available temporal resolutions:

* `hourly` (`time_coverage_resolution`: `P0000-00-00T01:00:00`)
* `3hourly` (`time_coverage_resolution`: `P0000-00-00T03:00:00`)
* `daily` (`time_coverage_resolution`: `P0000-00-01T00:00:00`)
* `monthly` (`time_coverage_resolution`: `P0000-01-00T00:00:00`)

## Example configuration file

```ini
[ACTRIS]
network = actris/actris
species = sconco3
resolution = hourly
start_date = 20180101
end_date = 20190101
observations_data_label = Ozone Mass Concentration
```

## GHOST–ACTRIS species mapping

Here is the mapping between GHOST variable names and the ACTRIS vocabulary. The species that have an asterisk at the end can only be passed to the configuration file using the BSC conventions (name to the left).

```
absco370: aerosol particle light absorption coefficient*
absco467: aerosol particle light absorption coefficient*
absco470: aerosol particle light absorption coefficient*
absco520: aerosol particle light absorption coefficient*
absco522: aerosol particle light absorption coefficient*
absco525: aerosol particle light absorption coefficient*
absco528: aerosol particle light absorption coefficient*
absco530: aerosol particle light absorption coefficient*
absco550: aerosol particle light absorption coefficient*
absco565: aerosol particle light absorption coefficient*
absco590: aerosol particle light absorption coefficient*
absco637: aerosol particle light absorption coefficient*
absco652: aerosol particle light absorption coefficient*
absco660: aerosol particle light absorption coefficient*
absco670: aerosol particle light absorption coefficient*
absco880: aerosol particle light absorption coefficient*
absco950: aerosol particle light absorption coefficient*
acprec: precipitation depth
lbsco450: aerosol particle light hemispheric backscatter coefficient*
lbsco520: aerosol particle light hemispheric backscatter coefficient*
lbsco525: aerosol particle light hemispheric backscatter coefficient*
lbsco530: aerosol particle light hemispheric backscatter coefficient*
lbsco532: aerosol particle light hemispheric backscatter coefficient*
lbsco550: aerosol particle light hemispheric backscatter coefficient*
lbsco635: aerosol particle light hemispheric backscatter coefficient*
lbsco700: aerosol particle light hemispheric backscatter coefficient*
lbsco850: aerosol particle light hemispheric backscatter coefficient*
lsco450: aerosol particle light scattering coefficient*
lsco520: aerosol particle light scattering coefficient*
lsco525: aerosol particle light scattering coefficient*
lsco530: aerosol particle light scattering coefficient*
lsco532: aerosol particle light scattering coefficient*
lsco550: aerosol particle light scattering coefficient*
lsco635: aerosol particle light scattering coefficient*
lsco700: aerosol particle light scattering coefficient*
lsco850: aerosol particle light scattering coefficient*
od1020aero: aerosol particle optical depth*
od380aero: aerosol particle optical depth*
od440aero: aerosol particle optical depth*
od500aero: aerosol particle optical depth*
od550aero: aerosol particle optical depth*
od675aero: aerosol particle optical depth*
od870aero: aerosol particle optical depth*
precal: precipitation aluminium mass concentration
precas: precipitation arsenic mass concentration
precca: precipitation calcium mass concentration
preccd: precipitation cadmium mass concentration
preccl: precipitation chloride mass concentration
preccobalt: precipitation cobalt mass concentration
preccr: precipitation chromium mass concentration
preccu: precipitation copper mass concentration
precfe: precipitation iron mass concentration
prechg: precipitation mercury mass concentration
preck: precipitation potassium mass concentration
precmg: precipitation magnesium mass concentration
precmn: precipitation manganese mass concentration
precmsa: precipitation methanesulfonic acid mass concentration
precna: precipitation sodium mass concentration
precnh4: precipitation ammonium mass concentration
precni: precipitation nickel mass concentration
precno3: precipitation nitrate mass concentration
precpb: precipitation lead mass concentration
precse: precipitation selenium mass concentration
precso4: precipitation sulphate mass concentration
precv: precipitation vanadium mass concentration
preczn: precipitation zinc mass concentration
pshltr: air pressure
sconcal: aerosol particle aluminium mass concentration
sconcald2: acetaldehyde mass concentration
sconcas: aerosol particle arsenic mass concentration
sconcbap: gas and particle benzo(a)pyrene mass concentration
sconcbappm: aerosol particle benzo(a)pyrene mass concentration
sconcbc: aerosol particle equivalent black carbon mass concentration
sconcc: aerosol particle total carbon mass concentration
sconcc10h16: monoterpenes amount fraction
sconcc2h4: ethene amount fraction
sconcc2h6: ethane amount fraction
sconcc3h6: propene amount fraction
sconcc3h8: propane amount fraction
sconcc4h6: 1,3-butadiene amount fraction
sconcc4h8: 1-butene amount fraction
sconcc5h12: n-pentane amount fraction
sconcc6h14: n-hexane amount fraction
sconcc6h6: benzene amount fraction
sconcc7h8: toluene
sconcc9h12: 1,2,4-trimethylbenzene amount fraction
sconcca: aerosol particle calcium mass concentration
sconccd: aerosol particle cadmium mass concentration
sconcch4: methane amount fraction
sconccl: aerosol particle chloride mass concentration
sconcco: carbon monoxide amount fraction
sconccobalt: aerosol particle cobalt mass concentration
sconccr: aerosol particle chromium mass concentration
sconccu: aerosol particle copper mass concentration
sconcdu: aerosol particle mineral dust mass concentration
sconcec: aerosol particle elemental carbon mass concentration
sconcetoh: ethanol mass concentration
sconcfe: aerosol particle iron mass concentration
sconcglyox: glyoxal mass concentration
sconchcho: formaldehyde mass concentration
sconchcl: hydrochloric acid mass concentration
sconchg: aerosol particle mercury mass concentration
sconchggem: gaseous elemental mercury mass concentration
sconchggom: gaseous mercury oxides mass concentration
sconchgtgm: gaseous elemental mercury and mercury oxides mass concentration
sconchno3: nitric acid mass concentration
sconcisop: isoprene amount fraction
sconck: aerosol particle potassium mass concentration
sconcmg: aerosol particle magnesium mass concentration
sconcmn: aerosol particle manganese mass concentration
sconcmpxyl: m/p-xylenes amount fraction
sconcmsa: aerosol particle methanesulfonic acid mass concentration
sconcmxyl: m-xylene amount fraction
sconcna: aerosol particle sodium mass concentration
sconcnh3: ammonia mass concentration
sconcnh4: aerosol particle ammonium mass concentration
sconcnh4no3: aerosol particle ammonium nitrate mass concentration
sconcni: aerosol particle nickel mass concentration
sconcno2: nitrogen dioxide mass concentration
sconcno3: aerosol particle nitrate mass concentration
sconco3: ozone mass concentration
sconcoc: aerosol particle organic carbon mass concentration
sconcoxyl: o-xylene amount fraction
sconcpb: aerosol particle lead mass concentration
sconcse: aerosol particle selenium mass concentration
sconcso2: sulfur dioxide mass concentration
sconcso4: aerosol particle sulphate mass concentration
sconcv: aerosol particle vanadium mass concentration
sconczn: aerosol particle zinc mass concentration
t2: air temperature
wetal: deposition aluminium mass flux
wetas: deposition arsenic mass flux
wetbap: deposition benzo(a)pyrene mass flux
wetcd: deposition cadmium mass flux
wetcobalt: deposition cobalt mass flux
wetcr: deposition chromium mass flux
wetcu: deposition copper mass flux
wetfe: deposition iron mass flux
wethg: deposition mercury mass flux
wetmn: deposition manganese mass flux
wetpb: deposition lead mass flux
wetv: deposition vanadium mass flux
wetzn: deposition zinc mass flux
```