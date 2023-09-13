#7DS-Forecast-Photoz
<p align="center"><img src="./images/svn_dim_photoz.gif"></p>

<!-- ![7DS photo-z](/images/svn_dim_photoz.gif) -->
<!-- ![7D photoz](https://piskel-imgstore-b.appspot.com/img/40e0dbe8-3c1e-11ee-8bd4-95e893aea127.gif) -->

A prediction study of photometric redshifts in upcoming 7-Dimensional Sky Survey

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8251967.svg)](https://doi.org/10.5281/zenodo.8251967)



**This repo is under sconstuction!**

* [Test page](https://eunheeko.github.io/7DS-Forecast/)
* [Read the Docs](https://7ds-forecast.readthedocs.io/)
* [pip](https://pypi.org/project/svnds/0.0.1/)

## Overview
Using `svnds`, you can generate mock data that would be observed in the future 7DS and probe feasibility of science cases.

### 7-Dimensional Telescope

Chile, El Sauce Observatory (near CTIO/Cerro Pachon)

Altitude: 1700 m

320 clear nights

Median seeing ~ 1.5”


Operated by ObsTech

```geojson
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "id": 1,
      "properties": {
        "ID": 0
      },
      "geometry": {
        "type": "Polygon",
        "coordinates": [
          [
              [-90,35],
              [-90,30],
              [-85,30],
              [-85,35],
              [-90,35]
          ]
        ]
      }
    }
  ]
}
```

* features of 40 medium band filters

### Generate mock data
* `7DS`, `SPHEREx`, `LSST`, `Euclid`, `VIKING`, `PANSTARRS` as of Aug. 26, 2023

### Customize SED models
* to be updated

## Installlation
### pip

`pip install svnds`

### git
`git clone https://github.com/eunheeko/7DS-Forecast.git`

### To be updated

|Input|prior|zstep|Description|
|---|---|---|---|
|7DS Y5|$m_{6250\AA}$|0.01|photoz_7DS_WFS_Y5_eff_adderr|
<!-- |7DS Y5 + VIKING|$m_{6250\AA}$|photoz_sds40_WFS_5yr_eff_VIKING_adderr.coeff| -->

<!-- ## Data Specifiaction
![data_spec](/images/data_specification.png)

## Survey Plan
=======
# Notes


## EAZY
<!-- - read binary files: https://github.com/gbrammer/eazy-photoz/tree/f8b84a20f8e781d1f8244a24cd347a24a40f1558/inputs -->
<!-- - source code: https://eazy-py.readthedocs.io/en/latest/_modules/eazy/photoz.html -->


<!-- ## data/filters
* 7DS, SPHERE: customized
* LSST: 
* [EUCLID](https://ui.adsabs.harvard.edu/abs/2022A%26A...662A..92E/abstract): [data](https://euclid.esac.esa.int/msp/refdata/data/)
* VIKIG:  -->


<!-- ## For large data: -->
<!-- https://docs.github.com/en/repositories/working-with-files/managing-large-files -->
