# AssetForJwstNirspecFs 
**AssetForJwstNirspecFs** is a Julia package to extract spectrum from JWST/NIRSpec Fixed Slit data. 

This package is a top layer of the more general package [`ASSET`][https://github.com/SlitSpectroscopyBuddies/ASSET], specifficaly for the JWST/NIRSpec Fixed Slit data. It contains calibration methods to obtain the spatial distribution maps required in [`ASSET`][https://github.com/SlitSpectroscopyBuddies/ASSET]. 

 [`ASSET`][https://github.com/SlitSpectroscopyBuddies/ASSET] should be installed prior to the installation of `AssetForJwstNirspecFs`. See the documentation of `ASSET` for its installation.

## Installation

In the package manager:

```julia
pkg> add https://github.com/LaurenceDenneulin/AssetForJwstNirspecFs
```
if you use HTTPS or:

```julia
pkg> add git@github.com:LaurenceDenneulin/AssetForJwstNirspecFs
```
if you use SSH.

