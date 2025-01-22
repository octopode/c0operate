# c0operate

### Tools for comparative biophysical interpretation of lipidomic data

`WORK IN PROGRESS; STAY TUNED - 20250112`

The code here implements a simplistic but conservative model for the effect of phospholipid composition on membrane monolayer curvature and thus phase behavior. 

Currently, that model is the same as in [this paper](https://www.science.org/stoken/author-tokens/ST-1955/full), which has its own [repo](https://github.com/octopode/homeocurvature-2024). New lipid classes (ether lipids! cardiolipin!) are coming soon in a more refined model.

I have also thrown in some simple math from [seelipids](https://github.com/octopode/seelipids) in to calculate double bond and chain length indices (DBI and CLI), measures that relate to membrane fluidity. Perhaps someday I will try for a more comprehensive fluidity index that incorporates headgroup and backbone structure.

## Model assumptions

_All models are wrong..._
This one is no exception, but I hope you will find it useful!

0. **Your lipids are correctly annotated.** 

1. **Curvature is additive.** [This is simply wrong on a nano-scale](), but it can be close enough on average. MD simulations can help test the assumption for any given purpose.

2. **Sterol effects are constant/negligible.** This is only a good assumption in comparative systems with very low fractions of total sterol or low variation in sterol concentration. The model does not take sterols into account because their effect on curvature is famously 

3. **Double bonds in hydrocarbon chains increase curvature.** This appears to be true for all _cis_-unsaturations.

4. **Additional carbons in hydrocarbon chains increase curvature.**

## How _not_ to use

The phospholipid curvature index is a comparative index, not a biophysical quantity. Interpret it as such!

### Debugging tips

0. It is recommended to keep an unmodified version of the repo so you can verify that example plots can be reproduced on your system and so that errors caused by modifying the example scripts can be isolated.

1. Lipidomics services tend to make minor changes in datafile formatting over time, which can cause parsing errors. Some parameters to look at and tweak if your raw data fails to parse:

- Argument `skip =`: the number of rows to skip at the beginning of a datasheet. These often contain unformatted notes and comments.

- `cols_meta_ltyp`: list of the metadata columns provided by LipoType. If they remove or change any of these names, `parse_lipotype.R` will try to handle the metadata as mole percentages and there will be trouble. Make sure it captures all the metadata columns in _your_ files!

> Feel free to open an issue if you cannot get your data into tidy format.

## License

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.