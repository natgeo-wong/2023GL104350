# **<div align="center">2023GL104350</div>**

This repository contains the scripts needed to extract the relevant data from the
[ExploreWTGSpace](https://github.com/natgeo-wong/ExploreWTGSpace) project, and the
notebooks needed to create the figures in the paper submission to GRL with the ID
`2023GL104350`.

This project will also contain the necessary experimental setup files, so that
anyone who wants will be able to run the model to get all the output.  They can
then explore the full range of the project themselves.

The data for this project can be found at [here](https://doi.org/10.7910/DVN/YPXNPG) at the Harvard Dataverse.

## Instructions (to be updated later)

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "2023GL104350"
```
which auto-activate the project and enable local path handling from DrWatson.
