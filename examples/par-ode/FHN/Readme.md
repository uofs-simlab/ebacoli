# FHN directory summary

## Simulation code

Driver file: [driver-FHN_monodomain-frameoutput.f95](./driver-FHN_monodomain-frameoutput.f95)
- sets parameters for modified FHN model
- sets `split` where the IC reaches equilibrium
- sets solver parameters (tolerance, order, etc)

System file: [FHN_monodomain.f](/FHN_monodomain.f)
- specifies model equations and Jacobian
- specifies IC (uses `split`)
- specifies BC (Neumann)

Once compiled, simply run the executable as `./driver-FHN_monodomain-frameoutput` to generate the files `Bsplines??????`.

## Post-processing code

Assumes `PYTHONPATH` env variable is set properly to find ebacoli python files.

[generate_sample_spatial_points.py](./generate_sample_spatial_points.py)
- generates the locations of sample points for ref solution generation
- call as `./generate_sample_points X_LEFT X_RIGHT N_POINTS > spatial_points.txt`

[sample_bsplines_at_discrete_points.py](./sample_bsplines_at_discrete_points.py)
- evaluate the set of Bspline files at desired spatial points
- call as `./sample_bsplines_at_discrete_points.py spatial_points.txt Bsplines?????? > ref_sol.txt`

[plot_all_bsplines.py](./plot_all_bsplines.py)
- generates plots of a set of Bspline files
- note: hard-coded locations for text overlays
  - i.e. they need to be changed if simulation domain changes
- call as `plot_all_bsplines.py png Bsplines??????`
