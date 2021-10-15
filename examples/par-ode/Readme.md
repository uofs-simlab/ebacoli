Brief description of extended examples found in this directory.

[./CDM](./CDM):

-   Ecological competition between 2 species with differing habitat requirements
-   NPDE=3, NU=2
-   front solution can move either direction dependent on parameters
-   epsilon sets steepness of the front (and its propagation speed)

[./FHN](./FHN):

-   Modified FitzHugh--Nagumo cell model
-   NPDE=2, NU=1
-   Initial conditions are smooth, and consistent with BC
-   Propagating pulse solution

[./FHN2](./FHN2):

-   Original Fitzhugh--Nagumo monodomain problem
-   NPDE=2, NU=1
-   Initial conditions such that a travelling pulse is generated

[./rabid](./rabid)

-   Model of rabid fox population (Murray et al. Proc R Soc Lond B Biol Sci. 1986 Nov 22;229(1255):111-50)
-   NPDE=3, NU=1
-   Small perturbation of rabid population at left end produces wave front with an oscillatory tail
-   Due to unstable dynamics of the R=0 solution, eBACOLI solutions can blow up if spatial domain is too large

[./TenTusscher](./TenTusscher):

-   ten Tusscher monodomain model directory
-   NPDE=19, NU=1
-   split to its own directory because of the many files needed for Jacobian computation
-   Initial conditions such that a travelling pulse is generated
