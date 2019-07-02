Brief description of extended examples found in this directory.

[./toy3.f](./toy3.f):

-   Toy model with NPDE=4, NU=2.
-   Obtained by duplicating all of toy2.f, but replacing one of the boundary conditions to use spatial derivative of v
-   Has the same solution as toy2.f

[./FHN\_monodomain.f](./FHN_monodomain.f):

-   modified FitzHugh--Nagumo cell model from Elham's thesis (results in Table 4.3)
-   NPDE=2, NU=1

[./FHN2\_monodomain.f](./FHN2_monodomain.f):

-   Original Fitzhugh--Nagumo monodomain problem
-   NPDE=2, NU=1
-   Initial conditions such that a travelling pulse is generated

[./TenTusscher](./TenTusscher):

-   ten Tusscher monodomain model directory
-   NPDE=19, NU=1
-   split to its own directory because of the many files needed for Jacobian computation
-   Initial conditions such that a travelling pulse is generated

[./CDM.f](./CDM.f):

-   Ecological competition between 2 species with differing habitat requirements
-   NPDE=3, NU=2
-   front solution can move either direction dependent on parameters
-   epsilon sets steepness of the front (and its propagation speed)

[./rabid.f](./rabid.f)

-   Model of rabid fox population (Murray et al. Proc R Soc Lond B Biol Sci. 1986 Nov 22;229(1255):111-50)
-   NPDE=3, NU=1
-   Small perturbation of rabid population at left end produces wave front with an oscillatory tail
-   Due to unstable dynamics of the R=0 solution, eBACOLI solutions can blow up if spatial domain is too large
