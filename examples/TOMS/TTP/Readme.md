Solution to the monodomain equation with the ten Tusscher--Panfilov model of human epicardial cells. Source code for evaluating the TTP model is found in the [./TTP\_src](./TTP_src) directory.

Running `make` will compile, run, generate the frames, and produce a movie of output (if you have `avconv`). It will then prompt for which directory to save all of this data in.

To avoid the prompt, you can just run `echo <dirName> | make`.
