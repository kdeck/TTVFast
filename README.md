TTVFast
=======

TTVFast efficiently calculates transit times for n-planet systems and the corresponding radial velocities.


Available Versions
=======

The C version of the code is in the directory c_version, the Fortran version is in fortran_version. Both versions have specific README files.

There is a Julia interface in the directory jl_version.  First, you must create libttvfast.so (e.g., "cd jl_version; source compile_libttvfast.cmd" on Linux).  Once you've done this, start julia and type
'include("demo_julia.jl")' to see an example of calling ttvfast and accessing the outputs from Julia.  There's another example of using the same input file foramts as TTVFast.

Citations
=======
If you use this code, please cite Deck, Agol, Holman, & Nesvorny (2014), ApJ, 787, 132, arXiv:1403.1895. 



Please check back for updates to ensure that you are using the latest version.

-Katherine Deck, Eric Agol, Matt Holman, & David Nesvorny