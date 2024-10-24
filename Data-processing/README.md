The 2dVar approach is a sophisticated method for processing High Frequency (HF) radar data to map sea surface currents, which was developed by Max Yaremchuk and Sentchev [1] and was updated by Tran Manh Cuong during his PhD [2] . The algorithm is written in FORTRAN. The main function of the algorithm is ass100.f. 

To run the program, the user has to prepare the grid fil in ASCII: mask.dat and mask11.dat, which cover the region of study. 
To compile the program, a FORTRAN compiler is required. The example of compiling the software is shown below: 

Ex: 
# With g95 compiler
g95 -ffixed-line-length-132 -fendian=big -r8 -freal=zero ass100.f m1qn3.f -o 2dVar.exe
# With ifort compiler
ifort -O3 -extend-source 132 -real-size 64 -zero -mcmodel=medium ass100.f m1qn3.f -o 2dVar
# With gfortran compiler
gfortran -ffixed-line-length-132 -finit-real=zero -fdefault-real-8 ass100.f m1qn3.f -o 2dVar
# gfortran 10
gfortran -ffixed-line-length-132 -finit-real=zero -fallow-argument-mismatch -fdefault-real-8 ass100.f m1qn3.f -o 2dVar

# Ref: 
[1] Yaremchuk, M., Sentchev, A., 2009. Mapping radar-derived sea surface currents with a variational method. Continental Shelf Research 29,
1711–1722. doi:https://doi.org/10.1016/j.csr.2009.05.016
[2] Tran, M. C. (2022). Caractérisation de la dynamique côtière et de la dispersion turbulente dans le golfe du Tonkin à partir de la courantographie radar HF et de la modélisation : effet de la dynamique à fine-échelle sur la structuration spatiale du phytoplancton http://www.theses.fr/2022DUNK0626/document