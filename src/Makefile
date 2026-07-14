
CMINPACK=lmder.c dpmpar.c enorm.c qrfac.c lmpar.c qrsolv.c

LMSOLVER=lmsolver.C objectives.C

LTA=lta_file.C myltasub.c mycorrsub.c sort.c

MISC=gaintable.C models.C logsolve.C


# test_lmsolver

test_lmsolver:
	g++ -g -std=c++11 -D__cminpack_float__ test_lmsolver.C $(LMSOLVER) $(CMINPACK)  -o $@

test_lta:
	g++ -g -std=c++11 -D__cminpack_float__ test_lta.C $(LTA) -o $@

clean:
	rm -f test_lmsolver 

test_solving:
	g++ -g -std=c++11 -D__cminpack_float__ test_solving.C $(LMSOLVER) $(CMINPACK)  -o $@

test_packs:
	g++ -g -o $@ -std=c++11 -D__cminpack_float__ -DUSE_LAPACK -DUSE_BLAS \
		-L/usr/lib \
		test_solving.C $(LMSOLVER) $(CMINPACK) -llapack -lblas 

#-llapack -lblas -lgfortran

## will be compiling everything at once
## because ~i am nuts~
## -DTIMING only when needed
polphase: polphase.C $(LMSOLVER) $(CMINPACK)
	g++ -O2 -o $@ -g -fopenmp -std=c++11 -DTIMING -DUSE_LAPACK -DUSE_BLAS -D__cminpack_float__ polphase.C $(LMSOLVER) $(CMINPACK) $(LTA) $(MISC) -llapack -lblas

polphaseomp: polphase_omp.C $(LMSOLVER) $(CMINPACK)
	g++ -O2 -o $@ -g -fopenmp -std=c++11 -DTIMING -DUSE_LAPACK -DUSE_BLAS -D__cminpack_float__ polphase_omp.C $(LMSOLVER) $(CMINPACK) $(LTA) $(MISC) -llapack -lblas

write_baselines: write_baselines.C
	g++ -O0 -o $@ -g -std=c++11 write_baselines.C $(LTA)

