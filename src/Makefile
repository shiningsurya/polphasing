
CMINPACK=lmder.c dpmpar.c enorm.c qrfac.c lmpar.c qrsolv.c

LMSOLVER=lmsolver.C objectives.C

LTA=lta_file.C myltasub.c mycorrsub.c sort.c

MISC=gaintable.C models.C

# test_lmsolver

test_lmsolver:
	g++ -g -std=c++11 -D__cminpack_float__ test_lmsolver.C $(LMSOLVER) $(CMINPACK)  -o $@

test_lta:
	g++ -g -std=c++11 -D__cminpack_float__ test_lta.C $(LTA) -o $@

clean:
	rm -f test_lmsolver 

test_solving:
	g++ -g -std=c++11 -D__cminpack_float__ test_solving.C $(LMSOLVER) $(CMINPACK)  -o $@

polphase: polphase.C $(LMSOLVER) $(CMINPACK)
	g++ -g -std=c++11 -D__cminpack_float__ polphase.C $(LMSOLVER) $(CMINPACK) $(LTA) $(MISC) -o $@
