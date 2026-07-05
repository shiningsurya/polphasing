
CMINPACK=lmder.c dpmpar.c enorm.c qrfac.c lmpar.c qrsolv.c

LMSOLVER=lmsolver.C objectives.C

LTA=lta_file.C myltasub.c mycorrsub.c sort.c

# test_lmsolver

test_lmsolver:
	g++ -g -std=c++11 -D__cminpack_float__ test_lmsolver.C $(LMSOLVER) $(CMINPACK)  -o $@

test_lta:
	g++ -g -std=c++11 -D__cminpack_float__ test_lta.C $(LTA) -o $@

clean:
	rm -f test_lmsolver 
