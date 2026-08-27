"""

polphase saves complex gains in "{tag}_r.gains" and "{tag}_l.gains"

this script writes amplitude and phase solutions in 
"amp_{130,175}_.dat"

why shouldn't polphase write it itself?
idk

this script will probably run in a place with no packages
this should be barebones

by default, the solutions are put in a separate directory because the filename convention 
of the solutions read by GWB is the same throughout
"""

import os
import sys
import math

def get_args():
    import argparse
    agp = argparse.ArgumentParser("write_gwb_solutions", description="Writing solutions in GWB format")
    add = agp.add_argument
    add ('tag', help="Stem/tag that was passed to polphase")
    add ('-O', '--outdir', help='Create an output directory here', default='./', dest='odir')
    return agp.parse_args()


def process_amp ( c ):
    """ complex to amp in string """
    mag = math.sqrt ( c.real*c.real + c.imag*c.imag )
    return "%.1f" % mag

def process_phs ( c ):
    """ complex to phs in string """
    phs = math.degrees ( math.atan2 ( c.imag, c.real ) )
    return "%.1f" % phs

def action ( cgains_file, amp_file, phs_file, modifier ):
    """
    read from cgains_file
    write to amp and phs file
    """
    cgains    = dict()
    ## amp, phs : antenna : [channels]
    amp       = dict()
    phs       = dict()
    ichan     = 0
    with open(cgains_file, 'r') as f:
        lines = [a.strip() for a in f.readlines()]
        ##
        ## first line has antennas
        ants  = [a.strip()[:3] for a in lines[0].split()]
        nants = len(ants)
        #print (" read ants", ants, sep='\t')
        ## initilaize cgains
        for ant in ants: 
            cgains[ant] = []
            amp[ant]    = []
            phs[ant]    = []
        ## read rest of lines
        for l in lines[1:]:
            _toks = l.split()
            # print ( "line ", _toks, sep='\t' )
            if len(_toks) != nants:
                raise RuntimeError("cgains_file is off")
            for ant,cg in zip (ants, _toks):
                __cg        = modifier ( complex ( cg ) )
                """
                in GWB correlator, we need to apply the inverse of the gain we compute from polphase.
                Which in case of diagonal complex is simply the reciprocal.
                That is, exp(-j phi_pp)

                in GWB correlator, the conjugate of the phase table is applied
                exp(-j phi_tt)

                That means phi_pp = phi_tt

                But when i compare rantsol and polphase, 
                there is a phase difference.

                so we do trial and error.
                """
                cgains[ant].append ( __cg )
                amp[ant].append ( process_amp ( __cg ) )
                phs[ant].append ( process_phs ( __cg ) )
                # print ( __cg, process_amp(__cg), process_phs(__cg), sep='\t' )
            ichan += 1
        # print (" read nchans = ", ichan)
    #######################
    ## the amp,phs need to go until 16K
    ## the max channel limit
    for _ in range(ichan,16384):
        __cg   = complex('1.00+0.00j')
        for ant in ants:
            amp[ant].append ( process_amp ( __cg ) )
            phs[ant].append ( process_phs ( __cg ) )
    #######################
    ## write amp and phs_file
    with open(amp_file, 'w') as f:
        ## write header
        line = "# "
        for ant in ants:
            line += ant
            line += " "
        line += "\n"
        f.write ( line )
        ## write solutions
        for _i in range(16384):
            line = ""
            for ant in ants:
                line += amp[ant][_i]
                line += " "
            line += "\n"
            f.write ( line )

    with open(phs_file, 'w') as f:
        ## write header
        line = "# "
        for ant in ants:
            line += ant
            line += " "
        line += "\n"
        f.write ( line )
        ## write solutions
        for _i in range(16384):
            line = ""
            for ant in ants:
                line += phs[ant][_i]
                line += " "
            line += "\n"
            f.write ( line )

    return cgains, amp, phs


if __name__ == "__main__":
    args = get_args ()

    cgains_130   = args.tag + "_r.gains"
    cgains_175   = args.tag + "_l.gains"

    # ODIR = os.path.join ( args.odir, "testtest" )
    ODIR = os.path.join ( args.odir, args.tag )

    if not os.path.exists(ODIR):
        os.mkdir (ODIR)

    ##
    ## modifier
    modifier = lambda c : c.conjugate()
    #modifier = lambda c : 1.0 / c

    sol_amp_130  = os.path.join ( ODIR, "amp.130.dat" )
    sol_amp_175  = os.path.join ( ODIR, "amp.175.dat" )
    sol_phs_130  = os.path.join ( ODIR, "phas.130.dat" )
    sol_phs_175  = os.path.join ( ODIR, "phas.175.dat" )

    rgains, ramps, rphs = action ( cgains_130, sol_amp_130, sol_phs_130, modifier )
    action ( cgains_175, sol_amp_175, sol_phs_175, modifier )






