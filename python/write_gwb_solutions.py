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
import cmath

MAXNCHAN = 16384

class GainTable:
    """
    _??.gains file
    """
    def __init__ (self, cfile, default=complex('1.00+0.00j')):
        """
        gain_file 
        """
        self.ichan     = 0
        self.cgains    = dict()

        self.gain_file = cfile

        with open(cfile, 'r') as f:
            lines = [a.strip() for a in f.readlines()]
        ## first line has antennas
        self.ants  = [a.strip()[:3] for a in lines[0].split()]
        self.nants = len(self.ants)
        ## initilaize cgains
        for ant in self.ants: 
            self.cgains[ant] = []
        ## read rest of lines
        for l in lines[1:]:
            _toks = l.split()
            # print ( "line ", _toks, sep='\t' )
            if len(_toks) != self.nants:
                raise RuntimeError("cgains_file is off")
            for ant,cg in zip (self.ants, _toks):
                __cg        = complex ( cg )
                self.cgains[ant].append ( __cg )
            self.ichan += 1
        ##
        for i in range(self.ichan, MAXNCHAN):
            for ant,cg in zip (self.ants, _toks):
                __cg   = default
                self.cgains[ant].append ( __cg )

    def __getitem__(self, key):
        """
        key is (ant, chan)
        """
        ant, ichan = key
        return self.cgains[ant][ichan]

def get_args():
    import argparse
    agp = argparse.ArgumentParser("write_gwb_solutions", description="Writing solutions in GWB format")
    add = agp.add_argument
    add ('tag', help="Stem/tag that was passed to polphase")
    add ('-a','--amp', help='amplitude scaling, computed using :find_gainamp_scaling:', required=True, type=float, dest='ampscaling')
    add ('-O', '--outdir', help='Create an output directory here', default='./', dest='odir')
    add ('-x', '--xphase', help='Put xphase into GWB tables', action='store_true', dest='xphase')
    add ('-l', '--leakage', help='Put determinant of leakages into GWB table', action='store_true', dest='leakage')
    return agp.parse_args()

def process_amp ( c, amp ):
    """ complex to amp in string """
    mag = math.sqrt ( c.real*c.real + c.imag*c.imag )
    if mag <= 1E-3 : mag = 1.0
    return "%.5f" % (amp*mag)

def process_phs ( c ):
    """ complex to phs in string """
    phs = math.degrees ( math.atan2 ( c.imag, c.real ) )
    return "%.3f" % phs

def writer( wfile, wdict ):
    """
    writes to :wfile: given in :wdict:
    """
    ants  = list(wdict.keys())
    ## write amp and phs_file
    with open(wfile, 'w') as f:
        ## write header
        line = "# "
        for ant in ants:
            line += ant
            line += " "
        line += "\n"
        f.write ( line )
        ## write solutions
        for _i in range(MAXNCHAN):
            line = ""
            for ant in ants:
                line += wdict[ant][_i]
                line += " "
            line += "\n"
            f.write ( line )

if __name__ == "__main__":
    args = get_args ()

    grr  = GainTable ( args.tag + "_rr.gains" )
    gll  = GainTable ( args.tag + "_ll.gains" )

    gxp  = None
    if args.xphase:
        gxp  = GainTable ( args.tag + "_xp.gains" )

    drl  = None
    dlr  = None
    if args.leakage:
        drl  = GainTable ( args.tag + "_drl.gains", default=complex('0.00+0.00j') )
        dlr  = GainTable ( args.tag + "_dlr.gains", default=complex('0.00+0.00j') )


    # ODIR = os.path.join ( args.odir, "testtest" )
    ODIR = os.path.join ( args.odir, args.tag )

    if not os.path.exists(ODIR):
        os.mkdir (ODIR)

    """
    modifier = lambda c : c.conjugate()
    Just conjugation works. 

    but in the updated (20261112) logic,
    we do inversion

    so conjugation happens automatically
    """
    ##############################################
    ###    BUILD AMP AND PHASE
    ##############################################
    def action( gt, amp_file, phs_file, add_xphase=False ):
        amp     = dict()
        phs     = dict()
        for a in gt.ants:
            amp[a] = list()
            phs[a] = list()
        ##
        for ichan in range(MAXNCHAN):
            for ant in gt.ants:

                if cmath.isclose ( gt[ant,ichan], 0., abs_tol=1E-7 ):
                    gain = 1.0
                else:
                    gain = 1.0 / gt[ant,ichan]

                if args.leakage:
                    leakrl = drl [ ant, ichan ] 
                    leaklr = dlr [ ant, ichan ]
                    det    = 1.0 - leakrl*leaklr
                    gain   = gain / det

                if add_xphase and args.xphase:
                    _xp    = xp[ant, ichan]
                    gain   = gain / _xp

                amp[ant].append ( process_amp ( gain, args.ampscaling ) )
                phs[ant].append ( process_phs ( gain ) )

        writer ( amp_file, amp )
        writer ( phs_file, phs )
    ##############################################

    sol_amp_130  = os.path.join ( ODIR, "amp.130.dat" )
    sol_amp_175  = os.path.join ( ODIR, "amp.175.dat" )
    sol_phs_130  = os.path.join ( ODIR, "phas.130.dat" )
    sol_phs_175  = os.path.join ( ODIR, "phas.175.dat" )

    action ( grr, sol_amp_130, sol_phs_130, True )
    action ( gll, sol_amp_175, sol_phs_175 )

