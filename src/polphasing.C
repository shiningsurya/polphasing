#include "polphasing.hpp"


int PolPhasing::solve_chan ( int ichan ) {

	/* make data */
	vc_type  chan_data ( nbaselines );

	/* make indices mapping */
	vi_type  xb1 ( nbaselines );
	vi_type  xb2 ( nbaselines );

	auto modeli  = Imodel [ ichan ]; 
	auto modelq  = Qmodel [ ichan ]; 
	auto modelu  = Umodel [ ichan ]; 

	data_t  solve_pkg;

	solve_pkg.data  = chan_data;
	solve_pkg.model = model;

	/* these do not change with ichan */
	solve_pkg.index_b1  = xb1;
	solve_pkg.index_b2  = xb2;
	

	for (int ib = 0; ib < nbaselines; ib++) {

		/* data */
		chan_data [ ib ] = data [ ichan + ib*nchannels ];

		/* indices */
		xb1 [ ib ]       = ;
		xb2 [ ib ]       = ;

		/* model */
		model [ ib ]     = ...
	}

	/* initialize solver here */

	polphasing::full_polar_fcn ( 
		reinterpret_cast<void*>(&solve_pkg),
		ndata, 
		npar,

	);




}
