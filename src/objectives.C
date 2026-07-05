#include "objectives.hpp"


namespace polphasing {

	using data_t = polphasing_data_t;

	/*
	 * This is full real.
	 *
	 * m is ndata, n is npar
	 *
	 * m is 2 * nbaselines (polar baselines)
	 *
	 * polar baselines already consider the four correlations.
	 *
	 * some baseline `im` would have ant1_index[im] and ant2_index[im].
	 *
	 * antk_index[im] has four gains 4*antk_index[im] + {0,1,2,3}
	 *
	 *
	 */
	int full_polar_fcn (void* vd, int m, int n, const real_type* x, real_type* fvec, real_type* fjac, int ldfjac, int iflag) {
		/* m is ndata, n is npar */

		/* pull data */
		data_t *pkg = reinterpret_cast<data_t*> ( vd );

		/* model 
		 *
		 * RR = model_i
		 * RL = model_q + j*model_u
		 * LR = model_q - j*model_u
		 * LL = model_i
		 *
		 * This will be like an array to save some logic
		 */

		/* complex arrays */
		vc_type& data   = pkg->data;
		vc_type& model  = pkg->model;

		/* mapping */
		vi_type& index_b1 = pkg->index_b1;
		vi_type& index_b2 = pkg->index_b2;

		if ( iflag == 0 ) {

		} /* printing */
		else if ( iflag == 1 ) {
			/* fvec computation */

			for (int im = 0; im < m; im+=2) {
				/* get antenna indices */
				int ib1 = index_b1 [ im ];
				int ib2 = index_b2 [ im ];

				/* get gains */
				complex_type p1  ( x[2*ib1 + 0], x[2*ib1 + 1] );
				/* get gains - directly complement */
				complex_type q2c ( x[2*ib2 + 0],-x[2*ib2 + 1] );

				/* data, model */
				complex_type idata    = data  [ im ];
				complex_type imodel   = model [ im ];
				/* model prediction */
				complex_type observed = p1 * imodel * q2c;
				/* residual */
				complex_type res      = idata - observed;
				/* error = DATA - MODEL */
				/* saving as real and imaginary part */
				fvec [im]    = res.real();
				fvec [im+1]  = res.imag();

			} /* for every polar_baseline */

		} /* fvec computation */
		else if ( iflag == 2 ) {
			/* fjac computation */

			for (int im = 0; im < m; im+=2) {
				/* get antenna indices */
				int ib1   = index_b1 [ im ];
				int ib2   = index_b2 [ im ];

				int i1r    = 2*ib1 + 0;
				int i1i    = 2*ib1 + 1;
				int i2r    = 2*ib2 + 0;
				int i2i    = 2*ib2 + 1;

				/* get gains */
				real_type  t1r    = x[i1r];
				real_type  t1i    = x[i1i];
				real_type  s2r    = x[i2r];
				real_type  s2i    = x[i2i];

				/* model */
				complex_type imodel   = model [ im ];
				real_type  imodelr    = imodel.real();
				real_type  imodeli    = imodel.imag();

				/* error = DATA - MODEL */
				/* saving as real and imaginary part */
				/* there are 4 non zero elements per im */
				/* there are then 8 because we have complex64 */

				/*
				 * REMATH
				 *
				 * t,s = {r,l}
				 *
				 * res_t1s2 = data_t1s2 - t1 * model_ts * s2c
				 *
				 * (ignore data because derivative hides it)
				 * (assume minus sign.)
				 *
				 * (a,b)* (c,d) = ( ac - bd, ad + bc )
				 *
				 * (t1r, t1i) * (modelr, modeli) * (s2r, -s2i)
				 *
				 * (t1r, t1i) * (modelr*s2r + modeli*s2i, -modelr*s2i + modeli*s2r)
				 *
				 * (
				 * 		t1r*modelr*s2r + t1r*modeli*s2i + t1i*modelr*s2i - t1i*modeli*s2r,
				 * 	 -t1r*modelr*s2i + t1r*modeli*s2r + t1i*modelr*s2r + t1i*modeli*s2i
				 * )
				 *
				 * re(res) / t1r = modelr*s2r + modeli*s2i
				 * im(res) / t1r =-modelr*s2r + modeli*s2r
				 *
				 * re(res) / t1i = modelr*s2i - modeli*s2r
				 * im(res) / t1i = modelr*s2r + modeli*s2i
				 *
				 * re(res) / s2r = t1r*modelr - t1i*modeli
				 * im(res) / s2r = t1r*modeli + t1i*modelr
				 *
				 * re(res) / s2i = t1r*modeli + t1i*modelr
				 * im(res) / s2i =-t1r*modelr + t1i*modeli
				 *
				 */

				/* gpr|hpr / re|im */
				fjac [ ldfjac*i1r + im   ] = - 1.0 * ( imodelr*s2r + imodeli*s2i );
				fjac [ ldfjac*i1r + im+1 ] = - 1.0 * (-imodelr*s2r + imodeli*s2r );

				/* gpi|hpi */
				fjac [ ldfjac*i1i + im   ] = - 1.0 * ( imodelr*s2i - imodeli*s2r );
				fjac [ ldfjac*i1i + im+1 ] = - 1.0 * ( imodelr*s2r + imodeli*s2i ); 

				/* gqr|hqr */
				fjac [ ldfjac*i2r + im   ] = - 1.0 * ( imodelr*t1r - imodeli*t1i );
				fjac [ ldfjac*i2r + im+1 ] = - 1.0 * ( imodeli*t1r + imodelr*t1i );

				/* gqi|hqi */
				fjac [ ldfjac*i2i + im   ] = - 1.0 * ( imodeli*t1r + imodelr*t1i );
				fjac [ ldfjac*i2i + im+1 ] = - 1.0 * (-imodelr*t1r + imodeli*t1i );

			} /* for every polar_baseline */

		} /* fjac computation */

		return iflag;

	} /* full_polar_fcn */

}; /* polphasing */
