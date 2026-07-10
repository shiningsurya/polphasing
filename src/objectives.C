#include "objectives.hpp"

#include <cstring>


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


		/*
		 * m  is 2 * polar_baselines
		 * hm is polar_baselines
		 *
		 * data, model, index_b1, index_b2 all are complex length=hm
		 *
		 * but LMSolver sees real residuals
		 * we do this kind of for loop
		 */
		int hm     = m / 2;

		if ( iflag == 0 ) {
#if 0
			std::cout << " jac=" <<  fjac[0] << " " << fjac[1] << std::endl;
			// i do not get why jacobian structure is lost towards the end?
			{
				std::ofstream of("everyjac.bin", std::ios::binary | std::ios::app);
				of.write (reinterpret_cast<const char*>(fjac), m*n*sizeof(float));
			}
#endif
		} /* printing */
		else if ( iflag == 1 ) {
			/* fvec computation */

			//#pragma omp parallel for num_threads(4)
			//openmp here does not do much
			for (int im = 0; im < hm; im++) {
				/* get antenna indices */
				int ib1 = index_b1 [ im ];
				int ib2 = index_b2 [ im ];

				/* output */
				int om  = 2 * im;

				/* unrolling the math */
				/* see res_math_forC */
				/* hoping that it speeds up a bit */
				float pr  = x [ 2*ib1 + 0 ];
				float pi  = x [ 2*ib1 + 1 ];
				float qr  = x [ 2*ib2 + 0 ];
				float qi  = x [ 2*ib2 + 1 ];

				float dr  = data[im].real();
				float di  = data[im].imag();

				float mr  = model[im].real();
				float mi  = model[im].imag();

				/* real and imag */
				fvec[om]  = dr + (mi*pi*qr) - (mi*pr*qi) - (mr*pi*qi) - (mr*pr*qr);
				fvec[om+1]= di - (mi*pi*qi) - (mi*pr*qr) - (mr*pi*qr) + (mr*pr*qi);

				/*
				 * real(res) 	d^r + m^i*p^i*q^r - m^i*p^r*q^i - m^r*p^i*q^i - m^r*p^r*q^r
 					 imag(res) 	d^i - m^i*p^i*q^i - m^i*p^r*q^r - m^r*p^i*q^r + m^r*p^r*q^i
				 */
#if 0
				/* get gains */
				complex_type p1  ( x[2*ib1 + 0], x[2*ib1 + 1] );
				/* get gains - directly complement */
				complex_type q2c ( x[2*ib2 + 0],-x[2*ib2 + 1] );

				/* data, model */
				complex_type idata    = data  [ im ];
				complex_type imodel   = model [ im ];
				/* model prediction */
				complex_type omodel   = p1 * imodel * q2c;
				/* residual */
				complex_type res      = idata - omodel;
				/* error = DATA - MODEL */
				/* saving as real and imaginary part */
				fvec [om]    = res.real();
				fvec [om+1]  = res.imag();
#endif

			} /* for every polar_baseline */

		} /* fvec computation */
		else if ( iflag == 2 ) {
			/* fjac computation */

			/*
			 * i reiterate, 
			 * m  is 2*polar_baselines
			 * hm is polar_baselines
			 *
			 * data,model,indices are size=hm
			 */

			/* do i need to zero out the jacobian everytime? */
			/* yes */
			std::memset ( fjac, 0, m*n*sizeof(real_type) );

			//#pragma omp parallel for num_threads(4)
			//openmp here does not do much
			for (int im = 0; im < hm; im++) {

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
				 * getting this from :jac_math_forC.py:
				 * t is first, s is second
				 *
						 real (res) / tr 	m^i*s^i + m^r*s^r
						 imag (res) / tr 	m^i*s^r - m^r*s^i
						 real (res) / ti 	-m^i*s^r + m^r*s^i
						 imag (res) / ti 	m^i*s^i + m^r*s^r
						 real (res) / sr 	-m^i*t^i + m^r*t^r
						 imag (res) / sr 	m^i*t^r + m^r*t^i
						 real (res) / si 	m^i*t^r + m^r*t^i
						 imag (res) / si 	m^i*t^i - m^r*t^r
				 *
				 * Jacobian shape is (ndata, npar)
				 *
				 * But minpack examples and my test suggests it should be
				 * (npar, ndata)
				 * Even scipy.optimize.least_squares also suggest that 
				 * shape should be (n,m)
				 *
				 * `polarbaselines` axis is the fastest
				 * so ldfjac = m,
				 *
				 * assuming so, this should be the ordering
				 * that matches with above math
				 * 	(tr, re(res)) = real(res) / tr
				 * 	(tr, im(res)) = imag(res) / tr
				 *
				 *
				 */

				int idx = 0;
				/* re,im(res) / tr */
				idx  = ldfjac*i1r + 2*im;
				fjac [ idx     ] = - 1.0 * ( imodelr*s2r + imodeli*s2i );
				fjac [ idx + 1 ] = - 1.0 * (-imodelr*s2i + imodeli*s2r );
				//fjac [ ldfjac*i1r + 2*im     ] = - 1.0 * ( imodelr*s2r + imodeli*s2i );
				//fjac [ ldfjac*i1r + 2*im + 1 ] = - 1.0 * (-imodelr*s2i + imodeli*s2r );

				/* re,im(res) / ti */
				idx  = ldfjac*i1i + 2*im;
				fjac [ idx     ] = - 1.0 * ( imodelr*s2i - imodeli*s2r );
				fjac [ idx + 1 ] = - 1.0 * ( imodelr*s2r + imodeli*s2i ); 
				//fjac [ ldfjac*i1i + 2*im     ] = - 1.0 * ( imodelr*s2i - imodeli*s2r );
				//fjac [ ldfjac*i1i + 2*im + 1 ] = - 1.0 * ( imodelr*s2r + imodeli*s2i ); 

				/* re,im(res) / sr */
				idx  = ldfjac*i2r + 2*im;
				fjac [ idx     ] = - 1.0 * ( imodelr*t1r - imodeli*t1i );
				fjac [ idx + 1 ] = - 1.0 * ( imodeli*t1r + imodelr*t1i );
				//fjac [ ldfjac*i2r + 2*im     ] = - 1.0 * ( imodelr*t1r - imodeli*t1i );
				//fjac [ ldfjac*i2r + 2*im + 1 ] = - 1.0 * ( imodeli*t1r + imodelr*t1i );

				/* re,im(res) / si */
				idx  = ldfjac*i2i + 2*im;
				fjac [ idx     ] = - 1.0 * ( imodelr*t1i + imodeli*t1r );
				fjac [ idx + 1 ] = - 1.0 * (-imodelr*t1r + imodeli*t1i );

				/* these indices are correct */
			//std::cout << ldfjac*i1r + 2*im << "," << ldfjac*i1r + 2*im + 1 << ",";
			//std::cout << ldfjac*i1i + 2*im << "," << ldfjac*i1i + 2*im + 1 << ",";
			//std::cout << ldfjac*i2r + 2*im << "," << ldfjac*i2r + 2*im + 1 << ",";
			//std::cout << ldfjac*i2i + 2*im << "," << ldfjac*i2i + 2*im + 1 << std::endl;

			} /* for every polar_baseline */

		} /* fjac computation */

		return iflag;

	} /* full_polar_fcn */

}; /* polphasing */
