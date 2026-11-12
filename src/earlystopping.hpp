#pragma once

namespace LBFGS {

	struct EarlyStopping {
		using cost_type   = float;

		static constexpr int        nsteps   = 64;
		static constexpr cost_type  _epsilon = 1e-5;

		int                  idx;
		int                  past;
		std::vector<float>   costs;

		EarlyStopping ( int _past = 1) : idx(0), past (_past), costs (nsteps, 0.0) {}

		void push_back ( const cost_type& cost ) {
			costs[idx] = cost;
			idx++;
		}

		explicit operator bool() const {
			int pdx = std::max ( idx - past, idx - 1 );

			//  we do not have enough metrics to decide
			if (pdx < 0) return false;

			const cost_type diff ( costs[pdx] - costs[idx] );

			// we expect costs[idx] < costs[pdx]
			// so if diff is negative, quickly return true
			// we are probably overfitting then
			// AND
			// if the difference is smaller than some cutoff
			// return true
			if ( diff < _epsilon ) return true;

			// default
			return false;
		}

	}; // struct early stopping

}; // namespace LBFGS
