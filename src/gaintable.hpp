#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <complex>
#include <map>
#include <fstream>
#include <iomanip>

/* fmt library */
//#define FMT_HEADER_ONLY
//#include "fmt/base.h"
//#include "fmt/core.h"
//#include "fmt/format.h"
//#include "fmt/ostream.h"

namespace gaintable {
	using real_type   = float;
	using antname_t   = std::array<char,4>;
	using gain_type   = std::complex<real_type>;
	using gains_t     = std::vector<gain_type>;
	using gaintable_t = std::map< antname_t, gains_t>;

	/* ordering in which to write */
	//static const std::vector<antname_t> sol_ant_order (
			//antname_t {"C", "0", "0", ""},
			//antname_t {"C", "0", "1", ""},
			//antname_t {"C", "0", "2", ""},
			//antname_t {"C", "0", "3", ""},
			//antname_t {"C", "0", "4", ""},
			//antname_t {"C", "0", "5", ""},
			//antname_t {"C", "0", "6", ""},
			//antname_t {"C", "0", "8", ""},
			//antname_t {"C", "0", "9", ""},
			//antname_t {"C", "1", "0", ""},
			//antname_t {"C", "1", "1", ""},
			//antname_t {"C", "1", "2", ""},
			//antname_t {"C", "1", "3", ""},
			//antname_t {"C", "1", "4", ""},
			//antname_t {"E", "0", "2", ""},
			//antname_t {"E", "0", "3", ""},
			//antname_t {"E", "0", "4", ""},
			//antname_t {"E", "0", "5", ""},
			//antname_t {"E", "0", "6", ""},
			//antname_t {"S", "0", "1", ""},
			//antname_t {"S", "0", "2", ""},
			//antname_t {"S", "0", "3", ""},
			//antname_t {"S", "0", "4", ""},
			//antname_t {"S", "0", "6", ""},
			//antname_t {"W", "0", "1", ""},
			//antname_t {"W", "0", "2", ""},
			//antname_t {"W", "0", "3", ""},
			//antname_t {"W", "0", "4", ""},
			//antname_t {"W", "0", "5", ""},
			//antname_t {"W", "0", "6", ""},
			//antname_t {"C", "0", "7", ""},
			//antname_t {"S", "0", "5", ""}
	//);
	static const std::vector<antname_t> sol_ant_order {
			antname_t {"C00"},
			antname_t {"C01"},
			antname_t {"C02"},
			antname_t {"C03"},
			antname_t {"C04"},
			antname_t {"C05"},
			antname_t {"C06"},
			antname_t {"C08"},
			antname_t {"C09"},
			antname_t {"C10"},
			antname_t {"C11"},
			antname_t {"C12"},
			antname_t {"C13"},
			antname_t {"C14"},
			antname_t {"E02"},
			antname_t {"E03"},
			antname_t {"E04"},
			antname_t {"E05"},
			antname_t {"E06"},
			antname_t {"S01"},
			antname_t {"S02"},
			antname_t {"S03"},
			antname_t {"S04"},
			antname_t {"S06"},
			antname_t {"W01"},
			antname_t {"W02"},
			antname_t {"W03"},
			antname_t {"W04"},
			antname_t {"W05"},
			antname_t {"W06"},
			antname_t {"C07"},
			antname_t {"S05"} 
	};

	/* initialize gaintables */
	gaintable_t prepare_gaintables ( int nchan );

	gaintable_t read_complex_solutions ( const std::string& infile );

	int write_complex_solutions ( const gaintable_t& gt, const std::string& outfile );

	int write_phase_solutions ( const gaintable_t& gt );

	int write_amp_solutions ( const gaintable_t& gt );

	struct GainNormalizer {
		/**/
		using real_type   = double;
		using gain_type   = std::complex<real_type>;
		using jones_t     = std::array<gain_type,4>;

		const jones_t        jones_ref;
		const gain_type      det;
		const real_type      det_abs;
		const gain_type      det_phs;
		const real_type      unitary_scale;
		const jones_t        unitary_t;

		/* unitary_t is not full unitary matrix, unitary = unitary_t * unitary_scale
		 *
		 * It has already been transpose, conjugated.
		*/

		GainNormalizer ( const gain_type _rr, const gain_type _rl, const gain_type _lr, const gain_type _ll ) :
			jones_ref {_rr, _rl, _lr, _ll},
			det ( _rr*_ll - _rl*_lr ),
			det_abs ( std::abs(det) ),
			det_phs ( std::conj(det) / det_abs ),
			unitary_scale ( std::pow( 2.0*det_abs + std::norm(_rr) + std::norm(_rl) + std::norm(_lr) + std::norm(_ll), -0.5 ) ),
			unitary_t {
				std::conj(_rr) + det_phs*_ll,
				std::conj(_lr) - det_phs*_rl,
				std::conj(_rl) - det_phs*_lr,
				std::conj(_ll) + det_phs*_rr
			}
		{}

		/* apply  
		 *
		 * Multiply unitary_t from right side
		 * do not forget the unitary scale
		*/
		jones_t operator() ( const gain_type _rr, const gain_type _rl, const gain_type _lr, const gain_type _ll ) const {
			jones_t ret;

			// 00
			ret[0] = _rr*unitary_t[0] + _rl*unitary_t[2] ;
			// 01
			ret[1] = _rr*unitary_t[1] + _rl*unitary_t[3] ;
			// 10
			ret[2] = _lr*unitary_t[0] + _ll*unitary_t[2] ;
			// 11
			ret[3] = _lr*unitary_t[1] + _ll*unitary_t[3] ;

			// scaling
			ret[0] *= unitary_scale;
			ret[1] *= unitary_scale;
			ret[2] *= unitary_scale;
			ret[3] *= unitary_scale;

			return ret;
		}

	}; // struct gainnormalizer
	using normalizer_t = struct GainNormalizer;

}; /* gaintable */

/*

// custom formatter 
template<>
struct fmt::formatter<gaintable::gain_type> : fmt::formatter<gaintable::real_type> {

	template<typename FormatContext>
		auto format(const gaintable::gain_type& g, FormatContext& ctx) const -> decltype(ctx.out()) {
			auto out = ctx.out ();

			out      = fmt::formatter<gaintable::real_type>::format (g.real(), ctx);

			out      = fmt::formatter<gaintable::real_type>::format (g.imag(), ctx);

			return fmt::format_to ( out, "j" );
		}
};

template<>
struct fmt::formatter<gaintable::antname_t> : fmt::formatter<string_view> {

	template<typename FormatContext>
		auto format(const gaintable::antname_t& g, FormatContext& ctx) const -> decltype(ctx.out()) {
			return fmt::format_to (ctx.out(), "{}", string_view(g.data(), g.size()));
		}
};

*/

