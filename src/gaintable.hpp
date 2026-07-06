#pragma once

#include <array>
#include <vector>
#include <string>
#include <complex>
#include <map>
#include <fstream>

/* fmt library */
#define FMT_HEADER_ONLY
#include "fmt/base.h"
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/ostream.h"

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

	int write_complex_solutions ( const gaintable_t& gt, const std::string& outfile );

	int write_phase_solutions ( const gaintable_t& gt );

	int write_amp_solutions ( const gaintable_t& gt );

}; /* gaintable */

/* custom formatter */
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

