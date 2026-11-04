/*
 * LTA interface impl
 */

#include "lta_file.hpp"
#include "sources.hpp"
#include <stdexcept>


LTA::LTA (const std::string &file) : filepath (file), fp (nullptr) {

	/* open file */
	fp = fopen ( filepath.c_str(), "r" );
	if ( fp == NULL ) throw std::runtime_error ("Unable to open given LTA file.");

	/* read header */
	if ( get_lta_hdr ( fp, &linfo ) ) throw std::runtime_error ("Unable to read hdr file");

	/* make scan table */
	if ( make_scantab ( fp, &linfo ) ) throw std::runtime_error ("Unable to make scan table");

	recl        = linfo.recl;

	/* load header data */
	auto vinfo  = linfo.lhdr.vinfo;
	nbaselines  = vinfo.baselines;
	nchannels   = vinfo.channels;
	nantennas   = vinfo.antennas;
	nsamplers   = vinfo.samplers;

	/* load antennas */
	/* not nantennas here because we want to read from ScanHdr.antmask */
	/* we are adventurous. */
	for (int ia = 0; ia < MAX_ANTS; ia++) {
		const AntennaParType& iant = vinfo.antenna[ia];

		antname_t iname;
		std::copy ( std::begin(iant.name), std::end(iant.name), iname.begin() );

		/*just the ant names*/
		all_antnames.push_back ( iname );
	}

	/* load baselines */
	for (int ib = 0; ib < nbaselines; ib++) {
		const BaselineType& bib  = vinfo.base[ib];

		baselines.emplace_back ( bib );
	}

	/* reserve space */
	//drow.reserve (recl);
}

LTA::~LTA () {

	/* need to clean memory */
	free ( linfo.stab );

	/* close file */
	if (fp != NULL) { fclose ( fp ); }

}

void LTA::time_average (int iscan, vf_type& inout) {
	/* sanity check */
	if (iscan >= linfo.scans) throw std::runtime_error ("Chosen scan out of range.");

	int  data_off  = linfo.lhdr.data_off;
	int  data_size = linfo.lhdr.data_size;

	/* inout should be nbaselines*nchannels complex64 */
	/* inout is actually float32, so we double  */

	/* another sanity check*/
	/* nd = (nbaselines * nchannels) */
	int    nd    = data_size / sizeof(float);
	if ( inout.size() != nd ) throw std::runtime_error ("Vector size should be ltahdr.data_size.");
	inout.reserve ( 2 * nbaselines * nchannels );

	/* seek in file */
	rewind (fp);
	ltaseek (fp, linfo.stab[iscan].start_rec + linfo.srecs, recl);

	/* nrecs records for this scan */
	int nrecs  = linfo.stab[iscan].recs;

	/* space to read */
	vb_type dbuf ( recl );

	
	for (int irec = 0; irec < nrecs; irec++) {

		/* read record from file */
		if ( fread ( dbuf.data(), recl, 1, fp ) != 1 ) 
			throw std::runtime_error ("Incomplete read while time averaging.");

		/* data offset, cast as float */
		float *ddata = reinterpret_cast<float*>(dbuf.data() + data_off);

		/* add to inout */
		//std::transform ( inout.begin(), inout.end(), ddata );
		std::transform ( inout.data(), inout.data() + nd, ddata, inout.data(), std::plus<float>() );

	} /* rec loop */

	/* rescale or not */
	std::transform ( inout.begin(), inout.end(), inout.begin(), [nrecs] (float &i) { return i / nrecs; } );
}

std::ostream& operator<< (std::ostream& os, const LTA & l) {
	os << "i "  << "Scan  " << "Nrecs" << std::endl;

	for ( int iscan = 0; iscan < l.linfo.scans; iscan++ ) {

		unsigned int recs = l.linfo.stab[iscan].recs;
		
		ScanHdr& shdr     = l.linfo.stab[iscan].shdr;

		os << iscan <<"  " << shdr.object << " " << recs << std::endl;
		
	}

	return os;
}

LTA::scan_t LTA::get_scan ( int iscan ) {
	/*
	 * Get scan properties of :iscan:
	 *
	 * need source name, ra/dec, mjd
	 *
	 * get middle mjd
	*/
	/* sanity check */
	if (iscan >= linfo.scans) throw std::runtime_error ("Chosen scan out of range.");

	const ScanHdrTab& stab = linfo.stab[iscan];
	const ScanHdr& shdr    = stab.shdr;

	const std::string source ( shdr.object );
	const float       integ  ( shdr.integ );
	const int         nrecs  ( stab.recs );

	/* we want time at middle of the scan */
	const float       half_scan ( 0.5 * nrecs * integ );

	/* the reference mjd inside the correlator */
	const double      mjd_ref ( shdr.mjd_ref );

	/* seek in file */
	rewind (fp);
	ltaseek (fp, linfo.stab[iscan].start_rec + linfo.srecs, recl);

	/* offsets */
	const int time_off  ( linfo.lhdr.time_off );
	const int time_size ( linfo.lhdr.time_size );

	if ( time_size != 8 ) throw std::runtime_error("type of time not recognized.");

	/* space to read */
	vb_type dbuf ( recl );

	/* read record from file */
	if ( fread ( dbuf.data(), recl, 1, fp ) != 1 ) 
		throw std::runtime_error ("Incomplete read while reading scan properties.");

	double *tref = reinterpret_cast<double*>( dbuf.data() + time_off );

	const double secs_since_ref ( *tref );

	/* this MJD is in UTC */
	const double mjd ( mjd_ref + ((secs_since_ref+half_scan)/86400.0) );

	/* get coordinates from a table */
	const sources::coords_t coord = sources::table.at ( source );
	const double ra  ( coord.first );
	const double dec ( coord.second );

	return scan_t ( source, ra, dec, mjd );
}

void LTA::get_baseline (int iscan, const antname_t& ant1, const antname_t& ant2, vc_type& inout) {
	/* sanity check */
	if (iscan >= linfo.scans) throw std::runtime_error ("Chosen scan out of range.");

	/* nrecs records for this scan */
	int nrecs  = linfo.stab[iscan].recs;

	/* get baseline index */
	int baseline_rr    = -1;
	int baseline_rl    = -1;
	int baseline_lr    = -1;
	int baseline_ll    = -1;

	for (int ibl = 0; ibl < nbaselines; ibl++) {

		const auto& ib = baselines[ibl];

		int match = 1;
		for (int im = 0; im < 3; im++) {
			if (ib.ant1[im] != ant1[im] || ib.ant2[im] != ant2[im]) {
				match = -1;
				break;
			}
		}
		//std::cout << "inputs=" << ant1 << " " << ant2 << " testing=" << ib.ant1 << " " << ib.ant2 << " match=" << match << std::endl;

		if ( match > 0) {
			if ( ib.band1 == 0 && ib.band2 == 0 ) {
				baseline_rr = ibl;
			}
			else if ( ib.band1 == 0 && ib.band2 == 1 ) {
				baseline_rl = ibl;
			}
			else if ( ib.band1 == 1 && ib.band2 == 0 ) {
				baseline_lr = ibl;
			}
			else if ( ib.band1 == 1 && ib.band2 == 1 ) {
				baseline_ll = ibl;
			}
		}
	}
	if ( baseline_rr < 0 || baseline_rl < 0 || baseline_lr < 0 || baseline_ll < 0 ) {
		throw std::runtime_error("did not find baseline");
	}

	int  data_off  = linfo.lhdr.data_off;
	int  data_size = linfo.lhdr.data_size;

	/* inout should be nrecs*nchannels*4 complex64 */
	/* there are four correlations */
	inout.assign ( 4 * nrecs * nchannels, complex_type(0.0, 0.0));
	/* layout = (time, channel, corr) */

	/* seek in file */
	rewind (fp);
	ltaseek (fp, linfo.stab[iscan].start_rec + linfo.srecs, recl);

	/* space to read */
	vb_type dbuf ( recl );

	for (int irec = 0; irec < nrecs; irec++) {

		/* read record from file */
		if ( fread ( dbuf.data(), recl, 1, fp ) != 1 ) 
			throw std::runtime_error ("Incomplete read while time averaging.");

		/* data offset, cast as float */
		float *ddata = reinterpret_cast<float*>(dbuf.data() + data_off);

		for (int ichan = 0; ichan < nchannels; ichan++) {

				/* get data */
				/* index in (baseline, channel) complex<float> */
				const int _rr      = ichan + nchannels * baseline_rr;
				const int _rl      = ichan + nchannels * baseline_rl;
				const int _lr      = ichan + nchannels * baseline_lr;
				const int _ll      = ichan + nchannels * baseline_ll;

				const complex_type rr ( ddata[2*_rr], ddata[2*_rr + 1] );
				const complex_type rl ( ddata[2*_rl], ddata[2*_rl + 1] );
				const complex_type lr ( ddata[2*_lr], ddata[2*_lr + 1] );
				const complex_type ll ( ddata[2*_ll], ddata[2*_ll + 1] );

				/* layout = (time, channel, corr) */
				/* (irec, ichan, corr) */
				inout [ 0 + 4*ichan + 4*nchannels*irec ] = rr;
				inout [ 1 + 4*ichan + 4*nchannels*irec ] = rl;
				inout [ 2 + 4*ichan + 4*nchannels*irec ] = lr;
				inout [ 3 + 4*ichan + 4*nchannels*irec ] = ll;

		} // channel

	} /* rec loop */

}
