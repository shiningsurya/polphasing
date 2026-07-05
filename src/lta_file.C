/*
 * LTA interface impl
 */

#include "lta_file.hpp"
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
