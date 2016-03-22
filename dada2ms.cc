//
// Convert LWA-OVRO LEDA dada data to a Measurement Set.
//
// Most of the work is in reordering the data from its
// raw from of: [time][real/imag][freq][weirdBaselineOrder][pol]
//          to: [time][baseline][freq][pol][real/imag]
//
// It's supposed to be somewhat general; array position can be specified,
// antenna positions (relative to array position) are provided via a text file.
// Line order into correlator supplied as text file.
//
// Some deficiencies:
//   Not much in the way of error handling.
//   Does not calibrate weights in "cal" mode.
//
// g++ -O3 -I$CASACORE_INC_DIR -L$CASACORE_LIB_DIR -o dada2ms *.cc
//     -lcasa_casa -lcasa_measures -lcasa_ms -lcasa_tables -lcasa_scimath -lcasa_scimath_f
//     -lboost_program_options
//
// Stephen Bourke, Caltech
// March, 2014.
//

#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <complex>
#include <memory>

// casacore headers
#include <casa/Arrays.h>
#include <tables/Tables.h>
#include <ms/MeasurementSets.h>

#include "options.h"
#include "SortedDada.h"
#include "ms_funcs.h"
#include "MSUVWGenerator.h"

#include <iostream>
using namespace std;

using namespace casa;

int
main(int argc, char *argv[])
{
	dada2ms::options opts(argc, argv);

    std::vector<std::shared_ptr<dada::SortedDada> > dadas;
    for (int i = 0; i < opts.dadaFile.size(); ++i) {
        std::shared_ptr<dada::SortedDada> ptr(new dada::SortedDada(opts.dadaFile[i].c_str()));
        dadas.push_back(ptr);
    }
    std::shared_ptr<dada::SortedDada> dada = dadas[0];

    // Assigning to local variables to make code below more readable
    const int nAnt = dada->header.nAnt();
    int nTime = dada->header.nTime();
    int nFreq = dada->header.nFreq();
    const int nPol = dada->header.nPol();
    const int nCorr = dada->header.nCorr();
    const double intTime = dada->header.intTime(); // Integration time
    double cFreq = dada->header.cFreq();     // Center frequency
    double bw = dada->header.bandwidth();    // Bandwidth
    const double startTime = dada->header.startTimeMJD();
    const double finishTime = dada->header.finishTimeMJD();
    const int nBaseline = dada->header.nBaseline(); // Including autocorrelations

    if (opts.oneSPW) {
        // sum the number of frequency channels
        // average the center frequency (valid if subbands are contiguous)
        // sum the bandwidth
        for (int i = 1; i < dadas.size(); ++i) {
            nFreq += dadas[i]->header.nFreq();
            cFreq += dadas[i]->header.cFreq();
            bw += dadas[i]->header.bandwidth();
        }
        cFreq /= dadas.size();
    }

    if (opts.firstOnly) {
    	nTime = 1;
    }
    int outBaseline;
    if (opts.autosOnly) {
    	outBaseline = nAnt;
    } else {
    	outBaseline = nBaseline;
    }
    MeasurementSet ms;
    Matrix<Double> antPos = readAnts(opts.antFile.c_str(), nAnt);
    double array_longitude, array_latitude;
    utm2latlong(opts.utmzone,antPos(0,0)/1e3,antPos(1,0)/1e3,&array_latitude,&array_longitude);
    MPosition arrPos(Quantity(antPos(2,0), "m"), Quantity(array_longitude, "deg"), Quantity(array_latitude, "deg"), MPosition::WGS84);
    if (opts.append) {
        ms = MeasurementSet(opts.msName, Table::Update);
        updateObservationTab(ms.observation(), startTime, finishTime);
        updateSourceTab(ms.source(), startTime, finishTime);
        if (opts.addSPW) {
        	int setSPW = fillSpWindowTab(ms.spectralWindow(), nFreq, cFreq, bw);
        	opts.dataDescID = ms.dataDescription().nrow();
        	ms.dataDescription().addRow();
        	MSDataDescColumns cols(ms.dataDescription());
        	cols.spectralWindowId().put(opts.dataDescID, setSPW);
        }
    } else {
        // Create MS
        SetupNewTable newTab(opts.msName, MS::requiredTableDesc(), Table::New);
        ms = MeasurementSet(newTab);
        ms.createDefaultSubtables(Table::New);

        // Fill subtables
        if (opts.antsAreITRF) {
        	fillAntTab(ms.antenna(), nAnt, antPos);
        } else {
            Matrix<Double> itrf = itrfAnts(antPos, opts.utmzone);
           	fillAntTab(ms.antenna(), nAnt, itrf);
        }
        ms.dataDescription().addRow(); // One default row should do
        fillFeedTab(ms.feed(), nAnt);
        fillObservationTab(ms.observation(), startTime, finishTime);
        fillPolarizationTab(ms.polarization());
        fillProcessorTab(ms.processor());
        fillSpWindowTab(ms.spectralWindow(), nFreq, cFreq, bw);
        addSourceTab(ms);
        if (opts.azel) {
        	// NULL => zenith AZEL direction
        	// Single zenith field for full obs
        	addField(ms.field(), "Zenith", NULL);
        }
        fillSourceTab(ms.source(), startTime, finishTime, NULL);
        fillPointingTab(ms.pointing(), nAnt, startTime, NULL);

        // Add DATA column
        TiledShapeStMan dataStMan("dataHyperColumn", IPosition(2, nCorr, nFreq));
        ArrayColumnDesc<Complex> dataColDesc(MS::columnName(MS::DATA), "The data column", 2);
        ms.addColumn(dataColDesc, dataStMan);
        if (opts.addWtSpec) {
            TiledShapeStMan wtSpecStMan("weightSpecHyperColumn", IPosition(2, nCorr, nFreq));
            ArrayColumnDesc<Float> wtSpecColDesc(MS::columnName(MS::WEIGHT_SPECTRUM), "The weight spectrum column", 2);
            ms.addColumn(wtSpecColDesc, wtSpecStMan);
        }
    }

    MSColumns msCols(ms);
    int preexistingRows = ms.nrow();
    Int numFields = ms.field().nrow();
    Int firstScan = opts.startScan; // scan number of first scan in new data
    Int firstField = opts.startScan - 1;

    // Arrays for MS columns
    Vector<Double> timeVals;
    Vector<Int> ant1Vals(outBaseline);
    Vector<Int> ant2Vals(outBaseline);
    int bl = 0;
    for (int i=0; i<nAnt; ++i) {
    	for (int j=i; j<nAnt; ++j) {
    		ant1Vals[bl] = i;
    		ant2Vals[bl] = j;
    		++bl;
    	}
    }
    if (bl != outBaseline) {
    	throw std::logic_error("Too many baselines internal error");
    }

    // Optionally apply a simple (not time variable) CASA bandpass table
    std::vector<std::complex<float> > gain;
    std::vector<char> calFlag;
    if (opts.applyCal) {
        readCalTable(opts.calTable.c_str(), gain, calFlag);
        if (opts.oneSPW) {
            int offset = 0;
            for (int d = 0; d < dadas.size(); ++d) {
                std::vector<std::complex<float> > mygain;
                std::vector<char> mycalFlag;
                for (int antenna = 0; antenna < nAnt; ++antenna) {
                    for (int freq = 0; freq < dadas[d]->header.nFreq(); ++freq) {
                        for (int pol = 0; pol < nPol; ++pol) {
                            int idx = antenna*nFreq*nPol + (freq+offset)*nPol + pol;
                            int myidx = antenna*dadas[d]->header.nFreq()*nPol + freq*nPol + pol;
                            mygain[myidx] = gain[idx];
                            mycalFlag[myidx] = calFlag[idx];
                        }
                    }
                }
                offset += dadas[d]->header.nFreq();
                dadas[d]->applyGains(mygain, mycalFlag);
            }
        } else {
            dada->applyGains(gain, calFlag);
        }
    }

    if (!opts.remapFile.empty()) {
        dada->setLineMappingFromFile(opts.remapFile.c_str());
    }

    // We need to keep two copies of the flags due to the different storage (Bool vs char)
    std::vector<char> &charFlags = dada->rCurrentVisFlags();
    Cube<Bool> flag(nCorr, nFreq, outBaseline, false);

    // Arrays common to all integrations
    Matrix<Double> uvws;
    uvws = Matrix<Double>(3, outBaseline, 0); // All zeros initially
    Vector<Double> interval(outBaseline, intTime);
    Matrix<Float> unity2d(nCorr, outBaseline, 1.0);
    Cube<Float> unity3d(nCorr, nFreq, outBaseline, 1.0);
    Vector<Int> dataDescVals(outBaseline, opts.dataDescID);

    // opts.integrations holds the list of integrations to image
    if (opts.integrations.empty()) {
    	opts.integrations.resize(nTime);
    	for (int i=0; i<nTime; ++i) {
    	  	opts.integrations[i] = i;
    	}
    } else {
    	// Check integrations requested are valid
       	for (int i=0; i<opts.integrations.size(); ++i) {
       		if (opts.integrations[i] >= nTime) {
        		throw std::out_of_range("Invalid integration specified");
        	}
        }
    }

    // Add the integrations to the MS
    for (int i=0; i<opts.integrations.size(); ++i) {
    	int t = opts.integrations[i];
    	int currField;
    	if (opts.azel) {
    		currField = firstField;
    	} else {
    		currField = firstField + i;
    	}
    	Vector<Int> fieldVals(outBaseline, currField);
    	Double currTime = startTime + (t + 0.5) * intTime;
    	Vector<Double> timeVals(outBaseline, currTime);
    	Vector<Int> scanVals(outBaseline, firstScan + i);

    	ms.addRow(outBaseline);
        Array<Complex> data;
        if (opts.oneSPW) {
            // Note to Michael: the 3 indicates the number of dimensions in the array
            // (it is not the size of a 4th dimension)
            data = Array<Complex>(IPosition(3, nCorr, nFreq, outBaseline));
            int offset = 0;
            for (int d = 0; d < dadas.size(); ++d) {
                Array<Complex> mydata(IPosition(3, nCorr, dadas[d]->header.nFreq(), outBaseline), dadas[d]->rGetChunk(t).data(), SHARE);
                for (int baseline = 0; baseline < outBaseline; ++baseline) {
                    for (int freq = 0; freq < dadas[d]->header.nFreq(); ++freq) {
                        for (int corr = 0; corr < nCorr; ++corr) {
                            data(IPosition(3, corr, offset + freq, baseline)) = mydata(IPosition(3, corr, freq, baseline));
                        }
                    }
                }
                offset += dadas[d]->header.nFreq();
            }
        } else {
            data = Array<Complex>(IPosition(3, nCorr, nFreq, outBaseline), dada->rGetChunk(t).data(), SHARE);
        }
    	if (opts.applyCal) {
    		charVector2boolArray(charFlags, flag);
    	}
        // Create a Slicer for the current integration
        IPosition currIntStart(1, preexistingRows + i*outBaseline);
        IPosition currIntLength(1,outBaseline);
        IPosition currIntStride(1,1);
        Slicer currIntSlicer(currIntStart, currIntLength, currIntStride);
        msCols.uvw().putColumnRange(currIntSlicer, uvws);
        msCols.flag().putColumnRange(currIntSlicer, flag);
        msCols.weight().putColumnRange(currIntSlicer, unity2d);
        msCols.sigma().putColumnRange(currIntSlicer, unity2d);
        msCols.antenna1().putColumnRange(currIntSlicer, ant1Vals);
        msCols.antenna2().putColumnRange(currIntSlicer, ant2Vals);
        msCols.dataDescId().putColumnRange(currIntSlicer, dataDescVals);
        msCols.exposure().putColumnRange(currIntSlicer, interval);
        msCols.fieldId().putColumnRange(currIntSlicer, fieldVals);
        msCols.interval().putColumnRange(currIntSlicer, interval);
        msCols.scanNumber().putColumnRange(currIntSlicer, scanVals);
        msCols.time().putColumnRange(currIntSlicer, timeVals);
        msCols.timeCentroid().putColumnRange(currIntSlicer, timeVals);
        msCols.data().putColumnRange(currIntSlicer, data);
        if (opts.addWtSpec) {
            msCols.weightSpectrum().putColumnRange(currIntSlicer, unity3d);
        }
        if (!opts.azel && currField >= numFields) {
        	std::stringstream fieldName;
        	fieldName << "Zenith" << fixed << std::setprecision(2) << currTime;
        	MDirection dir = getZenith(arrPos, MEpoch(Quantity(currTime, "s"), MEpoch::UTC));
        	addField(ms.field(), fieldName.str(), &dir);
        }
    }

    // If we're only grabbing the autocorrelations, the uvws are already set
    // to zero so there is no need to calculate them.
    if (!opts.autosOnly) {
    	MSUVWGenerator uvwGen(msCols, MBaseline::J2000, Muvw::J2000);
    	Vector<Int> flds(nTime);
    	for (int i=0; i<nTime; i++) {
    		flds(i) = i;
    	}
    	uvwGen.make_uvws(flds);
    }

    return 0;
}
