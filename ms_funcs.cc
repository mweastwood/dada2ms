/*
 * ms_funcs.cc
 * Functions to write to MS sub-tables plus some misc functions
 */

#include "ms_funcs.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>

// casacore headers
#include <casa/Arrays.h>
#include <tables/Tables.h>
#include <ms/MeasurementSets.h>

using namespace casa;

// Take a date,time string (utc) in the from YYYY-MM-DD-HH:MM:SS.S and an
// offset in seconds and return a MEpoch object.
MEpoch
str2MEpoch(const char *time, double offset)
{
    int yy, MM, dd, hh, mm;
    double ss;
    if (sscanf(time, "%d-%d-%d-%d:%d:%lf", &yy, &MM, &dd, &hh, &mm, &ss) != 6)
        throw std::runtime_error("Invalid String enountered in str2MEpoch()");
    return MEpoch(MVEpoch(MVTime(Time(yy, MM, dd, hh, mm, ss)+offset)), MEpoch::UTC);
}

MDirection
getZenith(const MPosition &pos, const MEpoch &epoch)
{
    MDirection up(Quantity(0.0, "deg"), Quantity(90.0, "deg"), MDirection::AZEL);
    MeasFrame mf(epoch, pos, up);
    MVDirection vzen;
    if (mf.getJ2000(vzen))
        return  MDirection(vzen, MDirection::J2000);
    else
        throw std::runtime_error("Error geting zenith postion in getZenith()");

}

inline double
radians(double degrees)
{
    return degrees * C::pi / 180.0;
}

// Get the WGS84 sea level from the center of the earth at the give latitude.
// May be an over simplification.
double
seaLevel(double latitude)
{
    const double majorAxis = 6378137.0; // Equator
    const double minorAxis = 6356752.3142; // Poles

    double rLat = radians(latitude);
    double x = majorAxis * cos(rLat);
    double y = minorAxis * sin(rLat);
    return hypot(x, y);
}

// Read a text file of 3 numbers per line
// and store as a 3xNANT Matrix. I suppose a Vector of Vectors
// would be more ideologically correct.
Matrix<Double>
readAnts(const char *filename, int nAnt)
{
    ifstream inf(filename);
    Matrix<Double> antPos(3, nAnt);
    for (int i=0; i<nAnt; i++) {
        inf >> antPos(0,i) >> antPos(1,i) >> antPos(2,i);
    }
    return antPos;
}

// Take a set of antenna positions and return the set of baselines
Matrix<Double>
zenithUVWs(Matrix<Double> antPos)
{
    int nAnt = antPos.ncolumn();
    Matrix<Double> uvw(3, (nAnt+1)*nAnt/2);
    int c=0;
    for (int i=0; i<nAnt; i++) {
        for (int j=i; j<nAnt; j++) {
        uvw(0,c) = antPos(0,i) - antPos(0,j);
        uvw(1,c) = antPos(1,i) - antPos(1,j);
        uvw(2,c) = antPos(2,i) - antPos(2,j);
        c++;
        }
    }
    return uvw;
}

// Calculate ITRF positions for antennas.
// antPos is offsets in meters from the array lon,lat,alt.
// Done by putting the antenna offset from (0,0) lon,lat and rotating
// to the array postion.
// We then have to jump through a few casacore hoops to convert to ITRF.
Matrix<Double>
itrfAnts(Matrix<Double> antPos, double longitude, double latitude, double altitude)
{
    int nAnt = antPos.ncolumn();
    // Origin is the centre of the earth.
    // X is towards 0 Lat, O Lon.
    // Y is towards 0 Lat, 90E Lon.
    // Z is towards North Pole.
    Matrix<Double> rLat = Rot3D(1, -radians(latitude)); // Negation to rotate anti-clockwise. i.e. towards porth pole
    Matrix<Double> rLon = Rot3D(2, radians(longitude));
    double seaLev = seaLevel(latitude);
    MPosition::Convert wgs2itrf(MPosition::WGS84, MPosition::ITRF); // conversion machine
    Matrix<Double> itrf(3,nAnt);

    Vector<Double> currPos(3), rotPos(3), pos(3);
    for (int i=0; i<nAnt; i++) {
        // Ant offsets are E,N,Alt. We want XYZ as defined above
        currPos(0) = antPos(2,i);
        currPos(1) = antPos(0,i);
        currPos(2) = antPos(1,i);

        currPos(0) += seaLev + altitude;
        rotPos = product(rLon, product(rLat, currPos));
        MVPosition posVectWGS(rotPos(0), rotPos(1), rotPos(2));
        // MPosition doesn't have a cartesian vector constructor, need lenght, angles.
        // WGS84 height is relative to sea level so subtract it out
        MPosition posWGS(posVectWGS.getLength("m") - Quantity(seaLev,"m"), posVectWGS.getAngle(), MPosition::WGS84);
        MPosition posITRF = wgs2itrf(posWGS);
        Vector<Double> posVectITRF = posITRF.get("m").getValue();
        itrf(0,i) = posVectITRF(0);
        itrf(1,i) = posVectITRF(1);
        itrf(2,i) = posVectITRF(2);
    }
    return itrf;
}


// The rest of the functions are to fill the extension tables.

int
fillAntTab(MSAntenna &ant, int nAnt, Matrix<Double> itrfPos)
{
    ant.addRow(nAnt);
    MSAntennaColumns cols(ant);
    int antNum;
    for (antNum=0; antNum<nAnt; antNum++) {
        std::stringstream outStr;
        outStr << ANT_NAME << std::setw(3) << std::setfill('0') << antNum+1;
        cols.name().put(antNum, outStr.str());
        cols.station().put(antNum, STATION_NAME);
        cols.type().put(antNum, ANTENNA_TYPE);
        cols.mount().put(antNum, ANTENNA_MOUNT);
        cols.offset().put(antNum, Vector<Double>(3, 0));
        cols.dishDiameter().put(antNum, ANTENNA_DISH_DIAMETER);
    }
    cols.position().putColumn(itrfPos);
    return antNum;
}

int
fillFeedTab(MSFeed &feed, int nAnt)
{
    feed.addRow(nAnt);
    MSFeedColumns cols(feed);

    Vector<Double> position(3, 0.0);
    Matrix<Double> beamOffset(2, 2, 0.0);
    Vector<String> polType(2);
    polType(0) = "X";
    polType(1) = "Y";
    Matrix<Complex> polResponse(2, 2, 0.0);
    polResponse(0, 0) = 1;
    polResponse(1, 1) = 1;
    Vector<Double> receptorAngle(2, 0.0);
    int antNum;
    for (antNum=0; antNum<nAnt; antNum++) {
        cols.position().put(antNum, position);
        cols.beamOffset().put(antNum, beamOffset);
        cols.polarizationType().put(antNum, polType);
        cols.polResponse().put(antNum, polResponse);
        cols.receptorAngle().put(antNum, receptorAngle);
        cols.antennaId().put(antNum, antNum);
        cols.beamId().put(antNum, -1);
        cols.feedId().put(antNum, 0);
        cols.interval().put(antNum, 1e30);
        cols.numReceptors().put(antNum, 2);
        cols.spectralWindowId().put(antNum, -1);
        cols.time().put(antNum, 0);
    }

    return antNum;
}

int
fillFieldTab(MSField &field, const MDirection *dir)
{
    MSFieldColumns cols(field);
    Vector<MDirection> zenith;
    if (dir == NULL) {
        cols.setDirectionRef(MDirection::AZEL);
        zenith = Vector<MDirection>(1, MDirection(Quantity(0.0, "deg"), Quantity(90.0, "deg"), MDirection::AZEL));
    } else {
        zenith = Vector<MDirection>(1, *dir);
    }
    field.addRow();
    cols.name().put(0, FIELD_NAME);
    cols.delayDirMeasCol().put(0, zenith);
    cols.phaseDirMeasCol().put(0, zenith);
    cols.referenceDirMeasCol().put(0, zenith);
    cols.sourceId().put(0, 0);

    return 0;
}

int
addField(MSField &field, const String &name, MDirection *dir)
{
	MSFieldColumns cols(field);
	int nrow = field.nrow();

	Vector<MDirection> dirVector;
	if (dir == NULL) {
	    cols.setDirectionRef(MDirection::AZEL);
	    dirVector = Vector<MDirection>(1, MDirection(Quantity(0.0, "deg"), Quantity(90.0, "deg"), MDirection::AZEL));
	} else {
	    dirVector = Vector<MDirection>(1, *dir);
	}
	field.addRow();
	cols.name().put(nrow, name);
	cols.delayDirMeasCol().put(nrow, dirVector);
	cols.phaseDirMeasCol().put(nrow, dirVector);
	cols.referenceDirMeasCol().put(nrow, dirVector);
	cols.sourceId().put(0, 0);

	return 0;
}

int
fillObservationTab(MSObservation &observation, Double startTime, Double finishTime)
{
    observation.addRow();
    MSObservationColumns observationCols(observation);
    Vector<Double> timeRange(2);
    timeRange[0] = startTime;
    timeRange[1] = finishTime;
    observationCols.timeRange().put(0, timeRange);
    observationCols.observer().put(0, OBSERVER);
    observationCols.project().put(0, PROJECT);
    observationCols.telescopeName().put(0, TELESCOPE_NAME);

    return 0;
}

int
fillPointingTab(MSPointing &pointing, int  const nAnt, Double time, const MDirection *dir=NULL)
{
    MSPointingColumns cols(pointing);
    Vector<MDirection> zenith;
    if (dir == NULL) {
        cols.setDirectionRef(MDirection::AZEL);
        zenith = Vector<MDirection>(1, MDirection(Quantity(0.0, "deg"), Quantity(90.0, "deg"), MDirection::AZEL));
    } else {
        zenith = Vector<MDirection>(1, *dir);
    }
    pointing.addRow(nAnt);
    int antNum;
    for (antNum=0; antNum<nAnt; antNum++) {
        cols.directionMeasCol().put(antNum, zenith);
        cols.antennaId().put(antNum, antNum);
        cols.interval().put(antNum, 1e30);
        cols.numPoly().put(antNum, 0);
        cols.targetMeasCol().put(antNum, zenith);
        cols.time().put(antNum, 0);
        cols.timeOrigin().put(antNum, time);
        cols.tracking().put(antNum, 0);
    }

    return antNum;
}

int
fillPolarizationTab(MSPolarization &polarization)
{
    polarization.addRow();
    MSPolarizationColumns cols(polarization);
    cols.numCorr().put(0, 4);

    Vector<Int> corrType(4);
    corrType[0] = Stokes::XX;
    corrType[1] = Stokes::XY;
    corrType[2] = Stokes::YX;
    corrType[3] = Stokes::YY;
    cols.corrType().put(0, corrType);

    Matrix<Int> corrProduct(2, 4);
    corrProduct(0, 0) = 0;
    corrProduct(1, 0) = 0;
    corrProduct(0, 1) = 0;
    corrProduct(1, 1) = 1;
    corrProduct(0, 2) = 1;
    corrProduct(1, 2) = 0;
    corrProduct(0, 3) = 1;
    corrProduct(1, 3) = 1;
    cols.corrProduct().put(0, corrProduct);

    return 0;
}

int
fillProcessorTab(MSProcessor &processor)
{
    processor.addRow();
    MSProcessorColumns processorCols(processor);
    processorCols.type().put(0, "CORRELATOR");
    processorCols.subType().put(0, CORRELATOR_NAME);

    return 0;
}

int
fillSpWindowTab(MSSpectralWindow &spw, int nFreq, double cFreq, double bw)
{
	int currRow = spw.nrow();
    spw.addRow();
    MSSpWindowColumns cols(spw);

    double refFreq = cFreq - bw / 2;
    double chanBW = bw / nFreq;
    Vector<Double> chanFreq(nFreq);
    for (int i=0; i<nFreq; i++) {
        chanFreq[i] = refFreq + (i+0.5) * chanBW;
    }
    Vector<Double> vChanBW(nFreq, chanBW);
    cols.measFreqRef().put(currRow, 1);
    cols.chanFreq().put(currRow, chanFreq);
    cols.refFrequency().put(currRow, refFreq);
    cols.chanWidth().put(currRow, vChanBW);
    cols.effectiveBW().put(currRow, vChanBW);
    cols.resolution().put(currRow, vChanBW);
    cols.freqGroupName().put(currRow, "Group 1");
    std::stringstream s;
    s << cFreq;
    cols.name().put(currRow, s.str().c_str());
    cols.netSideband().put(currRow, 1);
    cols.numChan().put(currRow, nFreq);
    cols.totalBandwidth().put(currRow, bw);

    return currRow;
}

int
addSourceTab(MeasurementSet &ms)
{
    TableDesc sourceDesc = MSSource::requiredTableDesc();
    MSSource::addColumnToDesc(sourceDesc, MSSource::TRANSITION, 1);
    MSSource::addColumnToDesc(sourceDesc, MSSource::REST_FREQUENCY, 1);
    MSSource::addColumnToDesc(sourceDesc, MSSource::SYSVEL, 1);
    SetupNewTable tabSetup(ms.sourceTableName(), sourceDesc, Table::New);
    ms.rwKeywordSet().defineTable(MS::keywordName(MS::SOURCE), Table(tabSetup));
    ms.initRefs();

    return 0;
}

int
fillSourceTab(MSSource &source, Double startTime, Double finishTime, const MDirection *dir)
{
    MSSourceColumns cols(source);
    MDirection zenith;
    if (dir == NULL) {
        cols.setDirectionRef(MDirection::AZEL);
        zenith = MDirection(Quantity(0.0, "deg"), Quantity(90.0, "deg"), MDirection::AZEL);
    } else {
        zenith = *dir;
    }
    source.addRow();
    cols.sourceId().put(0, 0);
    cols.time().put(0, (finishTime + startTime) / 2);
    cols.interval().put(0, finishTime - startTime);
    cols.spectralWindowId().put(0, -1);
    cols.numLines().put(0, 0);
    cols.name().put(0, FIELD_NAME);
    cols.calibrationGroup().put(0, 0);
    cols.code().put(0, "");
    cols.directionMeas().put(0, zenith);
    cols.properMotion().put(0, Vector<Double>(2, 0));

    return 0;
}

int
updateSourceTab(MSSource &source, Double startTime, Double finishTime)
{
    MSSourceColumns cols(source);
    cols.time().put(0, (finishTime + startTime) / 2);
    cols.interval().put(0, finishTime - startTime);

    return 0;
}

int
updateObservationTab(MSObservation &observation, Double startTime, Double finishTime)
{
    MSObservationColumns observationCols(observation);
    Vector<Double> timeRange(2);
    timeRange[0] = startTime;
    timeRange[1] = finishTime;
    observationCols.timeRange().put(0, timeRange);

    return 0;
}

void
boolArray2charVector(Array<Bool> &boolArr, std::vector<char> &charVec)
{
	if (boolArr.shape().product() != charVec.size())
		throw std::length_error("array length mismatch in ms_funcs::boolArray2charVector()");
	int i = 0;
	for (Array<Bool>::iterator iter=boolArr.begin(); iter != boolArr.end(); ++iter)
		charVec[i++] = static_cast<char>(*iter);
}

void
charVector2boolArray(std::vector<char> &charVec, Array<Bool> &boolArr)
{
	if (boolArr.shape().product() != charVec.size())
		throw std::length_error("array length mismatch in ms_funcs::boolArray2charVector()");
	int i = 0;
	for (Array<Bool>::iterator iter=boolArr.begin(); iter != boolArr.end(); ++iter)
		*iter = static_cast<Bool>(charVec[i++]);
}

void
readCalTable(const char *calName, std::vector<std::complex<float> > &gain, std::vector<char> &flag)
{
    const Table cal(calName, Table::Old);
    const int nrow = cal.nrow();
    ROArrayColumn<Complex> gain_col(cal, "CPARAM");
    IPosition cal_shape = gain_col.shape(0); // shape of cell 0
    cal_shape.append(IPosition(1, nrow));    // shape of column data
    int num_gains = cal_shape.product();
    gain.resize(num_gains);
    flag.resize(num_gains);
    // Use same data as vector, data ordering must be the same
    Array<Complex> arr_gain(cal_shape, gain.data(), SHARE);
    gain_col.getColumn(arr_gain);
    Array<Bool> arr_flag = ROArrayColumn<Bool>(cal, "FLAG").getColumn();
    boolArray2charVector(arr_flag, flag);

    // Check for "simple" table, one row per antenna
    ROScalarColumn<Int> antenna1(cal, "ANTENNA1");
    for (int i=0; i<nrow; i++) {
        if(antenna1.get(i) != i)
		throw std::length_error("Cal table not expected shape (one row per ant)");
    }
}

