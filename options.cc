/*
 * options.cc
 */

#include "options.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>

namespace dada2ms {

static const char * const default_config_file = "dada2ms.cfg";

options::options(int argc, char *argv[]) :
    utmzone(0),
    append(false),
    firstOnly(false),
    autosOnly(false),
    azel(false),
    addWtSpec(false),
    addSPW(false),
    applyCal(false),
    antsAreITRF(false),
    oneSPW(false),
    dataDescID(0),
    startScan(1),
    configFile(default_config_file)
{
    namespace po = boost::program_options;
    po::options_description poGeneric("Generic options");
    poGeneric.add_options()
        ("help,h", "produce this help message")
        ("config,c", po::value<std::string>(), "configuration file")
        //("config,c", po::value<std::string>(&configFile), "configuration file")
        ("append", po::bool_switch(&append), "append to MS (otherwise overwrite)")
        ("first", po::bool_switch(&firstOnly), "output first integration only")
        ("ints", po::value<std::string>(), "Integrations to take")
        //Not longer implemented ("autos", po::bool_switch(&autosOnly), "output auto-correlations only")
        ("azel", po::bool_switch(&azel), "Use AZEL reference for directions. "
                  "Otherwise J2000 is used.")
        ("wtspec", po::bool_switch(&addWtSpec), "create a WEIGHT_SPECTRUM column.")
        ("addspw", po::bool_switch(&addSPW), "create and use a new SPW for these data. Only used with --append.")
        ("onespw", po::bool_switch(&oneSPW), "concatenate all the dada files together so that only one SPW is created.")
        ("ddid", po::value<int>(&dataDescID), "use the specified pre-existing DATA_DESC_ID for these data. Only used with --append. Overridden by --addSPW. Default: 0")
        ("startscan", po::value<int>(&startScan), "use this value as the first scan/field value. Default: 1")
    ;
    po::options_description poConfig("Configuration options");
    poConfig.add_options()
        ("utmzone", po::value<int>(&utmzone), "utm zone of array (NAD83)")
        ("antfile", po::value<std::string>(), "antenna positions (NAD83 northing, easting, and elevation all in m)")
        // Currently disabled as untested ("itrfant", po::value<std::string>(), "antenna ITRF positions (m)")
        ("remap", po::value<std::string>(&remapFile), "remap lines as per file")
        ("cal", po::value<std::string>(&calTable), "Calibrate with CASA Bandpass table")
    ;
    po::options_description poHidden("Hidden options");
    poHidden.add_options()
        ("file-list", po::value<std::vector<std::string> >(), "file list")
    ;
    po::positional_options_description poPos;
    poPos.add("file-list", -1);

    po::options_description poCmdlineOptions;
    poCmdlineOptions.add(poGeneric).add(poConfig).add(poHidden);
    po::options_description poConfigFileOptions;
    poConfigFileOptions.add(poConfig).add(poHidden);
    po::options_description poVisible; // options to show with --help
    poVisible.add(poGeneric).add(poConfig);

    po::variables_map args;
    po::store(po::command_line_parser(argc, argv).options(poCmdlineOptions).positional(poPos).run(),  args);
    if (args.count("help")) {
        std::cout << "Fill Measurement Set from LEDA dada file." << std::endl << std::endl;
        std::cout << "Usage: " << argv[0] << " [options] <dada(s)> <MS>" << std::endl;
        std::cout << poVisible << std::endl;
        exit(EXIT_SUCCESS);
    }
    if (args.count("config"))
        configFile = args["config"].as<std::string>();
    std::ifstream ifs(configFile.c_str());
    if (ifs) {
        store(po::parse_config_file(ifs, poConfigFileOptions), args);
    }
    po::notify(args);
    if (args.count("file-list")) {
        dadaFile = args["file-list"].as<std::vector<std::string> >();
        msName = dadaFile.back();
        dadaFile.pop_back();
    }
    if (!args.count("file-list") || dadaFile.empty()) {
        std::cerr << "Error: at least one input and exactly one output required" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (args.count("itrfant")) {
        antsAreITRF = true;
        std::cerr << "Warning: ITRF antennas hasn't been tested in a while and probably broken." << std::endl;
    }
    if (args.count("antfile")) {
        if (antsAreITRF) {
            std::cerr << "Error: antfile and itrfant specified" << std::endl;
            exit(EXIT_FAILURE);
        }
    } else {
        if (!antsAreITRF) {
            std::cerr << "Error: antfile or itrfant required" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    if (args.count("ints"))
        integrations = split<int>(args["ints"].as<std::string>(), ',');

    if (args.count("cal"))
        applyCal = true;

    if (antsAreITRF)
        antFile = args["itrfant"].as<std::string>();
    else
        antFile = args["antfile"].as<std::string>();
}

template <typename T>
T StringToNumber(const std::string &Text)
{
    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

template <typename T>
std::vector<T> split(const std::string &s, char delim)
{
    std::vector<T> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(StringToNumber<T>(item));
    }
    return elems;
}

} // namespace dada2ms
