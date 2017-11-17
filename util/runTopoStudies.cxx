//std/stl
#include <cstdlib>
#include <string>
#include <iostream>
using namespace std;

//ROOT
#include "TChain.h"
#include "TSystem.h"

//xAOD
#include "xAODRootAccess/Init.h"

//my
#include "TopoStudies/ChainHelper.h"
#include "TopoStudies/TopoTupler.h"

void help()
{
    cout << "___________________________________________________" << endl;
    cout << " ~TopoStudies~" << endl;
    cout << "" << endl;
    cout << " Usage: $ ./runTopoStudies -f [filelist/dir/file] [options]" << endl;
    cout << " Options:" << endl;
    cout << "   -f          input filelist,dir, or file [required]" << endl;
    cout << "   --em-filter" << endl;
    cout << "   --mu-filter" << endl;
    cout << "   -n          number of events to process (default: all)" << endl;
    cout << "   -d|--dbg    set debug level ON (default: OFF)" << endl;
    cout << "   -h|--help   print this help message" << endl;
    cout << "___________________________________________________" << endl;
}

int main(int argc, char** argv)
{
    string filelist = "";
    int nevents = -1;
    bool dbg = false;
    bool em_filter = false;
    bool mu_filter = false;


    int optin(1);
    while(optin < argc) {
        string opt = argv[optin];
        if      (opt == "-f") { filelist = argv[++optin]; }
        else if (opt == "--em-filter") { em_filter = true; }
        else if (opt == "--mu-filter") { mu_filter = true; }
        else if (opt == "-n") { nevents = atoi(argv[++optin]); }
        else if (opt == "-d" || opt == "--dbg") { dbg = true; }
        else if (opt == "-h" || opt == "--help") { help(); return 0; }
        else {
            cout << "Unknown command line argument : '" << opt << "'" << endl;
            return 1;
        }
        optin++;
    }; // while

    if(filelist == "") {
        cout << "You did not provide an input, exiting" << endl;
        return 1;
    }

    if(mu_filter && em_filter) {
        cout << "You requested both MU and EM filter. Can only do one. Not both." << endl;
        return 1;
    }

    cout << " ~TopoStudies~" << endl;

    TChain* chain = new TChain("CollectionTree");
    int file_err = ChainHelper::addInput(chain, filelist, true);
    if(file_err) return 1;
    xAOD::Init("TopoStudies");
    Long64_t n_entries = chain->GetEntries();
    chain->ls();
    if(nevents < 0) nevents = n_entries;
    topo::TopoTupler* ana = new topo::TopoTupler();
    ana->set_dbg(dbg);
    ana->set_filename(filelist);
    ana->set_mu_filter(mu_filter);
    ana->set_em_filter(em_filter);
    chain->Process(ana, "", nevents);


    delete chain;
    return 0;
}
