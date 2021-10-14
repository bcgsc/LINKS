#ifndef LINKS_INPUTPARSER_HPP
#define LINKS_INPUTPARSER_HPP

#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>

//Globals
#define BASE_TEN 10
std::string progname = "LINKS";
std::string version = "2.0";

class InputParser {
    private:
    std::vector <std::string> tokens;
    int argc;
    int c;
    int optindex = 0;
	int help = 0;
    char* end = nullptr;

    std::vector<uint32_t> splitDistanceInput(std::string input){
        std::vector<uint32_t> distances; 

        std::stringstream ss(input);

        while( ss.good() )
        {
            std::string substr;
            getline( ss, substr, ',' );
            distances.push_back(static_cast<uint32_t>(std::stoul(substr) ));
        }
        return distances;
    }

    std::vector<uint16_t> splitStepSizesInput(std::string input){
        std::vector<uint16_t> step_sizes; 

        std::stringstream ss(input);

        while( ss.good() )
        {
            std::string substr;
            getline( ss, substr, ',' );
            step_sizes.push_back(static_cast<uint16_t>(std::stoul(substr) ));
        }
        return step_sizes;
    }

    std::queue<std::string> parseFofInput(std::string fofFile){
        std::queue<std::string> longReads; 

        std::ifstream infile(fofFile);
        std::string readFile;
        while(infile >> readFile){
            longReads.push(readFile);
            std::cout << longReads.back() << std::endl;
        }

        return longReads;
    }

    public:
    std::string assemblyFile = "";
    std::string fofFile = "";
    bool arguments_satisfied = true;
    std::queue<std::string> longReads;
    std::string distances_text = "";
    std::vector<uint32_t> distances = {4000};
    std::string step_sizes_text = "";
    std::vector<uint16_t> step_sizes = {2};
    //uint64_t distances = 4000;
    uint64_t k = 15;
    bool verbose = false;
    uint64_t minLinks = 5;
    uint64_t minSize = 500;
    float maxLinkRatio = 0.3;
    //uint64_t step = 2;
    // Added for MPET
    uint64_t readLength;         // MPET
    float insertStdev = 0.1;      // MPET (Need to adjust to a wider-range of distances when dealing with MPET) 
    std::string baseName = "";   // When set, this will override the MPET-induced changes on -e
    uint64_t offset = 0;
    std::string bfFile = "";
    float fpr = 0.001;
    bool bfOff = false;
    uint thread = 4;

    static void
    printUsage() {
        std::cout << "Usage:  " << progname << " " << version << "\n"
                << "  -f  sequences to scaffold (Multi-FASTA format, required)\n"
                    "  -s  file-of-filenames, full path to long sequence reads or MPET pairs [see below] (Multi-FASTA/fastq format, required)\n"
                    "  -m  MPET reads (default -m 1 = yes, default = no, optional)\n"
                    "  \t! DO NOT SET IF NOT USING MPET. WHEN SET, LINKS WILL EXPECT A SPECIAL FORMAT UNDER -s\n"
                    "  \t! Paired MPET reads in their original outward orientation <- -> must be separated by \":\"\n"
                    "  \t  >template_name\n\t  ACGACACTATGCATAAGCAGACGAGCAGCGACGCAGCACG:ATATATAGCGCACGACGCAGCACAGCAGCAGACGAC\n"
                    "  -d  distance between k-mer pairs (ie. target distances to re-scaffold on. default -d 4000, optional)\n"
                    "  \tMultiple distances are separated by comma. eg. -d 500,1000,2000,3000\n"
                    "  -k  k-mer value (default -k 15, optional)\n"
                    "  -t  step of sliding window when extracting k-mer pairs from long reads (default -t 2, optional)\n"
                    "  \tMultiple steps are separated by comma. eg. -t 10,5\n"
                    "  -o  offset position for extracting k-mer pairs (default -o 0, optional)\n"
                    "  -e  error (%) allowed on -d distance   e.g. -e 0.1  == distance +/- 10% (default -e 0.1, optional)\n"
                    "  -l  minimum number of links (k-mer pairs) to compute scaffold (default -l 5, optional)\n"
                    "  -a  maximum link ratio between two best contig pairs (default -a 0.3, optional)\n"
                    "  \t *higher values lead to least accurate scaffolding*\n"
                    "  -z  minimum contig length to consider for scaffolding (default -z 500, optional)\n"
                    "  -b  base name for your output files (optional)\n"
                    "  -r  Bloom filter input file for sequences supplied in -s (optional, if none provided will output to .bloom)\n"
                    "  \t NOTE: BLOOM FILTER MUST BE DERIVED FROM THE SAME FILE SUPPLIED IN -f WITH SAME -k VALUE\n"
                    "  \t IF YOU DO NOT SUPPLY A BLOOM FILTER, ONE WILL BE CREATED (.bloom)\n"
                    "  -p  Bloom filter false positive rate (default -p 0.001, optional; increase to prevent memory allocation errors)\n"
                    "  -x  Turn off Bloom filter functionality (-x 1 = yes, default = no, optional)\n"
                    "  -v  Runs in verbose mode (-v 1 = yes, default = no, optional)\n";
                    
                    //  nError: Missing mandatory options -f and -s.\n\n";
    }
    void
    printOpts() {
        // Command to test all input parsing
        // -f inputFasta.fa -s fofLongReads.txt -m 10 -d 20 -k 30 -t 40 -o 50 -e 60.45 -l 70 -a 1.345 -z 80 -b mybase -r myBloomFilter -p 0.003 -x 34 -v

        std::cout   << "  -f " << assemblyFile << "\n" 
                    << "  -s " << fofFile << "\n" 
                    << "  -m " << readLength << "\n" 
                    //<< "  -d " << distances.str() << "\n" TODO: print distances
                    << "  -k " << k << "\n"
                    //<< "  -t " << step << "\n"
                    << "  -o " << offset << "\n"
                    << "  -e " << insertStdev << "\n"
                    << "  -l " << minLinks << "\n"
                    << "  -a " << maxLinkRatio << "\n"
                    << "  -z " << minSize << "\n"
                    << "  -b " << baseName << "\n"
                    << "  -r " << bfFile << "\n"
                    << "  -p " << fpr << "\n"
                    << "  -x " << bfOff << "\n"
                    << "  -v " << verbose << "\n";
    }




    // std::string bfout = $base_name . ".bloom";
    // $base_name = $assemblyfile . ".scaff_s-" . $fofFile . "_d" . $distances . "_k" . $k . "_e" . $insert_stdev . "_l" . $min_links . "_a" . $max_link_ratio . "_z" . $min_size . "_t" . $step . "_o" . $offset . "_r-" . $bf_file . "_p" . $fpr . "_x" . $bfoff . "_m" . $readlength;

    InputParser() = default;

    InputParser (int &argc, char **argv) {
            this->argc = argc;
            static const struct option longopts[] = { { "help", no_argument, &help, 1 },
		                                      { nullptr, 0, nullptr, 0 } };
            while ((c = getopt_long(argc, argv, "f:s:m:d:k:t:j:o:e:l:a:z:b:r:p:x:v:", longopts, &optindex)) != -1) {
                switch (c) {
                case 0:
                    break;
                case 'f':
                    assemblyFile.assign(optarg);
                    break;
                case 's':
                    // full path for fof
                    fofFile.assign(optarg);
                    longReads = parseFofInput(fofFile);
                    break;
                case 'm':
                    readLength = strtoul(optarg, &end, BASE_TEN);
                    insertStdev = 0.5;
                    break;
                case 'd':
                    //distances.clear();
                    distances_text.assign(optarg); 
                    distances = splitDistanceInput(distances_text);
                    //distances = strtoul(optarg, &end, BASE_TEN);
                    //std::cout << optarg << "\n";
                    break;
                case 'k':
                    k = strtoul(optarg, &end, BASE_TEN);
                    break;
                case 't':
                    step_sizes_text.assign(optarg);
                    step_sizes = splitStepSizesInput(step_sizes_text);
                    //step = strtoul(optarg, &end, BASE_TEN);
                    break;
                case 'j':
                    thread = strtoul(optarg, &end, BASE_TEN);
                    break;
                case 'o':
                    offset = strtoul(optarg, &end, BASE_TEN);
                    break;
                case 'e':
                    insertStdev = strtof(optarg, &end);
                    break;
                case 'l':
                    minLinks = strtoul(optarg, &end, BASE_TEN);
                    break;
                case 'a':
                    maxLinkRatio = strtof(optarg, &end);
                    break;
                case 'z':
                    minSize = strtoul(optarg, &end, BASE_TEN);
                    break;
                case 'b':
                    baseName.assign(optarg);
                    break;
                case 'r':
                    bfFile.assign(optarg);
                    break;
                case 'p':
                    fpr = strtof(optarg, &end);
                    break;
                case 'x':
                    if(strtoul(optarg, &end, BASE_TEN) == 1){
                        bfOff = true;
                    }
                    break;
                case 'v':
                    if(strtoul(optarg, &end, BASE_TEN) == 1){
                        verbose = true;
                    }
                    break;
                default:
                    exit(EXIT_FAILURE);
                }
            }
            if(distances.size() > step_sizes.size()){
                uint last_step = step_sizes.back();
                for(int i = step_sizes.size(); i < distances.size(); i++){
                    step_sizes.push_back(last_step);
                }
            }
            if(distances.size() < step_sizes.size()){
                std::cerr << "\n Number of provided distances can't be lower than number of step sizes provided.\n";
            }
            if(baseName == ""){
                // copied from v1.8.7 -- add pid to basename
                baseName = assemblyFile + "scaff_s-" + fofFile + "_d" + distances_text + "_k" + std::to_string(k) + "_e" + std::to_string(insertStdev) + "_l" + std::to_string(minLinks) + "_a" + std::to_string(maxLinkRatio) + "_z" + std::to_string(minSize) + "_t" + step_sizes_text + "_o" + std::to_string(offset) + "_r-" + (bfOff ? "-" : bfFile) + "_p" + std::to_string(fpr) + "_x" + std::to_string(bfOff) + "_m" + std::to_string(readLength);
            }
            if(assemblyFile == "" || fofFile == ""){
                printUsage();
                arguments_satisfied = false;
            }
        }
};
#endif