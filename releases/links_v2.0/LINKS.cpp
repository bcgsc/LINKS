// using namespace btllib;
#include <string>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <cmath>
#include <regex>
#include <unordered_map>
#include <unordered_set>
#include <ctime>
#include <sstream>
// #include <pair>


#include "btllib/bloom_filter.hpp"
#include "btllib/seq_reader.hpp"

//Globals
#define BASE_TEN 10
std::string version = "2.0";

unsigned readFastAFastq_debug_counter_1 = 0;
unsigned readFastAFastq_debug_counter_2 = 0;
unsigned readFastAFastq_debug_counter_3 = 0;
unsigned readFastAFastq_debug_counter_4 = 0;
unsigned readFastAFastq_debug_counter_5 = 0;

unsigned readContigs_debug_counter_1 = 0;
unsigned readContigs_debug_counter_2 = 0;
unsigned readContigs_debug_counter_3 = 0;
unsigned readContigs_debug_counter_4 = 0;
unsigned readContigs_debug_counter_5 = 0;
unsigned readContigs_debug_counter_6 = 0;
unsigned readContigs_debug_counter_7 = 0;

unsigned pairContigs_debug_counter_1 = 0;
unsigned pairContigs_debug_counter_2 = 0;
unsigned pairContigs_debug_counter_3 = 0;
unsigned pairContigs_debug_counter_4 = 0;
unsigned pairContigs_debug_counter_5 = 0;
unsigned pairContigs_debug_counter_6 = 0;
unsigned pairContigs_debug_counter_7 = 0;
unsigned pairContigs_debug_counter_8 = 0;
unsigned pairContigs_debug_counter_9 = 0;
unsigned pairContigs_debug_counter_10 = 0;
unsigned pairContigs_debug_counter_11 = 0;
unsigned pairContigs_debug_counter_12 = 0;
unsigned pairContigs_debug_counter_13 = 0;
unsigned pairContigs_debug_counter_14 = 0;

unsigned pairContigs_debug_counter_50 = 0;
unsigned pairContigs_debug_counter_51 = 0;
unsigned pairContigs_debug_counter_52 = 0;
unsigned pairContigs_debug_counter_53 = 0;
unsigned pairContigs_debug_counter_54 = 0;
unsigned pairContigs_debug_counter_55 = 0;

class BT_IS {
    private:
    bool bt;
    int64_t is;
    // change this to a map of hash to distance
    // uint64_t distance;
    
    public:
    // BT_IS(const BT_IS&) = default;
    // BT_IS(BT_IS&&) = default;
    BT_IS() {
        this->bt = 0;
        this-> is = 0;
    }
    BT_IS(bool bt, int64_t is) {
        this->bt = bt;
        this-> is = is;
        // secondLayerReads.insert({hash, distance});
    }

    void setBT(int64_t bt) {
        this->bt = bt;
    }
    void setIS(int64_t is) {
        this->is = is;
    }
    bool getBT() {
        return this->bt;
    }
    int64_t getIS() {
        return this->is;
    }
};

class Gaps_Links {
    private:
    int64_t gaps;
    uint64_t links;
    std::string type;
    
    public:
    Gaps_Links(int64_t gap, uint64_t link) {
        this->gaps = gap;
        this->links = link;
        // secondLayerReads.insert({hash, distance});
    }

    Gaps_Links() {
        this->gaps = 0;
        this->links = 0;
        // secondLayerReads.insert({hash, distance});
    }
    
    Gaps_Links(std::string type) {
        this->type = type;
    }

    void addToGap(int64_t distance) {
        this->gaps += distance;
    }

    void incrementLinks() {
        this->links++;
    }

    int64_t getGaps() {
        return this->gaps;
    }
    uint64_t getLinks() {
        return this->links;
    }
};

class KmerInfo {
    private:
    std::string tig;
    uint64_t start;
    uint64_t end;
    uint64_t multiple;
    bool orient;
    
    public:
    KmerInfo(){
        this->tig = "";
        this->start = 0;
        this->end = 0;
        this->multiple = 1;
        this->orient = 0;
    }
    KmerInfo(uint64_t start, uint64_t end){
        this->tig = "";
        this->start = start;
        this->end = end;
        this->multiple = 1;
        this->orient = 0;
    }
    KmerInfo(std::string tig, uint64_t start, uint64_t end, bool orient){
        this->tig = tig;
        this->start = start;
        this->end = end;
        this->multiple = 1;
        this->orient = orient;
    }

    void setTig(std::string tig) {
        this->tig = tig;
    }

    void setStart(uint64_t start) {
        this->start = start;
    }

    void setEnd(uint64_t end) {
        this->end = end;
    }

    uint64_t getMultiple() {
        return this->multiple;
    }

    void incrementMultiple() {
        this->multiple += 1;
    }

    std::string getTig() {
        return this->tig;
    }

    uint64_t getStart() {
        return this->start;
    }

    uint64_t getEnd() {
        return this->end;
    }
    void setOrient(bool f){
        this->orient = f;
    }
    bool getOrient(){
        return this->orient;
    }
    
};

class InputParser {
    private:
    std::vector <std::string> tokens;
    int argc;
    int c;
    int optindex = 0;
	int help = 0;
    char* end = nullptr;

    std::vector<uint64_t> splitDistanceInput(std::string input){
        std::vector<uint64_t> distances; 
        std::istringstream ss(input); // Turn the string into a stream. 
        std::string tok; 

        while (ss >> tok) 
        {
            distances.push_back(static_cast<unsigned int>(std::stoul(tok)));
            std::cout << distances.back() << "\n";
        }

        std::sort(distances.begin(),distances.end());
 
        return distances;
    }

    public:
    std::string assemblyFile;
    std::string longFile;
    std::vector<uint64_t> distances = {4000};
    //uint64_t distances = 4000;
    uint64_t k = 15;
    bool verbose = false;
    uint64_t minLinks = 5;
    uint64_t minSize = 500;
    float maxLinkRatio = 0.3;
    uint64_t step = 2;
    // Added for MPET
    uint64_t readLength;         // MPET
    float insertStdev = 0.1;      // MPET (Need to adjust to a wider-range of distances when dealing with MPET) 
    std::string baseName;   // When set, this will override the MPET-induced changes on -e
    uint64_t offset = 0;
    std::string bfFile;
    float fpr = 0.001;
    uint64_t bfoff = 0;

    // std::string bfout = $base_name . ".bloom";
    // $base_name = $assemblyfile . ".scaff_s-" . $longfile . "_d" . $distances . "_k" . $k . "_e" . $insert_stdev . "_l" . $min_links . "_a" . $max_link_ratio . "_z" . $min_size . "_t" . $step . "_o" . $offset . "_r-" . $bf_file . "_p" . $fpr . "_x" . $bfoff . "_m" . $readlength;



    InputParser (int &argc, char **argv) {
            this->argc = argc;
            static const struct option longopts[] = { { "help", no_argument, &help, 1 },
		                                      { nullptr, 0, nullptr, 0 } };
            while ((c = getopt_long(argc, argv, "f:s:m:d:k:t:o:e:l:a:z:b:r:p:x:v", longopts, &optindex)) != -1) {
                switch (c) {
                case 0:
                    break;
                case 'f':
                    assemblyFile.assign(optarg);
                    break;
                case 's':
                    // full path for fof
                    longFile.assign(optarg);
                    break;
                case 'm':
                    readLength = strtoul(optarg, &end, BASE_TEN);
                    insertStdev = 0.5;
                    break;
                case 'd':
                    //distances.clear();
                    distances = splitDistanceInput(optarg);
                    //distances = strtoul(optarg, &end, BASE_TEN);
                    //std::cout << optarg << "\n";
                    break;
                case 'k':
                    k = strtoul(optarg, &end, BASE_TEN);
                    break;
                case 't':
                    step = strtoul(optarg, &end, BASE_TEN);
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
                    bfoff = strtoul(optarg, &end, BASE_TEN);
                    break;
                case 'v':
                    verbose = true;
                    break;
                default:
                    exit(EXIT_FAILURE);
                }
            }
        }
    // class 
    static void
    printErrorMsg(const std::string& progname, const std::string& msg) {
        std::cerr << progname << ": " << msg << "\nTry 'physlr-makebf --help' for more information.\n";
    }

    static void
    printUsage(const std::string& progname) {
        std::cout << "Usage:  " << progname << " " << version << "\n"
                << "  -f  sequences to scaffold (Multi-FASTA format, required)\n"
                    "  -s  file-of-filenames, full path to long sequence reads or MPET pairs [see below] (Multi-FASTA/fastq format, required)\n"
                    "  -m  MPET reads (default -m 1 = yes, default = no, optional)\n"
                    "  \t! DO NOT SET IF NOT USING MPET. WHEN SET, LINKS WILL EXPECT A SPECIAL FORMAT UNDER -s\n"
                    "  \t! Paired MPET reads in their original outward orientation <- -> must be separated by \":\"\n"
                    "  \t  >template_name\n\t  ACGACACTATGCATAAGCAGACGAGCAGCGACGCAGCACG:ATATATAGCGCACGACGCAGCACAGCAGCAGACGAC\n"
                    "  -d  distance between k-mer pairs (ie. target distances to re-scaffold on. default -d $distances, optional)\n"
                    "  \tMultiple distances are separated by comma. eg. -d 500,1000,2000,3000\n"
                    "  -k  k-mer value (default -k $k, optional)\n"
                    "  -t  step of sliding window when extracting k-mer pairs from long reads (default -t $step, optional)\n"
                    "  \tMultiple steps are separated by comma. eg. -t 10,5\n"
                    "  -o  offset position for extracting k-mer pairs (default -o $offset, optional)\n"
                    "  -e  error (%) allowed on -d distance   e.g. -e 0.1  == distance +/- 10% (default -e $insert_stdev, optional)\n"
                    "  -l  minimum number of links (k-mer pairs) to compute scaffold (default -l $min_links, optional)\n"
                    "  -a  maximum link ratio between two best contig pairs (default -a $max_link_ratio, optional)\n"
                    "  \t *higher values lead to least accurate scaffolding*\n"
                    "  -z  minimum contig length to consider for scaffolding (default -z $min_size, optional)\n"
                    "  -b  base name for your output files (optional)\n"
                    "  -r  Bloom filter input file for sequences supplied in -s (optional, if none provided will output to .bloom)\n"
                    "  \t NOTE: BLOOM FILTER MUST BE DERIVED FROM THE SAME FILE SUPPLIED IN -f WITH SAME -k VALUE\n"
                    "  \t IF YOU DO NOT SUPPLY A BLOOM FILTER, ONE WILL BE CREATED (.bloom)\n"
                    "  -p  Bloom filter false positive rate (default -p $fpr, optional; increase to prevent memory allocation errors)\n"
                    "  -x  Turn off Bloom filter functionality (-x 1 = yes, default = no, optional)\n"
                    "  -v  Runs in verbose mode (-v 1 = yes, default = no, optional)\n";
                    
                    //  nError: Missing mandatory options -f and -s.\n\n";
    }
    void
    printOpts() {
        // Command to test all input parsing
        // -f inputFasta.fa -s fofLongReads.txt -m 10 -d 20 -k 30 -t 40 -o 50 -e 60.45 -l 70 -a 1.345 -z 80 -b mybase -r myBloomFilter -p 0.003 -x 34 -v

        std::cout   << "  -f " << assemblyFile << "\n" 
                    << "  -s " << longFile << "\n" 
                    << "  -m " << readLength << "\n" 
                    //<< "  -d " << distances.str() << "\n" TODO: print distances
                    << "  -k " << k << "\n"
                    << "  -t " << step << "\n"
                    << "  -o " << offset << "\n"
                    << "  -e " << insertStdev << "\n"
                    << "  -l " << minLinks << "\n"
                    << "  -a " << maxLinkRatio << "\n"
                    << "  -z " << minSize << "\n"
                    << "  -b " << baseName << "\n"
                    << "  -r " << bfFile << "\n"
                    << "  -p " << fpr << "\n"
                    << "  -x " << bfoff << "\n"
                    << "  -v " << verbose << "\n";
    }
};

// Helper Methods
void sortErr(std::unordered_map<std::string, Gaps_Links>& M);
void sortInsertSize(std::unordered_map<uint64_t, uint64_t>& M);
void printBloomStats(btllib::KmerBloomFilter& bloom, std::ostream& os);
bool does_file_exist(std::string fileName);
btllib::KmerBloomFilter *makeBF(uint64_t bfElements, InputParser linksArgParser);
uint64_t getFileSize(std::string filename);
void readFastaFastq(
                const std::string file,
                const btllib::BloomFilter& bloom, 
                std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS>>& matePair,
                std::unordered_set<uint64_t>& mates, 
                const uint64_t distance,
                const uint64_t k,
                const uint64_t step);
void kmerize(std::string *seq,
                uint64_t frag_dist,
                uint64_t k,
                std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > *matePair,
                uint64_t step,
                std::string head,
                uint64_t endPosition,
                uint64_t initPosition,
                uint64_t pairct,
                btllib::KmerBloomFilter *bloom,
                uint64_t delta,
                uint64_t readLength);
void readContigs(
        std::string assemblyFile,
        std::unordered_map<uint64_t, KmerInfo>& trackAll,
        std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > matePair,
        std::unordered_set<uint64_t>& mates,
        std::unordered_map<std::string, uint64_t>& tigLength,
        uint64_t k,
        uint64_t minSize,
        unsigned hashFcts,
        uint64_t step);
void pairContigs(
    std::string longReadsFile,
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > matePair,
    std::unordered_map<uint64_t, KmerInfo>& trackAll,
    std::unordered_map<std::string, uint64_t> tigLength,
    std::string issues,
    std::string distribution,
    uint64_t totalPairs,
    std::string tigpair_checkpoint,
    std::string simplepair_checkpoint,
    bool verbose,
    float insertStdev);
int getDistance(
    uint64_t insert_size,
    uint64_t length_i,
    uint64_t start_i,
    uint64_t start_j);
void addToPairMap(
    int& isz,
    std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, Gaps_Links>>>& pair,
    int& distance,
    std::string kmer1_name,
    std::string kmer2_name,
    unsigned orient_enum
);
inline int getDistanceBin(int distance)
{
    return distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 10000;
}


int main(int argc, char** argv) {
    // Get todays date
    time_t now = time(0);
   
    // convert now to string form
    char* dt = ctime(&now);
    // Parse command line arguments
    InputParser linksArgParser = InputParser(argc, argv);

    //Set Bloom Filter element number based on the size of the assembly file (1 byte = 1 character)
    std::string assemblyPath = linksArgParser.assemblyFile;
    uint64_t bfElements = getFileSize(assemblyPath);

    // Checking validity of input assemble file
    if(bfElements == -1){
        std::cout << std::to_string(bfElements);
        std::cout << "Invalid file: " << linksArgParser.assemblyFile << " -- fatal\n";
        return -1;
    }

    // Naming output files
    if (linksArgParser.baseName == "") {
        linksArgParser.baseName =  linksArgParser.assemblyFile + ".scaff_s-" + 
                                    linksArgParser.longFile + "_d" + 
                                    //std::to_string(linksArgParser.distances.str()) + TODO:: print distance
                                    "_k" + std::to_string(linksArgParser.k) + 
                                    "_e" + std::to_string(linksArgParser.insertStdev) +
                                    "_l" + std::to_string(linksArgParser.minLinks) +
                                    "_a" + std::to_string(linksArgParser.maxLinkRatio) +
                                    "_z" + std::to_string(linksArgParser.minSize) +
                                    "_t" + std::to_string(linksArgParser.step) +
                                    "_o" + std::to_string(linksArgParser.offset) +
                                    "_r-"+ linksArgParser.bfFile +
                                    "_p" + std::to_string(linksArgParser.fpr) +
                                    "_x" + std::to_string(linksArgParser.bfoff) +
                                    "_m" + std::to_string(linksArgParser.readLength);
        pid_t pid_num = getppid();
        linksArgParser.baseName += "_pid" + std::to_string(pid_num);
    }
    std::cout << linksArgParser.baseName;
    std::string outlog = linksArgParser.baseName + ".log";
    std::string scaffold = linksArgParser.baseName + ".scaffolds";
    std::string issues = linksArgParser.baseName + ".pairing_issues";
    std::string distribution = linksArgParser.baseName + ".pairing_distribution.csv";
    std::string bfout = linksArgParser.baseName + ".bloom";
    std::string graph = linksArgParser.baseName + ".gv";
    std::string numnamecorr = linksArgParser.baseName + ".assembly_correspondence.tsv";
    std::string tigpair_checkpoint = linksArgParser.baseName + ".tigpair_checkpoint.tsv"; // add a checkpoint file, prevent re-running LINKS from scratch if crash
    std::string simplepair_checkpoint = linksArgParser.baseName + ".simplepair_checkpoint.tsv"; // add a checkpoint file, prevent re-running LINKS from scratch if crash
    std::string assemblyruninfo = "";
    std::string reading_tigbloom_message = "";

    if(freopen(outlog.c_str(), "w", stderr ) == NULL) {
        std::cout << "Can't write to " << outlog << " -- fatal\n";
        return -1;
    }
    //---------------------------
    // Initialization message
    std::string initMessage = "\nRunning: " + version + 
                                "\n-f " + linksArgParser.assemblyFile +
                                "\n-s " + linksArgParser.longFile + 
                                "\n-m " + std::to_string(linksArgParser.readLength) + 
                                //"\n-d " + std::to_string(linksArgParser.distances.str()) + TODO: print distance
                                "\n-k " + std::to_string(linksArgParser.k) + 
                                "\n-e " + std::to_string(linksArgParser.insertStdev) +
                                "\n-l " + std::to_string(linksArgParser.minLinks) +
                                "\n-a " + std::to_string(linksArgParser.maxLinkRatio) +
                                "\n-t " + std::to_string(linksArgParser.step) + 
                                "\n-o " + std::to_string(linksArgParser.offset) +
                                "\n-z " + std::to_string(linksArgParser.minSize) +
                                "\n-b " + linksArgParser.baseName +
                                "\n-r " + linksArgParser.bfFile +
                                "\n-p " + std::to_string(linksArgParser.fpr) +
                                "\n-x " + std::to_string(linksArgParser.bfoff) +
                                "\n\n----------------- Verifying files -----------------\n\n";
    
    std::cout << initMessage;
    std::cerr << initMessage;
    assemblyruninfo += initMessage;
    //-----------------------------
    btllib::KmerBloomFilter * myFilter;
    std::cout << "\n\n=>Before makeBF in c++ " + std::to_string(time(0)) + "\n";
    myFilter = makeBF(bfElements, linksArgParser);
    std::cout << "\n\n=>After makeBF in c++ " + std::to_string(time(0)) + "\n";
    unsigned hashFct = myFilter->get_hash_num();

    // // k-merize long reads
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > matePair;

    btllib::BloomFilter& filtering = myFilter->get_bloom_filter();

    // store second mates in a set
    std::unordered_set<uint64_t> mates;
    
    
    uint64_t totalpairs = 0;
    for(uint64_t dist : linksArgParser.distances){
        readFastaFastq(linksArgParser.longFile,filtering,matePair,mates,dist,linksArgParser.k,linksArgParser.step);
    }

    // debug by MURATHAN --->
    for (auto& it: matePair) {
        readFastAFastq_debug_counter_4 += it.second.size();
    }

    std::cout << "readFastAFastq_debug_counter_1 " + std::to_string(readFastAFastq_debug_counter_1) + "\n"; 
    std::cout << "readFastAFastq_debug_counter_2 " + std::to_string(readFastAFastq_debug_counter_2) + "\n"; 
    std::cout << "readFastAFastq_debug_counter_3 " + std::to_string(readFastAFastq_debug_counter_3) + "\n";
    std::cout << "readFastAFastq_debug_counter_4 " + std::to_string(readFastAFastq_debug_counter_4) + "\n";
    std::cout << "mates size: " + std::to_string(mates.size()) + "\n";
    std::cout << "matepair size: " + std::to_string(matePair.size()) + "\n";
    // debug by MURATHAN <---

    std::unordered_map<uint64_t, KmerInfo> trackAll;
    //std::unordered_map<uint64_t, KmerInfo> trackFor;
    //std::unordered_map<uint64_t, KmerInfo> trackRev;
    std::unordered_map<std::string, uint64_t> tigLength;
    // std::cout << "\n\n=>Reading sequence contigs (to scaffold), tracking k-mer positions :" << dt << "\n";
    // // Read contigs to find where the long read kmers belong in
    //std::cout << "innercounter: " << innercounter << std::endl; 
    //std::cout << "After reading long reads hit size: " << totalpairs << std::endl; 
    //std::cout << "\n\n=>Before readContigs c++ " + std::to_string(time(0)) + "\n";
    readContigs(linksArgParser.assemblyFile, trackAll, matePair, mates, tigLength, linksArgParser.k, linksArgParser.minSize, hashFct, linksArgParser.step);
    // mates.clear();
    //std::cout << "\n\n=>After readContigs c++ " + std::to_string(time(0)) + "\n";
  

    std::cout << "\nreadContigs_debug_counter_1: total kmer read from contig\n";
    std::cout << "readContigs_debug_counter_2: where forward has found in matepair\n";
    std::cout << "readContigs_debug_counter_1: where reverse hash found in matepair\n";
    std::cout << "readContigs_debug_counter_4: trackFor size\n";
    std::cout << "readContigs_debug_counter_5: trackRev size\n";

    std::cout << "readContigs_debug_counter_1 " + std::to_string(readContigs_debug_counter_1) + "\n";
    std::cout << "readContigs_debug_counter_2 " + std::to_string(readContigs_debug_counter_2) + "\n";
    std::cout << "readContigs_debug_counter_3 " + std::to_string(readContigs_debug_counter_3) + "\n";
    //std::cout << "readContigs_debug_counter_4 " + std::to_string(readContigs_debug_counter_4) + "\n";
    std::cout << "readContigs_debug_counter_4 " + std::to_string(readContigs_debug_counter_4) + "\n";
    std::cout << "readContigs_debug_counter_5 " + std::to_string(readContigs_debug_counter_5) + "\n";

    unsigned forward_kmers = 0, reverse_kmers = 0;
    for (auto& it: trackAll) {
        if(it.second.getOrient()){
            reverse_kmers++;
        }
        else{
            forward_kmers++;
        }
    }
    std::cout << "reverse_kmers " + std::to_string(reverse_kmers) + "\n";
    std::cout << "forward_kmers " + std::to_string(forward_kmers) + "\n";

    //std::cout << " The resulting trackAll map size is: " << trackAll.size() << "\n\n";
    //std::cout << " pairContigs Parameter List : \n\n";
    //std::cout << " 1- LongFile " << linksArgParser.longFile <<"\n";
    //std::cout << " 2- matePair Size " << matePair.size() <<"\n";
    //std::cout << " 3- trackAll size " << trackAll.size() <<"\n";
    //std::cout << " 4- trackFor size " << trackFor.size() <<"\n";
    //std::cout << " 5- trackRev size " << trackRev.size() <<"\n";
    //std::cout << " 6- tigLength size " << tigLength.size() <<"\n";
    //std::cout << " 7- issuesName " << issues <<"\n";
    //std::cout << " 8- distributionName " << distribution <<"\n";
    //std::cout << " 9- totalPairs " << totalpairs <<"\n";
    //std::cout << " 10- tigpair_checkpoint " << tigpair_checkpoint <<"\n";
    //std::cout << " 11- verbose " << linksArgParser.verbose <<"\n";
    //std::cout << " 12- distributionName " << linksArgParser.insertStdev <<"\n";
    //std::cout << "\n\n=>Before pairContigs c++ " + std::to_string(time(0)) + "\n";
    pairContigs(
        linksArgParser.longFile,
        matePair,
        trackAll,
        tigLength,
        issues,
        distribution,
        totalpairs,
        tigpair_checkpoint,
        simplepair_checkpoint,
        linksArgParser.verbose,
        linksArgParser.insertStdev);

    std::cout << "pairContigs_debug_counter_1: all matepair iteration\n";
    std::cout << "pairContigs_debug_counter_2: matepair is not seen\n";
    std::cout << "pairContigs_debug_counter_3: second mate is found in assembly\n";
    std::cout << "pairContigs_debug_counter_4: second mate is singleton \n";
    std::cout << "pairContigs_debug_counter_5: first mate is in assembly\n";
    std::cout << "pairContigs_debug_counter_6: first mate is in singleton\n";
    std::cout << "pairContigs_debug_counter_7: both contigs found in trackall\n";

    std::cout << "pairContigs_debug_counter_1 " + std::to_string(pairContigs_debug_counter_1) + "\n";
    std::cout << "pairContigs_debug_counter_2 " + std::to_string(pairContigs_debug_counter_2) + "\n";
    std::cout << "pairContigs_debug_counter_3 " + std::to_string(pairContigs_debug_counter_3) + "\n";
    std::cout << "pairContigs_debug_counter_4 " + std::to_string(pairContigs_debug_counter_4) + "\n";
    std::cout << "pairContigs_debug_counter_5 " + std::to_string(pairContigs_debug_counter_5) + "\n";
    std::cout << "pairContigs_debug_counter_6 " + std::to_string(pairContigs_debug_counter_6) + "\n";
    std::cout << "pairContigs_debug_counter_7 " + std::to_string(pairContigs_debug_counter_7) + "\n";
    std::cout << "pairContigs_debug_counter_8 " + std::to_string(pairContigs_debug_counter_8) + "\n";
    std::cout << "pairContigs_debug_counter_9 " + std::to_string(pairContigs_debug_counter_9) + "\n";
    std::cout << "pairContigs_debug_counter_10 " + std::to_string(pairContigs_debug_counter_10) + "\n";
    std::cout << "pairContigs_debug_counter_11 " + std::to_string(pairContigs_debug_counter_11) + "\n";
    std::cout << "pairContigs_debug_counter_12 " + std::to_string(pairContigs_debug_counter_12) + "\n";
    std::cout << "pairContigs_debug_counter_13 " + std::to_string(pairContigs_debug_counter_13) + "\n";
    std::cout << "pairContigs_debug_counter_14 " + std::to_string(pairContigs_debug_counter_14) + "\n";

    std::cout << "pairContigs_debug_counter_50 " + std::to_string(pairContigs_debug_counter_50) + "\n";
    std::cout << "pairContigs_debug_counter_51 " + std::to_string(pairContigs_debug_counter_51) + "\n";
    std::cout << "pairContigs_debug_counter_52 " + std::to_string(pairContigs_debug_counter_52) + "\n";
    std::cout << "pairContigs_debug_counter_53 " + std::to_string(pairContigs_debug_counter_53) + "\n";
    // --  std::cout << "\n\n=>After pairContigs c++ " + std::to_string(time(0)) + "\n";
    return 0;
}

btllib::KmerBloomFilter *makeBF(uint64_t bfElements, InputParser linksArgParser) {
    btllib::KmerBloomFilter * assemblyBF;
    if(linksArgParser.bfFile != "") {
        std::cout << "A Bloom filter was supplied (" << linksArgParser.bfFile << ") and will be used instead of building a new one from -f " << linksArgParser.assemblyFile << "\n";
        if(!does_file_exist(linksArgParser.bfFile)) {
            std::cout << "\nInvalid file: " << linksArgParser.bfFile <<  " -- fatal\n";
            exit;
        } else {
            std::cout << "Checking Bloom filter file " << linksArgParser.bfFile <<"...ok\n";
        }
        // std::cout << "Loading bloom filter of size " << getFileSize(linksArgParser.bfFile) << " from " << linksArgParser.bfFile << "\n";
        assemblyBF = new btllib::KmerBloomFilter(linksArgParser.bfFile);
    } else {
        uint64_t m = ceil((-1 * (double)bfElements * log(linksArgParser.fpr)) / (log(2) * log(2)));
        uint64_t rem = 64 - (m % 64);
        m = ((uint64_t)(m / 8) + 1) * 8;
        std::cout << "HASHES CALC: " << std::to_string(((double)m / bfElements)) << " second: " << std::to_string(((double)m / bfElements) * log(2)) << "\n";
        unsigned hashFct = floor(((double)m / bfElements) * log(2));
        std::cout << "- Number of bfElements: " << bfElements << "\n"
                    << "- Input file path: " << linksArgParser.bfFile << "\n"
                    << "- Input file: " << linksArgParser.assemblyFile << "\n"
                    << "- kmersize: " << linksArgParser.k << "\n"
                    << "- m: " << m << "\n"
                    << "- fpr: " << linksArgParser.fpr << "\n"
                    << "- hashFct: " << hashFct << "\n";

        std::string reading_tigbloom_message = "\n\n=>Reading contig/sequence assembly file : " + std::to_string(time(0)) + "\n";
        // assemblyruninfo += reading_tigbloom_message;
        // std::cout << "- Filter output file : " << outFileBf << "\n";
        std::cout << "- Filter output file : " << linksArgParser.k << "\n";
        assemblyBF = new btllib::KmerBloomFilter(m/8, hashFct, linksArgParser.k);
        btllib::SeqReader assemblyReader(linksArgParser.assemblyFile, 8);
        // int builder = 0;
        for (btllib::SeqReader::Record record; (record = assemblyReader.read());) {
            // if(builder % 100 == 0) {
            //     std::cout << "reading... builder: " << builder << "\n";
            // }
            // builder++;

            assemblyBF->insert(record.seq);
        }
        std::string bfmsg = "\n\nWriting Bloom filter to disk (" + linksArgParser.bfFile + ") : " + std::to_string(time(0)) + "\n";
        // assemblyruninfo += bfmsg;
        std::cout << bfmsg;
        assemblyBF->write("bftest.out");
        // std::cout << "Done mybf, printing stats...\n";
        printBloomStats(*assemblyBF, std::cout);
        
    }
    return assemblyBF;
}
void readFastaFastq(
                const std::string file,
                const btllib::BloomFilter& bloom, 
                std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS>>& matePair,
                std::unordered_set<uint64_t>& mates, 
                const uint64_t distance,
                const uint64_t k,
                const uint64_t step) {
        btllib::SeqReader longReader(file, 8); // CHECK FOR FLAG MODES

        //uint64_t delta = distance - (2 * k);
        uint64_t delta = distance;

        int breakFlag = 0;

        for (btllib::SeqReader::Record record; (record = longReader.read());) {
            btllib::NtHash nthash(record.seq, k, bloom.get_hash_num());
            btllib::NtHash nthashLead(record.seq, k, bloom.get_hash_num(), delta);

            for (size_t i = 0; nthash.roll() && nthashLead.roll(); i+=step) {
                // roll for the number of steps
                readFastAFastq_debug_counter_1++;

                // for step ----
                for(int j = 1; j < step; j++) {
                    std::cout << "Rolling one more time\n";
                    if(!nthashLead.roll() || !nthash.roll()) {
                        breakFlag = 1;
                        std::cout << "setting break flag\n";
                    }
                }   
                if(breakFlag){break;}
                // for step ----

                if(bloom.contains(nthash.hashes()) && bloom.contains(nthashLead.hashes())) { // May need to change with forward reverse hashes
                    readFastAFastq_debug_counter_2++;
                    
                    if(matePair.find(nthash.get_forward_hash()) == matePair.end()) {
                        matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()] = BT_IS(false, distance);
                    } else {
                        matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()].setIS(distance);
                    }
                    mates.insert(nthashLead.get_forward_hash());
                }else{
                    readFastAFastq_debug_counter_3++;
                }
            }
        }
}
void kmerize(std::string *seq,
                uint64_t frag_dist,
                uint64_t k,
                std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > *matePair,
                std::unordered_set<uint64_t>& mates,
                uint64_t step,
                std::string head,
                uint64_t endPosition,
                uint64_t initPosition,
                uint64_t pairct,
                btllib::KmerBloomFilter *bloom,
                uint64_t delta,
                uint64_t readLength) {
//TODO: change from main method to here
}
// matepair is from long reads Forward & reverse duos
// trackAll is from contigs Forward and reverse for both
void inline kmerizeContig( std::string *seq, 
                    std::unordered_map<uint64_t, KmerInfo>& trackAll,
                    std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > *matePair,
                    std::unordered_set<uint64_t>& mates,
                    uint64_t k,
                    std::string head,
                    unsigned hashFcts,
                    uint64_t step,
                    uint64_t &tmpCounter) {
    btllib::NtHash ntHashContig(*seq, k, hashFcts); // hashFunc can be 1 after first step

    int breakFlag = 0;
    unsigned seq_length = seq->length();
    for (size_t i = 0; ntHashContig.roll(); i+=step) {
        // for rolling step
        for(int j = 1; j < step; j++) {
            if(!ntHashContig.roll()) {
                breakFlag = 1;
            }
        }
        if(breakFlag) {break;}
        // for rolling step
        readContigs_debug_counter_1++;
        
        // Forward part
	    if(matePair->find(ntHashContig.get_forward_hash()) != matePair->end() || 
            mates.find(ntHashContig.get_forward_hash()) != mates.end()) {
            readContigs_debug_counter_2++;

            if(trackAll.find(ntHashContig.get_forward_hash()) == trackAll.end()) {
                trackAll[ntHashContig.get_forward_hash()] = KmerInfo(head, i, i + k, 0);
            } else {
                trackAll[ntHashContig.get_forward_hash()].incrementMultiple();
            }
        }

        // Reverse part
        if(matePair->find(ntHashContig.get_reverse_hash()) != matePair->end() || 
            mates.find(ntHashContig.get_reverse_hash()) != mates.end()) {
            readContigs_debug_counter_3++;

            if(trackAll.find(ntHashContig.get_reverse_hash()) == trackAll.end()) {
                trackAll[ntHashContig.get_reverse_hash()] = KmerInfo(head, i, i + k, 1);
            } else {
                trackAll[ntHashContig.get_reverse_hash()].incrementMultiple();
            }
        }
    }
}

void readContigs(
        std::string assemblyFile,
        std::unordered_map<uint64_t, KmerInfo>& trackAll,
        std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > matePair,
        std::unordered_set<uint64_t>& mates,
        std::unordered_map<std::string, uint64_t>& tigLength,
        uint64_t k,
        uint64_t minSize,
        unsigned hashFcts, 
        uint64_t step) {
    // std::cout << "hashFct in readContig: " << hashFcts << "\n";
    uint64_t cttig = 0;
    btllib::SeqReader contigReader(assemblyFile, 8);// CHANGE TO FLAGS LATER
    uint64_t tmpCounter = 0;
    for (btllib::SeqReader::Record record; (record = contigReader.read());) {
        cttig++;
        tigLength.insert({std::to_string(cttig), record.seq.length()});
        //std::cout << "\r" << cttig;
        if(record.seq.length() >= minSize) {
            kmerizeContig(&record.seq, trackAll, &matePair, mates, k, std::to_string(cttig), hashFcts, step, tmpCounter);
        }
    }
}

void pairContigs(
    std::string longReadsFile,
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > matePair,
    std::unordered_map<uint64_t, KmerInfo>& trackAll,
    std::unordered_map<std::string, uint64_t> tigLength,
    std::string issues,
    std::string distribution,
    uint64_t totalPairs,
    std::string tigpair_checkpoint,
    std::string simplepair_checkpoint,
    bool verbose,
    float insertStdev) {

    uint64_t ct_illogical = 0, ct_ok_contig = 0, ct_ok_pairs = 0, ct_problem_pairs = 0, ct_iz_issues = 0, ct_single = 0, ct_multiple = 0, ct_both = 0, trackInsert = 0;
    std::unordered_map<uint64_t, uint64_t> ct_single_hash, ct_both_hash, ct_illogical_hash, ct_ok_contig_hash, ct_ok_pairs_hash, ct_problem_pairs_hash, ct_iz_issues_hash;
    // Mapping of tiga_head -> insertSize -> tigb_head -> links & gaps
    std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, Gaps_Links> > > pair;
    std::unordered_map<std::string, std::unordered_map<std::string, Gaps_Links> >simplePair;
    std::unordered_map<std::string, Gaps_Links> err;
    std::string order1;
    std::string order2;
    if(verbose) std::cout << "Pairing contigs...\n";
    //******************
    int CheckCounterBase = 0;
    int filter1 = 0;
    int filter2 = 0;
    int filter3 = 0;
    int filter4 = 0;
    int Check0Counter = 0;
    int Check1Counter = 0;
    int Check2Counter = 0;
    int Check3Counter = 0;
    int Check4Counter = 0;
    int Check5Counter = 0;
    int Check6Counter = 0;
    int Check7Counter = 0;
    int Check8Counter = 0;
    int Check9Counter = 0;
    int Check10Counter = 0;
    int Check11Counter = 0;
    int Check12Counter = 0;
    int Check13Counter = 0;
    int Check14Counter = 0;
    int Check15Counter = 0;
    int Check16Counter = 0;
    int Check17Counter = 0;
    int Check18Counter = 0;
    int Check19Counter = 0;
    int Check20Counter = 0;
    int Check21Counter = 0;
    int Check22Counter = 0;
    int Check23Counter = 0;
    int Check24Counter = 0;
    int Check25Counter = 0;
    int Check26Counter = 0;
    //*****************

    std::ofstream issuesFile;
    issuesFile.open (issues);
    uint64_t insert_size = 0;
    int min_allowed = 0;
    int low_iz = 0;
    int up_iz = 0;

    int distance;
    int isz;

    std::string tig_a, tig_b, ftig_a, ftig_b, rtig_a, rtig_b;
    uint64_t A_length, A_start, A_end, B_length, B_start, B_end;

    KmerInfo kmer1, kmer2;

    std::unordered_map<uint64_t, BT_IS>::iterator mateListItr;
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> >::iterator matePairItr;

    for(matePairItr = matePair.begin(); matePairItr != matePair.end(); matePairItr++) {
        for(mateListItr = matePairItr->second.begin(); mateListItr != matePairItr->second.end(); mateListItr++) {
            CheckCounterBase++;

            if( mateListItr->second.getBT() == false &&                 //matepair is not seen
                trackAll.find(matePairItr->first) != trackAll.end() &&  //first mate is tracked
                trackAll[matePairItr->first].getMultiple() == 1 &&      //first mate seen once
                trackAll.find(mateListItr->first) != trackAll.end() &&  //second mate is tracked
                trackAll[mateListItr->first].getMultiple() == 1) {      //second mate is seen once
                mateListItr->second.setBT(true);

                insert_size = matePair[matePairItr->first][mateListItr->first].getIS();
                min_allowed = -1 * (insertStdev * insert_size); // check int
                low_iz = insert_size + min_allowed; // check int
                up_iz = insert_size - min_allowed; // check int

                if(verbose) std::cout << "Pair read1Hash=" << matePairItr->first << " read2Hash=" << mateListItr->first << "\n";

                if(trackAll[matePairItr->first].getTig() != "" && trackAll[mateListItr->first].getTig() != "") { //double check if tig names not null
                    ct_both++;
                    if(ct_both_hash.find(insert_size) == ct_both_hash.end()) {
                        ct_both_hash[insert_size] = 1;
                    } else {
                        ct_both_hash[insert_size] = ct_both_hash[insert_size] + 1;
                    }

                    /*
                    tig_a = trackAll[matePairItr->first].getTig();
                    tig_b = trackAll[mateListItr->first].getTig();

                    ftig_a = "f" + tig_a;
                    ftig_b = "f" + tig_b;

                    rtig_a = "r" + tig_a;
                    rtig_b = "r" + tig_b;

                    A_length = tigLength[tig_a];
                    A_start = trackAll[matePairItr->first].getStart();
                    A_end = trackAll[matePairItr->first].getEnd();

                    B_length = tigLength[tig_b];
                    B_start = trackAll[mateListItr->first].getStart();
                    B_end = trackAll[mateListItr->first].getEnd();
                    */
                    kmer1 = trackAll[matePairItr->first];
                    kmer2 = trackAll[mateListItr->first]; 

                    if(kmer1.getTig() != kmer2.getTig()) { // paired reads located on <> contigs
                        // MURATHAN DEBUG 11.5.21
                        pairContigs_debug_counter_8++;
                        if(!kmer1.getOrient()){             // if kmer1 is forward
                            pairContigs_debug_counter_9++;
                            if(!kmer2.getOrient()){         // if kmer2 is forward
                                pairContigs_debug_counter_10++;
                                distance = getDistance(insert_size, tigLength[kmer1.getTig()], kmer1.getStart(), kmer2.getStart());
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.getTig(), kmer2.getTig(), 0);
                                }
                            }else{                          // if kmer2 is reverse
                                pairContigs_debug_counter_11++;
                                distance = getDistance(insert_size, tigLength[kmer1.getTig()], kmer1.getStart(), tigLength[kmer2.getTig()] - kmer2.getEnd());
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.getTig(), kmer2.getTig(), 1);
                                }
                            }
                        }else{                              // if kmer1 is reverse
                            pairContigs_debug_counter_12++;
                            if(!kmer2.getOrient()){         // if kmer2 is forward
                                pairContigs_debug_counter_13++;
                                distance = getDistance(insert_size, tigLength[kmer2.getTig()], tigLength[kmer1.getTig()] - kmer1.getEnd(), kmer2.getStart());
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.getTig(), kmer2.getTig(), 2);
                                }
                            }else{                          // if kmer2 is reverse
                                pairContigs_debug_counter_14++;
                                distance = getDistance(insert_size, tigLength[kmer2.getTig()], kmer2.getEnd(), kmer1.getEnd());
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.getTig(), kmer2.getTig(), 3);
                                }
                            }
                        }
                        std::cout << "tig1: " << kmer1.getTig()  << " tig2: " << kmer2.getTig() << " distance: " << distance << " isz: " << isz << std::endl;
                        std::cout << "insert_size: " << insert_size << "distance: " << distance << " tig_a: " << kmer1.getTig() << " A_length: " << tigLength[kmer1.getTig()] << " A_start: " << kmer1.getStart() << " A_end: " << kmer1.getEnd() << " tig_b: " << kmer2.getTig()  << " B_start: " << kmer2.getStart() << " B_end: " << kmer2.getEnd() << std::endl; 

                        //Check2Counter++;
                        // std::cout << "Checkpoint 4 (if tigs are different)\n";
                        //Determine most likely possibility
                        // Checking if forward
                        // MURATHAN FIX 11.5.21
                        /*
                        if(A_start < A_end){
                        //if(trackFor.find(matePairItr->first) != trackFor.end()) {//trackAll[matePairItr->first].getStart() < trackAll[matePairItr->first].getEnd()) {
                        
                            Check3Counter++;
                            // std::cout << "Checkpoint 5 (A.start < A.end)\n";
                            // std::cout << "End: " << std::to_string(trackAll[mateListItr->first].getEnd()) << "Start: " << std::to_string(trackAll[mateListItr->first].getStart()) << "\n";
                            // Checking if reverse
                            // MURATHAN FIX 11.5.21
                            if(B_start > B_end){
                            //if(trackRev.find(matePairItr->first) != trackRev.end()) {//trackAll[mateListItr->first].getEnd() < trackAll[mateListItr->first].getStart()) { // -> <- :::  A-> <-B  /  rB -> <- rA
                                Check4Counter++;
                                // std::cout << "Checkpoint 6 (B.end < B.start)\n";
                                //std::cout << "insert_size: " << insert_size << " A_length: " << A_length << " B_start: " << B_start << std::endl; 
                                int distance = getDistance(insert_size, A_length, A_start, B_start);
                                //std::cout << "4 distance: " << distance << " tig_a: " << tig_a << " A_length: " << A_length << " A_start: " << A_start << " A_end: " << A_end << " tig_b: " << tig_b << " B_start: " << B_start << " B_end: " << B_end << std::endl; 
                                //std::cout << "distance: " << distance << " min allowed: " << min_allowed << std::endl;
                                if(verbose) std::cout << "A-> <-B  WITH " << tig_a << "-> <- " << tig_b << " GAP " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " " << std::to_string(B_start-B_end) << " Alen, Astart,Bstart\n";
                                if(distance > min_allowed) {
                                    Check5Counter++;
                                    // std::cout << "Checkpoint 7 distance > min allowed\n";
                                    int isz = distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 10000; // distance categories
                                    if(tig_a == "22" && tig_b == "11"){
                                        std::cout << "4 insert_size: " << insert_size << "distance: " << distance << " tig_a: " << tig_a << " A_length: " << A_length << " A_start: " << A_start << " A_end: " << A_end << " tig_b: " << tig_b << " B_start: " << B_start << " B_end: " << B_end << std::endl; 
                                    }
                                    if(pair.find(ftig_a) == pair.end() || pair[ftig_a].find(isz) == pair[ftig_a].end() || pair[ftig_a][isz].find(rtig_b) == pair[ftig_a][isz].end()) {
                                        // std::cout << "Checkpoint 7.1 adding to pair new GAPSLINKS\n";
                                        pair[ftig_a][isz][rtig_b] = Gaps_Links(distance,1);
                                        pairContigs_debug_counter_50++;
                                    } else {
                                        // std::cout << "Checkpoint 7.2 adding to pair existing gapslings\n";
                                        pair[ftig_a][isz][rtig_b].addToGap(distance);
                                        pair[ftig_a][isz][rtig_b].incrementLinks();
                                        pairContigs_debug_counter_51++;
                                    }
                                    if(pair.find(rtig_b) == pair.end() || pair[rtig_b].find(isz) == pair[rtig_b].end() || pair[rtig_b][isz].find(rtig_a) == pair[rtig_b][isz].end()) {
                                        // std::cout << "Checkpoint 7.3 adding to pair new GAPSLINKSs\n";
                                        pair[rtig_b][isz][rtig_a] = Gaps_Links(distance,1);
                                        pairContigs_debug_counter_52++;
                                    } else {
                                        pair[rtig_b][isz][rtig_a].addToGap(distance);
                                        pair[rtig_b][isz][rtig_a].incrementLinks();
                                        pairContigs_debug_counter_53++;
                                    }
                                    if(tig_a < tig_b) {
                                        order1 = tig_a;
                                        order2 = tig_b;
                                    } else {
                                        order1 = tig_b;
                                        order2 = tig_a;
                                    }
                                    // Check if exists
                                    simplePair[order1][order2] = Gaps_Links("11");
                                    simplePair[order1][order2].incrementLinks();
                                    simplePair[order1][order2].addToGap(distance);
                                    
                                    ct_ok_pairs++;
                                    if(ct_ok_pairs_hash.find(insert_size) == ct_ok_pairs_hash.end()) {
                                        ct_ok_pairs_hash[insert_size] = 1;
                                    } else {
                                        ct_ok_pairs_hash[insert_size] = ct_ok_pairs_hash[insert_size] + 1;
                                    }
                                } else {
                                    Check6Counter++;
                                    std::string err_pair = ftig_a + "-" + ftig_b;
                                    if(err.find(err_pair) == err.end()) {
                                        err[err_pair] = Gaps_Links(distance, 1);
                                    } else {
                                        err[err_pair].addToGap(distance);
                                        err[err_pair].incrementLinks();
                                    }
                                    ct_problem_pairs++;
                                    if(ct_problem_pairs_hash.find(insert_size) == ct_problem_pairs_hash.end()) {
                                        ct_problem_pairs_hash[insert_size] = 1;
                                    } else {
                                        ct_problem_pairs_hash[insert_size] = ct_problem_pairs_hash[insert_size] + 1;
                                    }
                                    // THIS IS NOT A DEBUGGING OUTPUT
                                    issuesFile << "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#" << tig_a << " -> " << std::to_string(distance) << " <- tig#"<< tig_b << ", A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") CALCULATED DISTANCE APART: " << distance << " < " << min_allowed << "\n";
                                }
                            } else { // -> -> ::: A-> <-rB  / B-> <-rA 
                                Check7Counter++;
                                uint64_t rB_start = B_length - B_start;
                                int distance = getDistance(insert_size, A_length, A_start, rB_start);
                                //std::cout << "7 distance: " << distance << " A_length: " << A_length << " A_start: " << A_start << " A_end: " << A_end << " B_start: " << B_start << " B_end: " << B_end << std::endl; 
                                //std::cout << "7 distance: " << distance << " tig_a: " << tig_a << " A_length: " << A_length << " A_start: " << A_start << " A_end: " << A_end << " tig_b: " << tig_b << " B_start: " << B_start << " B_end: " << B_end << std::endl; 

                                if(verbose) std::cout << "A-> <-rB  WITH " << tig_a << "-> <- " << tig_b << " GAP " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " " << std::to_string(B_start-B_end) << " Alen, Astart,rBstart\n";
                                //std::cout << "7 A-> <-rB  WITH " << tig_a << "-> <- " << tig_b << " GAP(dist) " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " rbstart" << std::to_string(rB_start) << " Alen, Astart,rBstart\n";
                                if(distance >= min_allowed) {
                                    Check8Counter++;
                                    // std::cout << "Checkpoint 10.1\n";
                                    int isz = distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 10000; // distance categories
                                    if(pair.find(ftig_a) == pair.end() || pair[ftig_a].find(isz) == pair[ftig_a].end() || pair[ftig_a][isz].find(rtig_b) == pair[ftig_a][isz].end()) {
                                        // std::cout << "Checkpoint 10.2\n";
                                        pair[ftig_a][isz][rtig_b] = Gaps_Links(distance,1);
                                    } else {
                                        // std::cout << "Checkpoint 10.3\n";
                                        pair[ftig_a][isz][rtig_b].addToGap(distance);
                                        pair[ftig_a][isz][rtig_b].incrementLinks();
                                    }
                                    if(pair.find(ftig_b) == pair.end() || pair[ftig_b].find(isz) == pair[ftig_b].end() || pair[ftig_b][isz].find(rtig_a) == pair[ftig_b][isz].end()) {
                                        // std::cout << "Checkpoint 10.4\n";
                                        pair[ftig_b][isz][rtig_a] = Gaps_Links(distance,1);
                                    } else {
                                        // std::cout << "Checkpoint 10.5\n";
                                        pair[ftig_b][isz][rtig_a].addToGap(distance);
                                        pair[ftig_b][isz][rtig_a].incrementLinks();
                                    }
                                    if(tig_a < tig_b) {
                                        order1 = tig_a;
                                        order2 = tig_b;
                                    } else {
                                        order1 = tig_b;
                                        order2 = tig_a;
                                    }
                                    simplePair[order1][order2] = Gaps_Links("10");
                                    simplePair[order1][order2].incrementLinks();
                                    simplePair[order1][order2].addToGap(distance);
                                    
                                    ct_ok_pairs++;
                                    if(ct_ok_pairs_hash.find(insert_size) == ct_ok_pairs_hash.end()) {
                                        ct_ok_pairs_hash[insert_size] = 1;
                                    } else {
                                        ct_ok_pairs_hash[insert_size] = ct_ok_pairs_hash[insert_size] + 1;
                                    }
                                } else {
                                    Check9Counter++;
                                    std::string err_pair = ftig_a + "-" + rtig_b;
                                    if(err.find(err_pair) == err.end()) {
                                        err[err_pair] = Gaps_Links(distance, 1);
                                    } else {
                                        err[err_pair].addToGap(distance);
                                        err[err_pair].incrementLinks();
                                    }
                                    ct_problem_pairs++;
                                    if(ct_problem_pairs_hash.find(insert_size) == ct_problem_pairs_hash.end()) {
                                        ct_problem_pairs_hash[insert_size] = 1;
                                    } else {
                                        ct_problem_pairs_hash[insert_size] = ct_problem_pairs_hash[insert_size] + 1;
                                    }
                                    issuesFile << "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#" << tig_a << " -> " << std::to_string(distance) << " <- tig#r."<< tig_b << ", A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") CALCULATED DISTANCE APART: " << distance << " < " << min_allowed << "\n";
                                }
                            }
                        } else {
                            Check10Counter++;
                            // if ({read_b}{'end'} > {$read_b}{'start'} (forward)
                            // MURATHAN FIX 11.5.21
                            if(B_end > B_start){
                            //if(trackFor.find(mateListItr->first) != trackFor.end()) {//trackAll[mateListItr->first].getEnd() > trackAll[mateListItr->first].getStart()) {
                                Check11Counter++;
                                int distance = getDistance(insert_size, B_length, B_start, A_start);
                                //std::cout << "11 distance: " << distance << " tig_a: " << tig_a << " A_length: " << A_length << " A_start: " << A_start << " A_end: " << A_end << " tig_b: " << tig_b << " B_start: " << B_start << " B_end: " << B_end << std::endl; 
                                if(tig_a == "11" && tig_b == "22"){
                                    std::cout << "11 insert_size: " << insert_size << "distance: " << distance << " tig_a: " << tig_a << " A_length: " << A_length << " A_start: " << A_start << " A_end: " << A_end << " tig_b: " << tig_b << " B_start: " << B_start << " B_end: " << B_end << std::endl; 
                                }
                                //std::cout << "11 distance: " << distance << " A_length: " << A_length << " A_start: " << A_start << " A_end: " << A_end << " B_start: " << B_start << " B_end: " << B_end << std::endl; 
                                if(verbose) std::cout << "B-> <-A  WITH " << tig_b << "-> <- " << tig_a << " GAP " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " " << std::to_string(B_start-B_end) << " Blen, Bstart,Astart\n";
                                if(distance >= min_allowed) {
                                    Check12Counter++;
                                    int isz = distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 10000; // distance categories
                                    if(pair.find(ftig_b) == pair.end() || pair[ftig_b].find(isz) == pair[ftig_b].end() || pair[ftig_b][isz].find(ftig_a) == pair[ftig_b][isz].end()) {
                                        // std::cout << "Checkpoint 11.1\n";
                                        pair[ftig_b][isz][ftig_a] = Gaps_Links(distance,1);
                                    } else {
                                        // std::cout << "Checkpoint 11.2\n";
                                        pair[ftig_b][isz][ftig_a].addToGap(distance);
                                        pair[ftig_b][isz][ftig_a].incrementLinks();
                                    }
                                    if(pair.find(rtig_a) == pair.end() || pair[rtig_a].find(isz) == pair[rtig_a].end() || pair[rtig_a][isz].find(rtig_b) == pair[rtig_a][isz].end()) {
                                        // std::cout << "Checkpoint 11.3\n";
                                        pair[rtig_a][isz][rtig_b] = Gaps_Links(distance,1);
                                    } else {
                                        // std::cout << "Checkpoint 11.4\n";
                                        pair[rtig_a][isz][rtig_b].addToGap(distance);
                                        pair[rtig_a][isz][rtig_b].incrementLinks();
                                    }
                                    if(tig_a < tig_b) {
                                        order1 = tig_a;
                                        order2 = tig_b;
                                    } else {
                                        order1 = tig_b;
                                        order2 = tig_a;
                                    }
                                    simplePair[order1][order2] = Gaps_Links("11");
                                    simplePair[order1][order2].incrementLinks();
                                    simplePair[order1][order2].addToGap(distance);
                                    
                                    ct_ok_pairs++;
                                    if(ct_ok_pairs_hash.find(insert_size) == ct_ok_pairs_hash.end()) {
                                        ct_ok_pairs_hash[insert_size] = 1;
                                    } else {
                                        ct_ok_pairs_hash[insert_size] = ct_ok_pairs_hash[insert_size] + 1;
                                    }
                                } else {
                                    Check13Counter++;
                                    std::string err_pair = ftig_b + "-" + ftig_a;
                                    if(err.find(err_pair) == err.end()) {
                                        err[err_pair] = Gaps_Links(distance, 1);
                                    } else {
                                        err[err_pair].addToGap(distance);
                                        err[err_pair].incrementLinks();
                                    }
                                    ct_problem_pairs++;
                                    if(ct_problem_pairs_hash.find(insert_size) == ct_problem_pairs_hash.end()) {
                                        ct_problem_pairs_hash[insert_size] = 1;
                                    } else {
                                        ct_problem_pairs_hash[insert_size] = ct_problem_pairs_hash[insert_size] + 1;
                                    }
                                    issuesFile << "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#" << tig_b << " -> " << std::to_string(distance) << " <- tig#"<< tig_a << ", B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") CALCULATED DISTANCE APART: " << distance << " < " << min_allowed << "\n";
                                }
                            }  else { // <- <-  :::  rB-> <-A / rA-> <-B
                                Check14Counter++;
                                uint64_t rB_start = B_length - B_start;
                                int distance = getDistance(insert_size, B_length, rB_start, A_start);
                                //std::cout << "14 distance: " << distance << " tig_a: " << tig_a << " A_length: " << A_length << " A_start: " << A_start << " A_end: " << A_end << " tig_b: " << tig_b << " B_start: " << B_start << " B_end: " << B_end << std::endl; 

                                //std::cout << "14 distance: " << distance << " A_length: " << A_length << " A_start: " << A_start << " A_end: " << A_end << " B_start: " << B_start << " B_end: " << B_end << std::endl; 
                                if(verbose) std::cout << "rB-> <-A  WITH r." << tig_b << "-> <- " << tig_a << " GAP " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " " << std::to_string(B_start-B_end) << " Blen, rBstart,Astart\n";
                                //std::cout << "14 rB-> <-A  WITH r." << tig_b << "-> <- " << tig_a << " GAP " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " rb start: " << std::to_string(rB_start) << " Blen, rBstart,Astart\n";
                                if(distance >= min_allowed) {
                                    Check15Counter++;
                                    int isz = distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 10000; // distance categories
                                    if(pair.find(rtig_b) == pair.end() || pair[rtig_b].find(isz) == pair[rtig_b].end() || pair[rtig_b][isz].find(rtig_b) == pair[ftig_a][isz].end()) {
                                        // std::cout << "Checkpoint 12.1\n";
                                        pair[rtig_b][isz][ftig_a] = Gaps_Links(distance,1);
                                    } else {
                                        // std::cout << "Checkpoint 12.2\n";
                                        pair[rtig_b][isz][ftig_a].addToGap(distance);
                                        pair[rtig_b][isz][ftig_a].incrementLinks();
                                    }
                                    if(pair.find(rtig_a) == pair.end() || pair[rtig_a].find(isz) == pair[rtig_a].end() || pair[rtig_a][isz].find(ftig_b) == pair[rtig_a][isz].end()) {
                                        // std::cout << "Checkpoint 12.3\n";
                                        pair[rtig_a][isz][ftig_b] = Gaps_Links(distance,1);
                                    } else {
                                        // std::cout << "Checkpoint 12.4\n";
                                        pair[rtig_a][isz][ftig_b].addToGap(distance);
                                        pair[rtig_a][isz][ftig_b].incrementLinks();
                                    }
                                    if(tig_a < tig_b) {
                                        order1 = tig_a;
                                        order2 = tig_b;
                                    } else {
                                        order1 = tig_b;
                                        order2 = tig_a;
                                    }
                                    simplePair[order1][order2] = Gaps_Links("01");
                                    simplePair[order1][order2].incrementLinks();
                                    simplePair[order1][order2].addToGap(distance);
                                    
                                    ct_ok_pairs++;
                                    if(ct_ok_pairs_hash.find(insert_size) == ct_ok_pairs_hash.end()) {
                                        ct_ok_pairs_hash[insert_size] = 1;
                                    } else {
                                        ct_ok_pairs_hash[insert_size] = ct_ok_pairs_hash[insert_size] + 1;
                                    }
                                } else {
                                    Check16Counter++;
                                    std::string err_pair = rtig_b + "-" + ftig_a;
                                    if(err.find(err_pair) == err.end()) {
                                        err[err_pair] = Gaps_Links(distance, 1);
                                    } else {
                                        err[err_pair].addToGap(distance);
                                        err[err_pair].incrementLinks();
                                    }
                                    ct_problem_pairs++;
                                    if(ct_problem_pairs_hash.find(insert_size) == ct_problem_pairs_hash.end()) {
                                        ct_problem_pairs_hash[insert_size] = 1;
                                    } else {
                                        ct_problem_pairs_hash[insert_size] = ct_problem_pairs_hash[insert_size] + 1;
                                    }
                                    issuesFile << "Pairs unsatisfied in distance within a contig pair.  rB-> <-A  WITH tig#r." << tig_b << " -> " << std::to_string(distance) << " <- tig#"<< tig_a << ", B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") CALCULATED DISTANCE APART: " << distance << " < " << min_allowed << "\n";
                                }
                            }
                        }
                        */
                    } else { // Clone, paired reads located on the same contig -- could be used to investigate misassemblies
                        Check17Counter++;
                        if (verbose) std::cout << "Pair (" << matePairItr->first << " and " << mateListItr->first << ") located on same contig " << tig_a << " (" << A_length << " nt)\n";
                        uint64_t pet_size = 0;
                        if(A_start > B_start && (B_start < B_end) && (A_start > A_end)) {   // B --> <-- A
                            Check18Counter++;
                            pet_size = A_start - B_start;
                            trackInsert += pet_size;
                            if(pet_size >= low_iz && pet_size <= up_iz) {
                                Check19Counter++;
                                ct_ok_contig++;
                                if(ct_ok_contig_hash.find(insert_size) == ct_ok_contig_hash.end()) {
                                    ct_ok_contig_hash[insert_size] = 1;
                                } else {
                                    ct_ok_contig_hash[insert_size] = ct_ok_contig_hash[insert_size] + 1;
                                }
                            } else {
                                Check20Counter++;
                                issuesFile << "Pairs unsatisfied in distance within a contig.  Pair (" << matePairItr->first << " - " << mateListItr->first << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << " CALCULATED DISTANCE APART: " << pet_size << "\n";
                                ct_iz_issues++;
                                if(ct_iz_issues_hash.find(insert_size) == ct_iz_issues_hash.end()) {
                                    ct_iz_issues_hash[insert_size] = 1;
                                } else {
                                    ct_iz_issues_hash[insert_size] = ct_iz_issues_hash[insert_size] + 1;
                                }
                            }
                        } else if(B_start > A_start && (B_start > B_end) && (A_start < A_end)) { // A --> <-- B
                            Check21Counter++;
                            pet_size = B_start - A_start;
                            trackInsert += pet_size;
                            if(pet_size >= low_iz && pet_size <= up_iz) {
                                Check22Counter++;
                                ct_ok_contig++;
                                if(ct_ok_contig_hash.find(insert_size) == ct_ok_contig_hash.end()) {
                                    ct_ok_contig_hash[insert_size] = 1;
                                } else {
                                    ct_ok_contig_hash[insert_size] = ct_ok_contig_hash[insert_size] + 1;
                                }
                            } else {
                                Check23Counter++;
                                issuesFile << "Pairs unsatisfied in distance within a contig.  Pair (" << matePairItr->first << " - " << mateListItr->first << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << "\n";
                                ct_iz_issues++;
                                if(ct_iz_issues_hash.find(insert_size) == ct_iz_issues_hash.end()) {
                                    ct_iz_issues_hash[insert_size] = 1;
                                } else {
                                    ct_iz_issues_hash[insert_size] = ct_iz_issues_hash[insert_size] + 1;
                                }
                            }
                        } else {
                            Check24Counter++;
                            ct_illogical++;
                            if(ct_illogical_hash.find(insert_size) == ct_illogical_hash.end()) {
                                ct_illogical_hash[insert_size] = 1;
                            } else {
                                ct_illogical_hash[insert_size] = ct_illogical_hash[insert_size] + 1;
                            }
                            // FOLLOWING IS NOT A DEBUGGING PRINT
                            issuesFile << "Pairs unsatisfied in pairing logic within a contig.  Pair (" << matePairItr->first << " - " << mateListItr->first << ") on contig" << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << "\n";
                        }
                    }
                } else { // both pairs assembled
                    Check25Counter++;
                    ct_single++;
                    if(ct_single_hash.find(insert_size) == ct_single_hash.end()) {
                        ct_single_hash[insert_size] = 1;
                    } else {
                        ct_single_hash[insert_size] = ct_single_hash[insert_size] + 1;
                    }
                }
            } else { // if unseen
                Check26Counter++;
                // std::cout << "UNSEEN\n";
                if(matePair[matePairItr->first][mateListItr->first].getBT() == false) {
                    std::cout << "UNSEEN getBT() increment ct\n";
                    ct_multiple++;
                }

            }
        } // pairing read b
    } // pairing read a

    // Summary of the contig pair issues
    //std::cout << "------------- Putative issues with contig pairing - Summary  ----------------\n";
    //std::cout << "err map size: " << err.size() << "\n"; 
    sortErr(err);
    std::unordered_map<std::string, Gaps_Links>::iterator errItr;
    for(errItr = err.begin(); errItr != err.end(); errItr++) {
        double mean_iz = 0;
        if(errItr->second.getLinks()) {
            mean_iz = errItr->second.getGaps() / errItr->second.getLinks();
        }
        // std::cout << "Pair " << errItr->first << " has " << errItr->second.getLinks() << " Links and mean distance of = " << mean_iz << "\n";
    }

    uint64_t satisfied = ct_ok_pairs + ct_ok_contig;
    uint64_t unsatisfied = ct_problem_pairs + ct_iz_issues + ct_illogical;
    uint64_t ct_both_reads = ct_both * 2;

    
    std::cout << "THESE ARE THE FILTERINGS:\n"<< "filter 1: "<< std::to_string(filter1) << "\n" << "filter 2: "<< std::to_string(filter2) << "\n" << "filter 3: "<< std::to_string(filter3) << "\n" << "filter 4: "<< std::to_string(filter4) << "\n";
    std::cout << "THESE ARE THE COUNTERS:\n" << std::to_string(CheckCounterBase) << "\n0 " << std::to_string(Check0Counter) << "\n1 " <<  std::to_string(Check1Counter) << "\n2 " <<  std::to_string(Check2Counter) << "\n3 " <<  std::to_string(Check3Counter) << "\n4 " <<  std::to_string(Check4Counter) << "\n5 " <<  std::to_string(Check5Counter) << "\n6 " <<  std::to_string(Check6Counter) << "\n7 " <<  std::to_string(Check7Counter) << "\n8 " <<  std::to_string(Check8Counter) << "\n9 " <<  std::to_string(Check9Counter) << "\n10 " << std::to_string(Check10Counter) << "\n11 " << std::to_string(Check11Counter) << "\n12 " << std::to_string(Check12Counter) << "\n13 " << std::to_string(Check13Counter) << "\n14 " << std::to_string(Check14Counter) << "\n15 " << std::to_string(Check15Counter) << "\n16 " << std::to_string(Check16Counter) << "\n17 " << std::to_string(Check17Counter) << "\n18 " << std::to_string(Check18Counter) << "\n19" << std::to_string(Check19Counter) << "\n20 " << std::to_string(Check20Counter) << "\n21 " << std::to_string(Check21Counter) << "\n22 " << std::to_string(Check22Counter) << "\n23 " << std::to_string(Check23Counter) << "\n24 " << std::to_string(Check24Counter) << "\n25 " << std::to_string(Check25Counter) << "\n26 " << std::to_string(Check26Counter) << "\n";
    std::cout << "\n===========PAIRED K-MER STATS===========\n";
    std::cout << "Total number of pairs extracted from -s " << longReadsFile << " " << totalPairs << "\n";
    std::cout << "At least one sequence/pair missing from contigs: " << ct_single << "\n";
    std::cout << "Ambiguous kmer pairs (both kmers are ambiguous): " << ct_multiple << "\n";
    std::cout << "Assembled pairs: " << ct_both << " (" << ct_both_reads << " sequences)\n";
    std::cout << "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: " << ct_ok_contig << "\n";
    std::cout << "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds): " << ct_iz_issues << "\n";
    std::cout << "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): " << ct_illogical << "\n";
    std::cout << "\t---\n";
    std::cout << "\tSatisfied in distance/logic within a given contig pair (pre-scaffold): " << ct_ok_pairs << "\n";
    std::cout << "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): " << ct_problem_pairs << "\n";
    std::cout << "\t---\n";
    std::cout << "Total satisfied: " << satisfied << "\tunsatisfied: " << unsatisfied << "\n\nBreakdown by distances (-d):\n";
    
    sortInsertSize(ct_both_hash);
    std::unordered_map<uint64_t, uint64_t>::iterator itrIS;
    std::cout << "ct_both_hash map size: " << err.size() << "\n"; 
    for(itrIS = ct_both_hash.begin(); itrIS != ct_both_hash.end(); itrIS++) {
        std::cout <<  "--------k-mers separated by "<< itrIS->first << " bp (outer distance)--------\n";
        int64_t maopt = -1 * (insertStdev * itrIS->first);
        int64_t low_izopt = itrIS->first + maopt;
        int64_t up_izopt = itrIS->first - maopt;
        /*
        std::cout <<  "MIN: " << low_izopt << " MAX: " << up_izopt << "  as defined by  " << itrIS->first << "  *  " << insertStdev << " \n";
        std::cout <<  "At least one sequence/pair missing:  " << ct_single_hash[itrIS->first] << " \n";
        std::cout <<  "Assembled pairs:  " << itrIS->second << " \n";
        std::cout <<  "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target:  " << ct_ok_contig_hash[itrIS->first] << " \n";
        std::cout <<  "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds):  " << ct_iz_issues_hash[itrIS->first] << " \n";
        std::cout <<  "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->):  " << ct_illogical_hash[itrIS->first] << " \n";
        std::cout <<  "\t---\n";
        std::cout <<  "\tSatisfied in distance/logic within a given contig pair (pre-scaffold):  " << ct_ok_pairs_hash[itrIS->first] << " \n";
        std::cout <<  "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds):  " << ct_problem_pairs_hash[itrIS->first] << " \n";
        */
    }
    std::cout << "============================================\n";
    std::ofstream distFile;
    distFile.open ("distfile.txt");
    // if (distFile.is_open()) {} else {}
    // foreach my $is (sort {$a<=>$b} keys %$track_insert){
    //     print CSV "$is,$track_insert->{$is}\n";
    // }
    distFile.close();

    // TIGPAIR CHECKPOINT
    std::ofstream tigpairCheckpointFile;
    tigpairCheckpointFile.open (tigpair_checkpoint);
    std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, Gaps_Links> > >::iterator pairItr;
    std::cout << "size of TIGPAIR: " << pair.size() << "\n";
    for(pairItr = pair.begin(); pairItr != pair.end(); pairItr++) {
        std::unordered_map<int64_t, std::unordered_map<std::string, Gaps_Links> >::iterator insertSizes;
        for(insertSizes = pairItr->second.begin(); insertSizes != pairItr->second.end(); insertSizes++) {
            std::unordered_map<std::string, Gaps_Links>::iterator secPairItr;
            for(secPairItr = insertSizes->second.begin(); secPairItr != insertSizes->second.end(); secPairItr++) {
                tigpairCheckpointFile << insertSizes->first /*distance*/ << "\t" << pairItr->first << "\t" << secPairItr->first  << "\t" << secPairItr->second.getLinks()  << "\t" << secPairItr->second.getGaps() << "\n";
            }
        }
    }
    // if (distFile.is_open()) {} else {}
    // foreach my $is (sort {$a<=>$b} keys %$track_insert){
    //     print CSV "$is,$track_insert->{$is}\n";
    // }
    tigpairCheckpointFile.close();

}

void addToPairMap(
    int& isz,
    std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, Gaps_Links>>>& pair,
    int& distance,
    std::string kmer1_name,
    std::string kmer2_name,
    unsigned orient_enum
    ){
        std::tuple<std::string, std::string> first_pair;
        std::tuple<std::string, std::string> second_pair;

        std::string ftig_a = "f" + kmer1_name;
        std::string ftig_b = "f" + kmer2_name;

        std::string rtig_a = "r" + kmer1_name;
        std::string rtig_b = "r" + kmer2_name;

        switch (orient_enum)
        {
        case 0:                 // A B or rB rA
            first_pair = std::make_tuple(ftig_a,ftig_b);
            second_pair = std::make_tuple(rtig_b,rtig_a);
            break;
        case 1:                 // A B or A rB
            first_pair = std::make_tuple(ftig_a,rtig_b);
            second_pair = std::make_tuple(ftig_b,rtig_a);
            break;
        case 2:                 // rA B or rB A
            first_pair = std::make_tuple(rtig_a,ftig_b);
            second_pair = std::make_tuple(rtig_b,ftig_a);
            break;
        case 3:                 // rA rB or B A
            first_pair = std::make_tuple(rtig_a,rtig_b);
            second_pair = std::make_tuple(ftig_b,ftig_a);
            break;
        default:
            break;
        }
        if(pair.find(std::get<0>(first_pair)) == pair.end() 
            || pair[std::get<0>(first_pair)].find(isz) == pair[std::get<0>(first_pair)].end() 
            || pair[std::get<0>(first_pair)][isz].find(std::get<1>(first_pair)) == pair[std::get<0>(first_pair)][isz].end()) {
            // std::cout << "Checkpoint 7.1 adding to pair new GAPSLINKS\n";
            pair[std::get<0>(first_pair)][isz][std::get<1>(first_pair)] = Gaps_Links(distance,1);
        } else {
            // std::cout << "Checkpoint 7.2 adding to pair existing gapslings\n";
            pair[std::get<0>(first_pair)][isz][std::get<1>(first_pair)].addToGap(distance);
            pair[std::get<0>(first_pair)][isz][std::get<1>(first_pair)].incrementLinks();
        }
        if(pair.find(std::get<0>(second_pair)) == pair.end() 
            || pair[std::get<0>(second_pair)].find(isz) == pair[std::get<0>(second_pair)].end() 
            || pair[std::get<0>(second_pair)][isz].find(std::get<1>(second_pair)) == pair[std::get<0>(second_pair)][isz].end()) {
            // std::cout << "Checkpoint 7.3 adding to pair new GAPSLINKSs\n";
            pair[std::get<0>(second_pair)][isz][std::get<1>(second_pair)] = Gaps_Links(distance,1);
        } else {
            pair[std::get<0>(second_pair)][isz][std::get<1>(second_pair)].addToGap(distance);
            pair[std::get<0>(second_pair)][isz][std::get<1>(second_pair)].incrementLinks();
        }
}

int getDistance(uint64_t insert_size, uint64_t length_i, uint64_t start_i, uint64_t start_j) {

   int insert_span = (length_i - start_i) + start_j;
   int gap_or_overlap = insert_size - insert_span;

   return gap_or_overlap;
}

void printBloomStats(btllib::KmerBloomFilter& bloom, std::ostream& os)
{
	os << "Bloom filter stats:"
	   << "\n\t#counters               = " << bloom.get_occupancy()//getFilterSize()
	   << "\n\t#size (B)               = " << bloom.get_bytes()
	   << "\n\tpopcount                = " << bloom.get_pop_cnt()
	   << "\n\tFPR                     = " << 100.f * bloom.get_fpr() << "%"
	   << "\n";
}

uint64_t getFileSize(std::string filename)
{
    // This buffer is a stat struct that the information is placed concerning the file.
    struct stat stat_buf;
    // Stat method returns true if successfully completed
    int rc = stat(filename.c_str(), &stat_buf);
    // st_size holds the total size of the file in bytes
    return rc == 0 ? stat_buf.st_size : -1;
}

inline bool exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}








//**********SORTER*********
bool cmpErr(std::pair<std::string, Gaps_Links>& a, 
         std::pair<std::string, Gaps_Links>& b) 
{ 
    return a.second.getLinks() < b.second.getLinks(); 
}

bool cmpIS(std::pair<uint64_t, uint64_t>& a, 
         std::pair<uint64_t, uint64_t>& b) 
{ 
    return a.second < b.second; 
} 
  
// // Function to sort the map according 
// // to value in a (key-value) pairs 
void sortErr(std::unordered_map<std::string, Gaps_Links>& M) 
{ 
  
    // Declare vector of pairs 
    std::vector<std::pair<std::string, Gaps_Links> > A; 
  
    // Copy key-value pair from Map 
    // to vector of pairs 
    // std::unordered_map<std::string, Gaps_Links>::iterator itr;
    // for(itr = M.begin(); itr != M.end(); itr++) {
    //     A.push_back(itr);
    // }
    for (auto& it : M) { 
        A.push_back(it); 
    } 
  
    // Sort using comparator function 
    sort(A.begin(), A.end(), cmpErr); 
  
    // Print the sorted value 
    for (auto& it : A) { 
  
        std::cout << it.first << ' '
             << it.second.getLinks(); 
    } 
} 

void sortInsertSize(std::unordered_map<uint64_t, uint64_t>& M) 
{ 
  
    // Declare vector of pairs 
    std::vector<std::pair<uint64_t, uint64_t> > A; 
  
    // Copy key-value pair from Map 
    // to vector of pairs 
    // std::unordered_map<std::string, Gaps_Links>::iterator itr;
    // for(itr = M.begin(); itr != M.end(); itr++) {
    //     A.push_back(itr);
    // }
    for (auto& it : M) { 
        A.push_back(it); 
    } 
  
    // Sort using comparator function 
    sort(A.begin(), A.end(), cmpIS); 
  
    // Print the sorted value 
    for (auto& it : A) { 
        std::cout << it.first << ' '
             << it.second; 
    } 
} 

bool does_file_exist(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}