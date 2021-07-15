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


#include "btllib/bloom_filter.hpp"
#include "btllib/seq_reader.hpp"

//Globals
#define BASE_TEN 10
std::string version = "2.0";

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

        std::stringstream ss(input);
        std::cout << input << std::endl;

        while( ss.good() )
        {
            std::cout << "here1" << std::endl;
            std::string substr;
            getline( ss, substr, ',' );
            distances.push_back(static_cast<unsigned int>(std::stoul(substr) ));
            std::cout << distances.back() << std::endl;
        }
/*         std::istringstream ss(input); // Turn the string into a stream. 
        std::string tok; 

        while (ss >> tok) 
        {
            distances.push_back(static_cast<unsigned int>(std::stoul(tok)));
            //std::cout << distances.back() << "\n";
        } */
        for(uint64_t d : distances){
            std::cout << "d: " << d << std::endl;
        }

        //std::sort(distances.begin(),distances.end());
        //std::cout << distances.back() << "\n";

        return distances;
    }

    public:
    std::string assemblyFile;
    std::string fofFile;
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
    // $base_name = $assemblyfile . ".scaff_s-" . $fofFile . "_d" . $distances . "_k" . $k . "_e" . $insert_stdev . "_l" . $min_links . "_a" . $max_link_ratio . "_z" . $min_size . "_t" . $step . "_o" . $offset . "_r-" . $bf_file . "_p" . $fpr . "_x" . $bfoff . "_m" . $readlength;



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
                    fofFile.assign(optarg);
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
                    << "  -s " << fofFile << "\n" 
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

//typedef
typedef std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > mate_pair;

// Helper Methods
void sortErr(std::unordered_map<std::string, Gaps_Links>& M);
void sortInsertSize(std::unordered_map<uint64_t, uint64_t>& M);
void printBloomStats(btllib::KmerBloomFilter& bloom, std::ostream& os);
bool does_file_exist(std::string fileName);
btllib::KmerBloomFilter *makeBF(uint64_t bfElements, InputParser linksArgParser);
uint64_t getFileSize(std::string filename);
void readFastaFastq_debug(
                const std::string file,
                const btllib::BloomFilter& bloom, 
                std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS>>& matePair,
                std::unordered_map<uint64_t, KmerInfo>& trackAll,
                std::unordered_set<uint64_t>& mates, 
                const uint64_t distance,
                const uint64_t k,
                const uint64_t step);
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
    //time_t now = time(0);
   
    // convert now to string form
    //char* dt = ctime(&now);
    // Parse command line arguments
    InputParser linksArgParser = InputParser(argc, argv);

    //Set Bloom Filter element number based on the size of the assembly file (1 byte = 1 character)
    std::string assemblyPath = linksArgParser.assemblyFile;
    int64_t bfElements = getFileSize(assemblyPath);

    // Checking validity of input assemble file
    if(bfElements == -1){
        std::cout << std::to_string(bfElements);
        std::cout << "Invalid file: " << linksArgParser.assemblyFile << " -- fatal\n";
        return -1;
    }

    // Naming output files
    if (linksArgParser.baseName == "") {
        linksArgParser.baseName =  linksArgParser.assemblyFile + ".scaff_s-" + 
                                    linksArgParser.fofFile + "_d" + 
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
                                "\n-s " + linksArgParser.fofFile + 
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
        // // k-merize long reads
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> > matePair;



    std::unordered_map<uint64_t, KmerInfo> trackAll;
    std::unordered_map<std::string, uint64_t> tigLength;

    // store second mates in a set
    std::unordered_set<uint64_t> mates;

        //test purpose
    matePair.reserve(50000000);
    mates.reserve(50000000);
    trackAll.reserve(50000000);
    

    // Stage 1 -- populate bloom filter with contigs
    btllib::KmerBloomFilter * myFilter;
    std::cout << "\n\n=>Before makeBF in c++ " + std::to_string(time(0)) + "\n";
    myFilter = makeBF(bfElements, linksArgParser);
    std::cout << "\n\n=>After makeBF in c++ " + std::to_string(time(0)) + "\n";
    unsigned hashFct = myFilter->get_hash_num();
    btllib::BloomFilter& filtering = myFilter->get_bloom_filter();


    // Stage 2 -- read long-reads and find mates
    std::ifstream infile(linksArgParser.fofFile);
    std::string readFile;
    while(infile >> readFile){
        for(uint64_t dist : linksArgParser.distances){
            std::cout << dist << std::endl;
            std::cout << readFile << std::endl;
            readFastaFastq(readFile,filtering,matePair,mates,dist,linksArgParser.k,linksArgParser.step);
        }
    }
    std::cout << "matepair size: " << matePair.size() << std::endl;
    std::cout << "mates size: " << mates.size() << std::endl;

    // Stage 3 -- read contigs to assign info to mates
    readContigs(linksArgParser.assemblyFile, trackAll, matePair, mates, tigLength, linksArgParser.k, linksArgParser.minSize, hashFct, linksArgParser.step);
   
    std::cout << "trackAll size: " << trackAll.size() << std::endl;
    // Stage 4 -- pair contigs based on mates
    uint64_t totalpairs = 0;
    pairContigs(
        linksArgParser.fofFile,
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

    return 0;
}

btllib::KmerBloomFilter *makeBF(uint64_t bfElements, InputParser linksArgParser) {
    btllib::KmerBloomFilter * assemblyBF;
    if(linksArgParser.bfFile != "") {
        std::cout << "A Bloom filter was supplied (" << linksArgParser.bfFile << ") and will be used instead of building a new one from -f " << linksArgParser.assemblyFile << "\n";
        if(!does_file_exist(linksArgParser.bfFile)) {
            std::cout << "\nInvalid file: " << linksArgParser.bfFile <<  " -- fatal\n";
            exit(1);
        } else {
            std::cout << "Checking Bloom filter file " << linksArgParser.bfFile <<"...ok\n";
        }
        // std::cout << "Loading bloom filter of size " << getFileSize(linksArgParser.bfFile) << " from " << linksArgParser.bfFile << "\n";
        assemblyBF = new btllib::KmerBloomFilter(linksArgParser.bfFile);
    } else {
        uint64_t m = ceil((-1 * (double)bfElements * log(linksArgParser.fpr)) / (log(2) * log(2)));
        //uint64_t rem = 64 - (m % 64);
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
        assemblyBF->save("bftest.out");
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
        std::cout << "here 1" << std::endl;
        //uint64_t delta = distance - (2 * k);
        uint64_t delta = distance - k;
        //uint64_t delta = distance;
        std::cout << "here 2" << std::endl;
        int breakFlag = 0;
        bool reverseExists = false;
        uint read_counter = 0;
        std::cout << "here 3" << std::endl;

        for (btllib::SeqReader::Record record; (record = longReader.read());) {
            if(read_counter % 1000 == 0){
                std::cout << "read counter: " << read_counter << std::endl;
            }
            //if(read_counter > 20000){
            //    break;
            //}
            ++read_counter;

            btllib::NtHash nthash(record.seq, bloom.get_hash_num(), k);
            btllib::NtHash nthashLead(record.seq, bloom.get_hash_num(), k, delta);

            for (size_t i = 0; nthash.roll() && nthashLead.roll(); i+=step) {
                // roll for the number of steps
                breakFlag = 0;
                reverseExists = false;

                // for step ----
                for(uint j = 1; j < step; j++) {
                    if(!nthashLead.roll() || !nthash.roll()) {
                        breakFlag = 1;
                    }
                }   
                if(breakFlag){break;}
                // for step ----

                // check if reverse pair exist
                mate_pair::iterator it = matePair.find(nthashLead.get_reverse_hash());
                if( it != matePair.end() ) {
                    std::unordered_map<uint64_t, BT_IS> &innerMap = it->second;
                    std::unordered_map<uint64_t, BT_IS>::iterator innerit = innerMap.find(nthash.get_reverse_hash());
                    if( innerit != innerMap.end() ){
                        innerit->second.setIS(distance);
                        reverseExists = true;
                    }
                }
               
                if(!reverseExists && bloom.contains(nthash.hashes()) && bloom.contains(nthashLead.hashes())) { // May need to change with forward reverse hashes
                    mate_pair::iterator it = matePair.find(nthash.get_forward_hash());
                    if( it != matePair.end() ) {
                        std::unordered_map<uint64_t, BT_IS> &innerMap = it->second;
                        std::unordered_map<uint64_t, BT_IS>::iterator innerit = innerMap.find(nthashLead.get_forward_hash());
                        if( innerit != innerMap.end() ){
                            matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()].setIS(distance);
                        }else{
                            matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()] = BT_IS(false, distance);
                        }
                    }else{
                        matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()] = BT_IS(false, distance);
                    }
                    
/*                     if(matePair.find(nthash.get_forward_hash()) == matePair.end()) {
                        matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()] = BT_IS(false, distance);
                        readFastAFastq_debug_counter_5++;
                    } else {
                        matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()].setIS(distance);
                        readFastAFastq_debug_counter_6++;
                    } */
                    mates.insert(nthashLead.get_forward_hash());
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
    btllib::NtHash ntHashContig(*seq, hashFcts, k); // hashFunc can be 1 after first step

    int breakFlag = 0;
    //unsigned seq_length = seq->length();
    for (size_t i = 0; ntHashContig.roll(); i+=step) {
        // for rolling step
        for(uint j = 1; j < step; j++) {
            if(!ntHashContig.roll()) {
                breakFlag = 1;
            }
        }
        if(breakFlag) {break;}
        // for rolling step

        i = ntHashContig.get_pos();
        
        // Forward part
	    if(matePair->find(ntHashContig.get_forward_hash()) != matePair->end() || 
            mates.find(ntHashContig.get_forward_hash()) != mates.end()) {

            if(trackAll.find(ntHashContig.get_forward_hash()) == trackAll.end()) {
                trackAll[ntHashContig.get_forward_hash()] = KmerInfo(head, i, i + k, 0);
            } else {
                trackAll[ntHashContig.get_forward_hash()].incrementMultiple();
            }
        }

        // Reverse part
        if(matePair->find(ntHashContig.get_reverse_hash()) != matePair->end() || 
            mates.find(ntHashContig.get_reverse_hash()) != mates.end()) {

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
        // for debug purposes
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

    std::ofstream issuesFile;
    issuesFile.open (issues);
    int64_t insert_size = 0;
    int min_allowed = 0;
    uint low_iz = 0;
    uint up_iz = 0;

    int distance;
    int isz;

    std::string tig_a, tig_b, ftig_a, ftig_b, rtig_a, rtig_b;
    uint64_t A_length = 0, A_start = 0, A_end = 0, B_length = 0, B_start = 0, B_end = 0;

    KmerInfo kmer1, kmer2;

    std::unordered_map<uint64_t, BT_IS>::iterator mateListItr;
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, BT_IS> >::iterator matePairItr;

  uint counter_1 = 0, counter_2 = 0, counter_3 = 0, counter_4 = 0, counter_5 = 0, counter_6 = 0, counter_7 = 0, counter_8 = 0;
    for(matePairItr = matePair.begin(); matePairItr != matePair.end(); matePairItr++) {
        for(mateListItr = matePairItr->second.begin(); mateListItr != matePairItr->second.end(); mateListItr++) {
            //std::cout << "*****" << std::endl;
            //std::cout << "trackAll[matePairItr->first].getMultiple(): " << trackAll[matePairItr->first].getMultiple() << std::endl;
            //std::cout << "trackAll[mateListItr->first].getMultiple(): " << trackAll[mateListItr->first].getMultiple() << std::endl;
            //std::cout << "*****" << std::endl;
            
            ++counter_1;
            if( mateListItr->second.getBT() == false &&                 //matepair is not seen
                trackAll.find(matePairItr->first) != trackAll.end() &&  //first mate is tracked
                trackAll[matePairItr->first].getMultiple() == 1 &&      //first mate seen once
                trackAll.find(mateListItr->first) != trackAll.end() &&  //second mate is tracked
                trackAll[mateListItr->first].getMultiple() == 1) {      //second mate is seen once
                ++counter_2;
                mateListItr->second.setBT(true);

                insert_size = matePair[matePairItr->first][mateListItr->first].getIS();
                min_allowed = -1 * (insertStdev * insert_size); // check int
                low_iz = insert_size + min_allowed; // check int
                up_iz = insert_size - min_allowed; // check int

                if(verbose) std::cout << "Pair read1Hash=" << matePairItr->first << " read2Hash=" << mateListItr->first << "\n";

                if(trackAll[matePairItr->first].getTig() != "" && trackAll[mateListItr->first].getTig() != "") { //double check if tig names not null
        ++counter_3;          
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
                        if(!kmer1.getOrient()){             // if kmer1 is forward
                            if(!kmer2.getOrient()){         // if kmer2 is forward
        			++counter_4;  
	                      distance = getDistance(insert_size, tigLength[kmer1.getTig()], kmer1.getStart(), kmer2.getStart());
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.getTig(), kmer2.getTig(), 0);
                                }
                            }else{                          // if kmer2 is reverse
				++counter_5;
                                distance = getDistance(insert_size, tigLength[kmer1.getTig()], kmer1.getStart(), tigLength[kmer2.getTig()] - kmer2.getEnd());
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.getTig(), kmer2.getTig(), 1);
                                }
                            }
                        }else{                              // if kmer1 is reverse
                            if(!kmer2.getOrient()){         // if kmer2 is forward
                                ++counter_6;
				distance = getDistance(insert_size, tigLength[kmer1.getTig()], tigLength[kmer1.getTig()] - kmer1.getEnd(), kmer2.getStart());
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.getTig(), kmer2.getTig(), 2);
                                }
                            }else{                          // if kmer2 is reverse
                                ++counter_7;
				distance = getDistance(insert_size, tigLength[kmer2.getTig()], kmer2.getEnd(), kmer1.getEnd());
                                if(distance > min_allowed  && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.getTig(), kmer2.getTig(), 3);
                                }
                            }
                        }
                    } else { // Clone, paired reads located on the same contig -- could be used to investigate misassemblies
                        if (verbose) std::cout << "Pair (" << matePairItr->first << " and " << mateListItr->first << ") located on same contig " << tig_a << " (" << A_length << " nt)\n";
                        uint64_t pet_size = 0;
                        if(A_start > B_start && (B_start < B_end) && (A_start > A_end)) {   // B --> <-- A
                            pet_size = A_start - B_start;
                            trackInsert += pet_size;
                            if(pet_size >= low_iz && pet_size <= up_iz) {
                                ct_ok_contig++;
                                if(ct_ok_contig_hash.find(insert_size) == ct_ok_contig_hash.end()) {
                                    ct_ok_contig_hash[insert_size] = 1;
                                } else {
                                    ct_ok_contig_hash[insert_size] = ct_ok_contig_hash[insert_size] + 1;
                                }
                            } else {
                                issuesFile << "Pairs unsatisfied in distance within a contig.  Pair (" << matePairItr->first << " - " << mateListItr->first << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << " CALCULATED DISTANCE APART: " << pet_size << "\n";
                                ct_iz_issues++;
                                if(ct_iz_issues_hash.find(insert_size) == ct_iz_issues_hash.end()) {
                                    ct_iz_issues_hash[insert_size] = 1;
                                } else {
                                    ct_iz_issues_hash[insert_size] = ct_iz_issues_hash[insert_size] + 1;
                                }
                            }
                        } else if(B_start > A_start && (B_start > B_end) && (A_start < A_end)) { // A --> <-- B
                            pet_size = B_start - A_start;
                            trackInsert += pet_size;
                            if(pet_size >= low_iz && pet_size <= up_iz) {
                                ct_ok_contig++;
                                if(ct_ok_contig_hash.find(insert_size) == ct_ok_contig_hash.end()) {
                                    ct_ok_contig_hash[insert_size] = 1;
                                } else {
                                    ct_ok_contig_hash[insert_size] = ct_ok_contig_hash[insert_size] + 1;
                                }
                            } else {
                                issuesFile << "Pairs unsatisfied in distance within a contig.  Pair (" << matePairItr->first << " - " << mateListItr->first << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << "\n";
                                ct_iz_issues++;
                                if(ct_iz_issues_hash.find(insert_size) == ct_iz_issues_hash.end()) {
                                    ct_iz_issues_hash[insert_size] = 1;
                                } else {
                                    ct_iz_issues_hash[insert_size] = ct_iz_issues_hash[insert_size] + 1;
                                }
                            }
                        } else {
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
                    ct_single++;
                    if(ct_single_hash.find(insert_size) == ct_single_hash.end()) {
                        ct_single_hash[insert_size] = 1;
                    } else {
                        ct_single_hash[insert_size] = ct_single_hash[insert_size] + 1;
                    }
                }
            } else { // if unseen
                // std::cout << "UNSEEN\n";
                if(matePair[matePairItr->first][mateListItr->first].getBT() == false) {
                    //std::cout << "UNSEEN getBT() increment ct\n";
                    ct_multiple++;
                }

            }
        } // pairing read b
    } // pairing read a

	std::cout << "counter_1: " << counter_1 << std::endl;
	std::cout << "counter_2: " << counter_2 << std::endl;
	std::cout << "counter_3: " << counter_3 << std::endl;
	std::cout << "counter_4: " << counter_4 << std::endl;
	std::cout << "counter_5: " << counter_5 << std::endl;
	std::cout << "counter_6: " << counter_6 << std::endl;
	std::cout << "counter_7: " << counter_7 << std::endl;
   
    
    // uint second_dim_size = 0;
    // uint third_dim_size = 0;
    // uint fourth_dim_size = 0;
    // for (auto it = matePair.begin(); it != matePair.end(); it++) {
    // //std::cout << *it << endl;
    //     second_dim_size += it->second.size();
    //     for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
    //         third_dim_size += it2->second.links;
    //         /*for (auto it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
    //             fourth_dim_size += it3->second.size(); 
    //         }*/
    //     }
    // } 

    // std::cout << "second dim size: " << second_dim_size << std::endl;
    // std::cout << "third dim size: " << third_dim_size << std::endl;
    // std::cout << "fourth dim size: " << fourth_dim_size << std::endl;
    

    uint pair_second_dim_size = 0;
    uint pair_third_dim_size = 0;
    uint pair_fourth_dim_size = 0;
    for (auto it = pair.begin(); it != pair.end(); it++) {
        //std::cout << *it << endl;
        //std::cout << "it->second.size(): " << it->second.size() << std::endl; 
        pair_second_dim_size += it->second.size();
        for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
        pair_third_dim_size += it2->second.size();
        for (auto it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
            pair_fourth_dim_size += it3->second.getLinks(); 
        }
        }
    } 
    std::cout << "second dim size: " << pair_second_dim_size << std::endl;
    std::cout << "third dim size: " << pair_third_dim_size << std::endl;
    std::cout << "fourth dim size: " << pair_fourth_dim_size << std::endl;

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

    
    //std::cout << "THESE ARE THE FILTERINGS:\n"<< "filter 1: "<< std::to_string(filter1) << "\n" << "filter 2: "<< std::to_string(filter2) << "\n" << "filter 3: "<< std::to_string(filter3) << "\n" << "filter 4: "<< std::to_string(filter4) << "\n";
    //std::cout << "THESE ARE THE COUNTERS:\n" << std::to_string(CheckCounterBase) << "\n0 " << std::to_string(Check0Counter) << "\n1 " <<  std::to_string(Check1Counter) << "\n2 " <<  std::to_string(Check2Counter) << "\n3 " <<  std::to_string(Check3Counter) << "\n4 " <<  std::to_string(Check4Counter) << "\n5 " <<  std::to_string(Check5Counter) << "\n6 " <<  std::to_string(Check6Counter) << "\n7 " <<  std::to_string(Check7Counter) << "\n8 " <<  std::to_string(Check8Counter) << "\n9 " <<  std::to_string(Check9Counter) << "\n10 " << std::to_string(Check10Counter) << "\n11 " << std::to_string(Check11Counter) << "\n12 " << std::to_string(Check12Counter) << "\n13 " << std::to_string(Check13Counter) << "\n14 " << std::to_string(Check14Counter) << "\n15 " << std::to_string(Check15Counter) << "\n16 " << std::to_string(Check16Counter) << "\n17 " << std::to_string(Check17Counter) << "\n18 " << std::to_string(Check18Counter) << "\n19" << std::to_string(Check19Counter) << "\n20 " << std::to_string(Check20Counter) << "\n21 " << std::to_string(Check21Counter) << "\n22 " << std::to_string(Check22Counter) << "\n23 " << std::to_string(Check23Counter) << "\n24 " << std::to_string(Check24Counter) << "\n25 " << std::to_string(Check25Counter) << "\n26 " << std::to_string(Check26Counter) << "\n";
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
    
    std::cout << "ct_both: " << ct_both << std::endl;
    sortInsertSize(ct_both_hash);
    std::unordered_map<uint64_t, uint64_t>::iterator itrIS;
    std::cout << "ct_both_hash map size: " << err.size() << "\n"; 
    for(itrIS = ct_both_hash.begin(); itrIS != ct_both_hash.end(); itrIS++) {
        std::cout <<  "--------k-mers separated by "<< itrIS->first << " bp (outer distance)--------\n";
        //int64_t maopt = -1 * (insertStdev * itrIS->first);
        //int64_t low_izopt = itrIS->first + maopt;
        //int64_t up_izopt = itrIS->first - maopt;
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
             << it.second << std::endl; 
    } 
} 

bool does_file_exist(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}
