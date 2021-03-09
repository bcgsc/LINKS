// using namespace btllib;
#include <string>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <cmath>
#include <regex>
#include <map>
#include <ctime>


#include "btllib/bloom_filter.hpp"
#include "btllib/seq_reader.hpp"

//Globals
#define BASE_TEN 10
std::string version = "2.0";


class LongReadKmer {
    private:
    uint64_t hash;
    uint64_t distance;
    
    public:
    LongReadKmer(const uint64_t hash, int distance){
        this->hash = hash;
        this->distance = distance;
    }
};

class kMerInfo {
    private:
    std::string head;
    uint64_t start;
    uint64_t end;
    uint64_t multiple;
    
    public:
    kMerInfo(){
        this->head = "";
        this->start = 0;
        this->end = 0;
        this->multiple = 0;
    }
    kMerInfo(uint64_t start, uint64_t end){
        this->head = "";
        this->start = start;
        this->end = end;
        this->multiple = 0;
    }
    kMerInfo(std::string head, uint64_t start, uint64_t end){
        this->head = head;
        this->start = start;
        this->end = end;
        this->multiple = 0;
    }

    void setHead(std::string head) {
        this->head = head;
    }

    void setStart(uint64_t start) {
        this->start = start;
    }

    void setEnd(uint64_t end) {
        this->end = end;
    }

    void incrementMultiple() {
        this->multiple += 1;
    }

    std::string getHead() {
        return this->head;
    }

    uint64_t getStart() {
        return this->start;
    }

    uint64_t getEnd() {
        return this->end;
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

    public:
    std::string assemblyFile;
    std::string longFile;
    unsigned distances = 4000;
    unsigned k = 15;
    bool verbose;
    unsigned minLinks = 5;
    unsigned minSize = 500;
    float maxLinkRatio = 0.3;
    unsigned step = 2;
    // Added for MPET
    unsigned readLength;         // MPET
    float insertStdev = 0.1;      // MPET (Need to adjust to a wider-range of distances when dealing with MPET) 
    std::string baseName;   // When set, this will override the MPET-induced changes on -e
    unsigned offset = 0;
    std::string bfFile;
    float fpr = 0.001;
    unsigned bfoff = 0;
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
                    distances = strtoul(optarg, &end, BASE_TEN);
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
                    << "  -d " << distances << "\n" 
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

void printBloomStats(btllib::KmerBloomFilter& bloom, std::ostream& os);
long getFileSize(std::string filename);
void readContigs(
        std::string assemblyFile,
        std::map<const uint64_t, kMerInfo *>& trackAll,
        std::map<const uint64_t, LongReadKmer *> matePair,
        unsigned k,
        unsigned minSize,
        unsigned hashFcts);

int main(int argc, char** argv) {
    // Get todays date
    time_t now = time(0);
   
    // convert now to string form
    char* dt = ctime(&now);
    // Parse command line arguments
    InputParser* linksArgParser = new InputParser(argc, argv);

    //Set Bloom Filter element number based on the size of the assembly file (1 byte = 1 character)
    std::string path = linksArgParser->assemblyFile;
    long bfElements = getFileSize(path);

    // Checking validity of input assemble file
    if(bfElements == -1){
        std::cout << std::to_string(bfElements);
        std::cout << "Invalid file: " << linksArgParser->assemblyFile << " -- fatal\n";
        return -1;
    }

    // Naming output files
    if (linksArgParser->baseName == "") {
        linksArgParser->baseName =  linksArgParser->assemblyFile + ".scaff_s-" + 
                                    linksArgParser->longFile + "_d" + 
                                    std::to_string(linksArgParser->distances) + 
                                    "_k" + std::to_string(linksArgParser->k) + 
                                    "_e" + std::to_string(linksArgParser->insertStdev) +
                                    "_l" + std::to_string(linksArgParser->minLinks) +
                                    "_a" + std::to_string(linksArgParser->maxLinkRatio) +
                                    "_z" + std::to_string(linksArgParser->minSize) +
                                    "_t" + std::to_string(linksArgParser->step) +
                                    "_o" + std::to_string(linksArgParser->offset) +
                                    "_r-"+ linksArgParser->bfFile +
                                    "_p" + std::to_string(linksArgParser->fpr) +
                                    "_x" + std::to_string(linksArgParser->bfoff) +
                                    "_m" + std::to_string(linksArgParser->readLength);
        pid_t pid_num = getppid();
        linksArgParser->baseName += "_pid" + std::to_string(pid_num);
    }
    std::cout << linksArgParser->baseName;
    std::string outlog = linksArgParser->baseName + ".log";
    std::string scaffold = linksArgParser->baseName + ".scaffolds";
    std::string issues = linksArgParser->baseName + ".pairing_issues";
    std::string distribution = linksArgParser->baseName + ".pairing_distribution.csv";
    std::string bfout = linksArgParser->baseName + ".bloom";
    std::string graph = linksArgParser->baseName + ".gv";
    std::string numnamecorr = linksArgParser->baseName + ".assembly_correspondence.tsv";
    std::string tigpair_checkpoint = linksArgParser->baseName + ".tigpair_checkpoint.tsv"; // add a checkpoint file, prevent re-running LINKS from scratch if crash
    std::string simplepair_checkpoint = linksArgParser->baseName + ".simplepair_checkpoint.tsv"; // add a checkpoint file, prevent re-running LINKS from scratch if crash
    
    if(freopen(outlog.c_str(), "w", stderr ) == NULL) {
        std::cout << "Can't write to " << outlog << " -- fatal\n";
        return -1;
    }
    //---------------------------
    // Initialization message
    std::string initMessage = "\nRunning: " + version + 
                                "\n-f " + linksArgParser->assemblyFile +
                                "\n-s " + linksArgParser->longFile + 
                                "\n-m " + std::to_string(linksArgParser->readLength) + 
                                "\n-d " + std::to_string(linksArgParser->distances) + 
                                "\n-k " + std::to_string(linksArgParser->k) + 
                                "\n-e " + std::to_string(linksArgParser->insertStdev) +
                                "\n-l " + std::to_string(linksArgParser->minLinks) +
                                "\n-a " + std::to_string(linksArgParser->maxLinkRatio) +
                                "\n-t " + std::to_string(linksArgParser->step) + 
                                "\n-o " + std::to_string(linksArgParser->offset) +
                                "\n-z " + std::to_string(linksArgParser->minSize) +
                                "\n-b " + linksArgParser->baseName +
                                "\n-r " + linksArgParser->bfFile +
                                "\n-p " + std::to_string(linksArgParser->fpr) +
                                "\n-x " + std::to_string(linksArgParser->bfoff) +
                                "\n\n----------------- Verifying files -----------------\n\n";
    
    std::cout << initMessage;
    std::cerr << initMessage;
    //-----------------------------



    unsigned long m = ceil((-1 * (double)bfElements * log(linksArgParser->fpr)) / (log(2) * log(2)));
    unsigned rem = 64 - (m % 64);
    m = ((unsigned)(m / 8) + 1) * 8;
    std::cout << "HASHES CALC: " << std::to_string(((double)m / bfElements)) << " second: " << std::to_string(((double)m / bfElements) * log(2)) << "\n";
    unsigned hashFct = floor(((double)m / bfElements) * log(2));
    std::cout << "- Number of bfElements: " << bfElements << "\n"
                << "- Input file path: " << path << "\n"
                << "- Input file: " << linksArgParser->assemblyFile << "\n"
                << "- kmersize: " << linksArgParser->k << "\n"
                << "- m: " << m << "\n"
                << "- fpr: " << linksArgParser->fpr << "\n"
                << "- hashFct: " << hashFct << "\n";

    // std::cout << "- Filter output file : " << outFileBf << "\n";
    std::cout << "- Filter output file : " << linksArgParser->k << "\n";
    btllib::KmerBloomFilter myFilter(m/8, hashFct, linksArgParser->k);
    btllib::SeqReader assemblyReader(linksArgParser->assemblyFile);
    int builder = 0;
    for (btllib::SeqReader::Record record; (record = assemblyReader.read());) {
        if(builder % 100000 == 0) {
            std::cout << "reading... builder: " << builder << "\n";
        }
                builder++;

        myFilter.insert(record.seq);
    }
    printBloomStats(myFilter, std::cout);
    // myFilter.storeFilter(outfile);

    // k-merize long reads
    std::map<const uint64_t, LongReadKmer *> matePair;

    btllib::BloomFilter& filtering = myFilter.get_bloom_filter();
    btllib::SeqReader longReader(linksArgParser->longFile);
    unsigned long counter = 0;
    unsigned long hits = 0;
    for (btllib::SeqReader::Record record; (record = longReader.read());) {
        btllib::NtHash nthash(record.seq, linksArgParser->k, hashFct);
        btllib::NtHash nthashLead(record.seq, linksArgParser->k, hashFct, linksArgParser->distances + linksArgParser->k);
        for (size_t i = 0; nthash.roll() && nthashLead.roll(); i+=linksArgParser->step) {
            counter++;
            if(filtering.contains(nthash.hashes()) && filtering.contains(nthashLead.hashes())) {
                hits++;
                LongReadKmer * leadPair = new LongReadKmer(nthashLead.hashes()[0], linksArgParser->distances);
                matePair.insert({nthash.hashes()[0], leadPair});
            }
        }
    }
    std::cout << hits << " match percentage: % " << "matePair size: " << (double)matePair.size()<< "   " << (double)matePair.size()/counter * 100.0 << " counter: " << counter << " \n";
    
    
    std::map<const uint64_t, kMerInfo *> trackAll;
    std::cout << "\n\n=>Reading sequence contigs (to scaffold), tracking k-mer positions :" << dt << "\n";

    std::cout << "hashFct in MAIN: " << hashFct << "\n";
    // Read contigs to find where the long read kmers belong in
    readContigs(linksArgParser->assemblyFile, trackAll, matePair, linksArgParser->k, linksArgParser->minSize, hashFct);

    std::cout << " The resulting trackAll map size is: " << trackAll.size() << "\n\n";

    return 0;
}

void kmerizeContig( std::string seq, 
                    std::map<const uint64_t, kMerInfo *>& trackAll,
                    std::map<const uint64_t, LongReadKmer *> matePair,
                    unsigned k,
                    std::string head,
                    unsigned hashFcts,
                    unsigned step) {

    btllib::NtHash ntHashContig(seq, k, hashFcts);
    int counter = 0;
    std::cout << "hashFct in kmerizeContig: " << hashFcts << "\n";
    for (size_t i = 0; ntHashContig.roll(); i+=step) {
        // roll for the steps
        // for (int j < step)
        if(matePair.find(ntHashContig.hashes()[0])!= matePair.end()) {
            // std::cout << "kmer not found in long reads! Insert it\n";
            trackAll.insert({ntHashContig.hashes()[0], new kMerInfo(head, i, i + k)});
            // std::cout << "size: " << trackAll.size() << "\n";
        } else {
            // std::cout << "kmer found in matePair! Increment multiple\n";
            trackAll.find(ntHashContig.hashes()[0])->second->incrementMultiple();
        }
        counter++;
    }
    std::cout << "matePair size:" << matePair.size();
    std::cout << "\n\n\n";
}

void readContigs(
        std::string assemblyFile,
        std::map<const uint64_t, kMerInfo *>& trackAll,
        std::map<const uint64_t, LongReadKmer *> matePair,
        unsigned k,
        unsigned minSize,
        unsigned hashFcts) {
    std::cout << "hashFct in readContig: " << hashFcts << "\n";
    uint64_t cttig = 0;
    btllib::SeqReader contigReader(assemblyFile);
    for (btllib::SeqReader::Record record; (record = contigReader.read());) {
        cttig++;
        std::cout << "\r" << cttig;
        if(record.seq.length() >= minSize) {
            std::cout << "Kmerizing contig\n";
            // std::cout << "seq:\n" << record.seq << "\n";
            std::cout << "k:\n" << k << "\n";
            std::cout << "hashFcts:\n" << hashFcts << "\n";
            // std::cout << "head:\n" << record.name << "\n";
            kmerizeContig(record.seq, trackAll, matePair, k, record.name, hashFcts, cttig);
        }
    }
}

//my ($longfile,$matepair,$track,$tig_length,$issues,$distribution,$totalpairs,$tigpair_checkpoint,$simplepair_checkpoint,$verbose) = @_;
// void pairContigs() {

// }

void printBloomStats(btllib::KmerBloomFilter& bloom, std::ostream& os)
{
	os << "Bloom filter stats:"
	   << "\n\t#counters               = " << bloom.get_occupancy()//getFilterSize()
	   << "\n\t#size (B)               = " << bloom.get_bytes()
	   << "\n\tpopcount                = " << bloom.get_pop_cnt()
	   << "\n\tFPR                     = " << 100.f * bloom.get_fpr() << "%"
	   << "\n";
}

long getFileSize(std::string filename)
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

