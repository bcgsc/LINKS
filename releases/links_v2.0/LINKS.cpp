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


class BT_IS {
    private:
    uint64_t bt;
    uint64_t is;
    // change this to a map of hash to distance
    // uint64_t distance;
    
    public:
    BT_IS(uint64_t bt, uint64_t is) {
        this->bt = bt;
        this-> is = is;
        // secondLayerReads.insert({hash, distance});
    }

    void setBT(uint64_t bt) {
        this->bt = bt;
    }
    void setIS(uint64_t is) {
        this->is = is;
    }
    uint64_t getBT() {
        return this->bt;
    }
    uint64_t getIS() {
        return this->is;
    }
};

class Gaps_Links {
    private:
    uint64_t gaps;
    uint64_t links;
    std::string type;
    
    public:
    Gaps_Links(uint64_t gap, uint64_t link) {
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

    void addToGap(uint64_t distance) {
        this->gaps += distance;
    }

    void incrementLinks() {
        this->links++;
    }

    uint64_t getGaps() {
        return this->gaps;
    }
    uint64_t getLinks() {
        return this->links;
    }
};

class kMerInfo {
    private:
    std::string tig;
    uint64_t start;
    uint64_t end;
    uint64_t multiple;
    
    public:
    kMerInfo(){
        this->tig = "";
        this->start = 0;
        this->end = 0;
        this->multiple = 0;
    }
    kMerInfo(uint64_t start, uint64_t end){
        this->tig = "";
        this->start = start;
        this->end = end;
        this->multiple = 0;
    }
    kMerInfo(std::string tig, uint64_t start, uint64_t end){
        this->tig = tig;
        this->start = start;
        this->end = end;
        this->multiple = 0;
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
    uint64_t distances = 4000;
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

// Helper Methods
void printBloomStats(btllib::KmerBloomFilter& bloom, std::ostream& os);
long getFileSize(std::string filename);
void readContigs(
        std::string assemblyFile,
        std::map<const uint64_t, kMerInfo *>& trackAll,
        std::map<const uint64_t, std::map<const uint64_t, BT_IS *> > matePair,
        std::map<std::string, uint64_t>& tigLength,
        uint64_t k,
        uint64_t minSize,
        uint64_t hashFcts);
void pairContigs(
    std::string longReadsFile,
    std::map<const uint64_t, std::map<const uint64_t, BT_IS *> > matePair,
    std::map<const uint64_t, kMerInfo *> trackAll,
    std::map<std::string, uint64_t> tigLength,
    std::string issues,
    std::string distribution,
    uint64_t totalPairs,
    std::string tigpair_checkpoint,
    std::string simplepair_checkpoint,
    bool verbose,
    float insertStdev);
uint64_t getDistance(
    uint64_t insert_size,
    uint64_t length_i,
    uint64_t start_i,
    uint64_t start_j);


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



    uint64_t m = ceil((-1 * (double)bfElements * log(linksArgParser->fpr)) / (log(2) * log(2)));
    uint64_t rem = 64 - (m % 64);
    m = ((uint64_t)(m / 8) + 1) * 8;
    std::cout << "HASHES CALC: " << std::to_string(((double)m / bfElements)) << " second: " << std::to_string(((double)m / bfElements) * log(2)) << "\n";
    uint64_t hashFct = floor(((double)m / bfElements) * log(2));
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
    std::map<const uint64_t, std::map<const uint64_t, BT_IS *> > matePair;

    btllib::BloomFilter& filtering = myFilter.get_bloom_filter();
    btllib::SeqReader longReader(linksArgParser->longFile);
    uint64_t counter = 0;
    uint64_t hits = 0;
    for (btllib::SeqReader::Record record; (record = longReader.read());) {
        btllib::NtHash nthash(record.seq, linksArgParser->k, hashFct);
        btllib::NtHash nthashLead(record.seq, linksArgParser->k, hashFct, linksArgParser->distances + linksArgParser->k);
        for (size_t i = 0; nthash.roll() && nthashLead.roll(); i+=linksArgParser->step) {
            // roll for the number of steps
            counter++;
            if(filtering.contains(nthash.hashes()) && filtering.contains(nthashLead.hashes())) {
                hits++;
                // If this hash exists in matePair, add the read to the second layer of instead of making a new entry
                if(matePair[nthash.hashes()[0]].find(nthashLead.hashes()[0]) == matePair[nthash.hashes()[0]].end()) {
                    matePair[nthash.hashes()[0]][nthashLead.hashes()[0]] = new BT_IS(0, linksArgParser->distances);
                } else {
                    // LongReadKmer * leadPair = new LongReadKmer(nthashLead.hashes()[0], linksArgParser->distances);
                    // Check for existence
                    // frag_dist is an array of distances
                    matePair[nthash.hashes()[0]][nthashLead.hashes()[0]]->setBT(1);
                    matePair[nthash.hashes()[0]][nthashLead.hashes()[0]]->setIS(linksArgParser->distances);
                }
            }
        }
    }
    std::cout << hits << " match percentage: % " << "matePair size: " << (double)matePair.size()<< "   " << (double)matePair.size()/counter * 100.0 << " counter: " << counter << " \n";
    
    
    std::map<const uint64_t, kMerInfo *> trackAll;
    std::map<std::string, uint64_t> tigLength;
    std::cout << "\n\n=>Reading sequence contigs (to scaffold), tracking k-mer positions :" << dt << "\n";
    // Read contigs to find where the long read kmers belong in
    readContigs(linksArgParser->assemblyFile, trackAll, matePair, tigLength, linksArgParser->k, linksArgParser->minSize, hashFct);
    std::cout << " The resulting trackAll map size is: " << trackAll.size() << "\n\n";
    pairContigs(
        linksArgParser->longFile,
        matePair,
        trackAll,
        tigLength,
        issues,
        distribution,
        tigLength.size(),
        tigpair_checkpoint,
        simplepair_checkpoint,
        linksArgParser->verbose,
        linksArgParser->insertStdev);
    return 0;
}

void kmerizeContig( std::string seq, 
                    std::map<const uint64_t, kMerInfo *>& trackAll,
                    std::map<const uint64_t, std::map<const uint64_t, BT_IS *> > matePair,
                    uint64_t k,
                    std::string head,
                    uint64_t hashFcts,
                    uint64_t step) {

    btllib::NtHash ntHashContig(seq, k, hashFcts);
    int counter = 0;
    std::cout << "hashFct in kmerizeContig: " << hashFcts << "\n";
    for (size_t i = 0; ntHashContig.roll(); i+=step) {
        // roll for the steps
        // for (int j < step)
        std::cout << "rollin\n";
        if(matePair.find(ntHashContig.hashes()[0]) == matePair.end()) {
            if(trackAll.find(ntHashContig.hashes()[0]) == trackAll.end()) {
                trackAll[ntHashContig.hashes()[0]] = new kMerInfo(head, i, i + k);
            } else {
                std::cout << "kmer found in trackall! Increment multiple\n";
                trackAll[ntHashContig.hashes()[0]]->incrementMultiple();
            }
        }
        counter++;
    }
    std::cout << "matePair size:" << matePair.size();
    std::cout << "\n\n\n";
}

void readContigs(
        std::string assemblyFile,
        std::map<const uint64_t, kMerInfo *>& trackAll,
        std::map<const uint64_t, std::map<const uint64_t, BT_IS *> > matePair,
        std::map<std::string, uint64_t>& tigLength,
        uint64_t k,
        uint64_t minSize,
        uint64_t hashFcts) {
    std::cout << "hashFct in readContig: " << hashFcts << "\n";
    uint64_t cttig = 0;
    btllib::SeqReader contigReader(assemblyFile);
    for (btllib::SeqReader::Record record; (record = contigReader.read());) {
        tigLength.insert({record.name, record.seq.length()});
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

void pairContigs(
    std::string longReadsFile,
    std::map<const uint64_t, std::map<const uint64_t, BT_IS *> > matePair,
    std::map<const uint64_t, kMerInfo *> trackAll,
    std::map<std::string, uint64_t> tigLength,
    std::string issues,
    std::string distribution,
    uint64_t totalPairs,
    std::string tigpair_checkpoint,
    std::string simplepair_checkpoint,
    bool verbose,
    float insertStdev) {

    uint64_t ct_illogical = 0, ct_illogical_hash = 0, ct_ok_contig = 0, ct_ok_contig_hash = 0, ct_ok_pairs = 0, ct_ok_pairs_hash = 0, ct_problem_pairs = 0, ct_problem_pairs_hash = 0, ct_iz_issues = 0, ct_iz_issues_hash = 0, ct_single = 0, ct_single_hash = 0, ct_multiple = 0, ct_both = 0, trackInsert = 0;

    // Mapping of tiga_head -> insertSize -> tigb_head -> links & gaps
    std::map<std::string, std::map<uint64_t, std::map<std::string, Gaps_Links *> > > pair;
    std::map<std::string, std::map<std::string, Gaps_Links*> >simplePair;
    std::map<std::string, Gaps_Links*> err;
    std::string order1;
    std::string order2;
    if(verbose) std::cout << "Pairing contigs...\n";
    
    std::map<const uint64_t, BT_IS *>::iterator mateListItr;
    std::map<const uint64_t, std::map<const uint64_t, BT_IS *> >::iterator matePairItr;
    for(matePairItr = matePair.begin(); matePairItr != matePair.end(); matePairItr++) {
        for(mateListItr = matePairItr->second.begin(); mateListItr != matePairItr->second.end(); mateListItr++) {
            if(mateListItr->second == 0 && trackAll[matePairItr->first]->getMultiple() == 1 && trackAll[mateListItr->first]->getMultiple() == 1) { // This has little if no effect, but negative for some odd reason
                // below indicates this specific pair has been seen (bt = 1)
                mateListItr->second->setBT(1);
                uint64_t insert_size = matePair[matePairItr->first][mateListItr->first]->getIS();
                int min_allowed = -1 * (insertStdev * insert_size); // check int
                int low_iz = insert_size + min_allowed; // check int
                int up_iz = insert_size - min_allowed; // check int
                if(verbose) std::cout << "Pair read1Hash=" << matePairItr->first << " read2Hash=" << mateListItr->first;

                if(trackAll[matePairItr->first]->getTig() != "" && trackAll[mateListItr->first]->getTig() != "") {
                    ct_both++;
                    // $ct_both_hash->{$insert_size}++;
                    std::string tig_a = trackAll[matePairItr->first]->getTig();
                    std::string tig_b = trackAll[mateListItr->first]->getTig();

                    std::string ftig_a = "f" + tig_a;
                    std::string ftig_b = "f" + tig_b;

                    std::string rtig_a = "r" + tig_a;
                    std::string rtig_b = "r" + tig_b;

                    uint64_t A_length = tigLength[tig_a];
                    uint64_t A_start = trackAll[matePairItr->first]->getStart();
                    uint64_t A_end = trackAll[matePairItr->first]->getEnd();

                    uint64_t B_length = tigLength[tig_b];
                    uint64_t B_start = trackAll[mateListItr->first]->getStart();
                    uint64_t B_end = trackAll[mateListItr->first]->getEnd();

                    if(tig_a != tig_b) { // paired reads located on <> contigs
                        //Determine most likely possibility
                        if(trackAll[matePairItr->first]->getStart() < trackAll[matePairItr->first]->getEnd()) {
                            if(trackAll[mateListItr->first]->getEnd() < trackAll[mateListItr->first]->getStart()) { // -> <- :::  A-> <-B  /  rB -> <- rA
                                uint64_t distance = getDistance(insert_size, A_length, A_start, B_start);
                                if(verbose) std::cout << "A-> <-B  WITH " << tig_a << "-> <- " << tig_b << " GAP " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " " << std::to_string(B_start-B_end) << " Alen, Astart,Bstart\n";
                                if(distance > min_allowed) {
                                    uint64_t isz = distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 1000; // distance categories
                                    if(pair.find(ftig_a) == pair.end() || pair[ftig_a].find(isz) == pair[ftig_a].end() || pair[ftig_a][isz].find(rtig_b) == pair[ftig_a][isz].end()) {
                                        pair[ftig_a][isz][rtig_b] = new Gaps_Links();
                                    } else {
                                        pair[ftig_a][isz][rtig_b]->addToGap(distance);
                                        pair[ftig_a][isz][rtig_b]->incrementLinks();
                                    }
                                    if(pair.find(rtig_b) == pair.end() || pair[rtig_b].find(isz) == pair[rtig_b].end() || pair[rtig_b][isz].find(rtig_a) == pair[rtig_b][isz].end()) {
                                        pair[rtig_b][isz][rtig_a] = new Gaps_Links();
                                    } else {
                                        pair[rtig_b][isz][rtig_a]->addToGap(distance);
                                        pair[rtig_b][isz][rtig_a]->incrementLinks();
                                    }
                                    if(tig_a < tig_b) {
                                        order1 = tig_a;
                                        order2 = tig_b;
                                    } else {
                                        order1 = tig_b;
                                        order2 = tig_a;
                                    }
                                    // Check if exists
                                    simplePair[order1][order2] = new Gaps_Links("11");
                                    simplePair[order1][order2]->incrementLinks();
                                    simplePair[order1][order2]->addToGap(distance);
                                    
                                    ct_ok_pairs++;
                                    ct_ok_pairs_hash += distance;
                                } else {
                                    std::string err_pair = ftig_a + "-" + ftig_b;
                                    if(err.find(err_pair) == err.end()) {
                                        err[err_pair] = new Gaps_Links(distance, 1);
                                    } else {
                                        err[err_pair]->addToGap(distance);
                                        err[err_pair]->incrementLinks();
                                    }
                                    ct_problem_pairs++;
                                    ct_problem_pairs_hash += distance;
                                    std::cout << "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#" << tig_a << " -> " << std::to_string(distance) << " <- tig#"<< tig_b << ", A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") CALCULATED DISTANCE APART: " << distance << " < " << min_allowed << "\n";
                                }
                            } else { // -> -> ::: A-> <-rB  / B-> <-rA 
                                uint64_t rB_start = B_length - B_start;
                                uint64_t distance = getDistance(insert_size, A_length, A_start, rB_start);
                                if(verbose) std::cout << "A-> <-rB  WITH " << tig_a << "-> <- " << tig_b << " GAP " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " " << std::to_string(B_start-B_end) << " Alen, Astart,rBstart\n";
                                if(distance >= min_allowed) {
                                    uint64_t isz = distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 1000; // distance categories
                                    if(pair.find(ftig_a) == pair.end() || pair[ftig_a].find(isz) == pair[ftig_a].end() || pair[ftig_a][isz].find(rtig_b) == pair[ftig_a][isz].end()) {
                                        pair[ftig_a][isz][rtig_b] = new Gaps_Links();
                                    } else {
                                        pair[ftig_a][isz][rtig_b]->addToGap(distance);
                                        pair[ftig_a][isz][rtig_b]->incrementLinks();
                                    }
                                    if(pair.find(ftig_b) == pair.end() || pair[ftig_b].find(isz) == pair[ftig_b].end() || pair[ftig_b][isz].find(rtig_a) == pair[ftig_b][isz].end()) {
                                        pair[ftig_b][isz][rtig_a] = new Gaps_Links();
                                    } else {
                                        pair[ftig_b][isz][rtig_a]->addToGap(distance);
                                        pair[ftig_b][isz][rtig_a]->incrementLinks();
                                    }
                                    if(tig_a < tig_b) {
                                        order1 = tig_a;
                                        order2 = tig_b;
                                    } else {
                                        order1 = tig_b;
                                        order2 = tig_a;
                                    }
                                    simplePair[order1][order2] = new Gaps_Links("10");
                                    simplePair[order1][order2]->incrementLinks();
                                    simplePair[order1][order2]->addToGap(distance);
                                    
                                    ct_ok_pairs++;
                                    ct_ok_pairs_hash += distance;
                                } else {
                                    std::string err_pair = ftig_a + "-" + rtig_b;
                                    if(err.find(err_pair) == err.end()) {
                                        err[err_pair] = new Gaps_Links(distance, 1);
                                    } else {
                                        err[err_pair]->addToGap(distance);
                                        err[err_pair]->incrementLinks();
                                    }
                                    ct_problem_pairs++;
                                    ct_problem_pairs_hash += distance;
                                    std::cout << "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#" << tig_a << " -> " << std::to_string(distance) << " <- tig#r."<< tig_b << ", A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") CALCULATED DISTANCE APART: " << distance << " < " << min_allowed << "\n";
                                }
                            }
                        } else {
                            // if ({read_b}{'end'} > {$read_b}{'start'}
                            if(trackAll[mateListItr->first]->getEnd() > trackAll[mateListItr->first]->getStart()) {
                                uint64_t distance = getDistance(insert_size, B_length, B_start, A_start);
                                if(verbose) std::cout << "B-> <-A  WITH " << tig_b << "-> <- " << tig_a << " GAP " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " " << std::to_string(B_start-B_end) << " Blen, Bstart,Astart\n";
                                if(distance >= min_allowed) {
                                    uint64_t isz = distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 1000; // distance categories
                                    if(pair.find(ftig_b) == pair.end() || pair[ftig_b].find(isz) == pair[ftig_b].end() || pair[ftig_b][isz].find(ftig_a) == pair[ftig_b][isz].end()) {
                                        pair[ftig_b][isz][ftig_a] = new Gaps_Links();
                                    } else {
                                        pair[ftig_b][isz][ftig_a]->addToGap(distance);
                                        pair[ftig_b][isz][ftig_a]->incrementLinks();
                                    }
                                    if(pair.find(rtig_a) == pair.end() || pair[rtig_a].find(isz) == pair[rtig_a].end() || pair[rtig_a][isz].find(rtig_b) == pair[rtig_a][isz].end()) {
                                        pair[rtig_a][isz][rtig_b] = new Gaps_Links();
                                    } else {
                                        pair[rtig_a][isz][rtig_b]->addToGap(distance);
                                        pair[rtig_a][isz][rtig_b]->incrementLinks();
                                    }
                                    if(tig_a < tig_b) {
                                        order1 = tig_a;
                                        order2 = tig_b;
                                    } else {
                                        order1 = tig_b;
                                        order2 = tig_a;
                                    }
                                    simplePair[order1][order2] = new Gaps_Links("11");
                                    simplePair[order1][order2]->incrementLinks();
                                    simplePair[order1][order2]->addToGap(distance);
                                    
                                    ct_ok_pairs++;
                                    ct_ok_pairs_hash += distance;
                                } else {
                                    std::string err_pair = ftig_b + "-" + ftig_a;
                                    if(err.find(err_pair) == err.end()) {
                                        err[err_pair] = new Gaps_Links(distance, 1);
                                    } else {
                                        err[err_pair]->addToGap(distance);
                                        err[err_pair]->incrementLinks();
                                    }
                                    ct_problem_pairs++;
                                    ct_problem_pairs_hash += distance;
                                    std::cout << "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#" << tig_b << " -> " << std::to_string(distance) << " <- tig#"<< tig_a << ", B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") CALCULATED DISTANCE APART: " << distance << " < " << min_allowed << "\n";
                                }
                            }  else { // <- <-  :::  rB-> <-A / rA-> <-B
                                uint64_t rB_start = B_length - B_start;
                                uint64_t distance = getDistance(insert_size, B_length, rB_start, A_start);
                                if(verbose) std::cout << "rB-> <-A  WITH r." << tig_b << "-> <- " << tig_a << " GAP " << std::to_string(distance) << " A=" << std::to_string(A_length) << " " << std::to_string(A_start - A_end) << " B= " << B_length << " " << std::to_string(B_start-B_end) << " Blen, rBstart,Astart\n";
                                if(distance >= min_allowed) {
                                    uint64_t isz = distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 1000; // distance categories
                                    if(pair.find(rtig_b) == pair.end() || pair[rtig_b].find(isz) == pair[rtig_b].end() || pair[rtig_b][isz].find(rtig_b) == pair[ftig_a][isz].end()) {
                                        pair[rtig_b][isz][ftig_a] = new Gaps_Links();
                                    } else {
                                        pair[rtig_b][isz][ftig_a]->addToGap(distance);
                                        pair[rtig_b][isz][ftig_a]->incrementLinks();
                                    }
                                    if(pair.find(rtig_a) == pair.end() || pair[rtig_a].find(isz) == pair[rtig_a].end() || pair[rtig_a][isz].find(ftig_b) == pair[rtig_a][isz].end()) {
                                        pair[rtig_a][isz][ftig_b] = new Gaps_Links();
                                    } else {
                                        pair[rtig_a][isz][ftig_b]->addToGap(distance);
                                        pair[rtig_a][isz][ftig_b]->incrementLinks();
                                    }
                                    if(tig_a < tig_b) {
                                        order1 = tig_a;
                                        order2 = tig_b;
                                    } else {
                                        order1 = tig_b;
                                        order2 = tig_a;
                                    }
                                    simplePair[order1][order2] = new Gaps_Links("01");
                                    simplePair[order1][order2]->incrementLinks();
                                    simplePair[order1][order2]->addToGap(distance);
                                    
                                    ct_ok_pairs++;
                                    ct_ok_pairs_hash += distance;
                                } else {
                                    std::string err_pair = rtig_b + "-" + ftig_a;
                                    if(err.find(err_pair) == err.end()) {
                                        err[err_pair] = new Gaps_Links(distance, 1);
                                    } else {
                                        err[err_pair]->addToGap(distance);
                                        err[err_pair]->incrementLinks();
                                    }
                                    ct_problem_pairs++;
                                    ct_problem_pairs_hash += distance;
                                    std::cout << "Pairs unsatisfied in distance within a contig pair.  rB-> <-A  WITH tig#r." << tig_b << " -> " << std::to_string(distance) << " <- tig#"<< tig_a << ", B=" << B_length << " nt (start:" << B_start << ", end:" << B_end << ") A=" << A_length << " nt (start:" << A_start << ", end:" << A_end << ") CALCULATED DISTANCE APART: " << distance << " < " << min_allowed << "\n";
                                }
                            }
                        }
                    } else { // Clone, paired reads located on the same contig -- could be used to investigate misassemblies
                        if (verbose) std::cout << "Pair (" << trackAll[matePairItr->first] << "-> and " << trackAll[mateListItr->first] << ") located on same contig " << tig_a << " (" << A_length << " nt)\n";
                        uint64_t pet_size = 0;
                        if(A_start > B_start && (B_start < B_end) && (A_start > A_end)) {   // B --> <-- A
                            pet_size = A_start - B_start;
                            trackInsert += pet_size;
                            if(pet_size >= low_iz && pet_size <= up_iz) {
                                ct_ok_contig++;
                                ct_ok_contig_hash += insert_size;
                            } else {
                                std::cout <<"Pairs unsatisfied in distance within a contig.  Pair (" << trackAll[matePairItr->first] << " - " << trackAll[mateListItr->first] << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << " CALCULATED DISTANCE APART: " << pet_size << "\n";
                                ct_iz_issues++;
                                ct_iz_issues_hash += insert_size;
                            }
                        } else if(B_start > A_start && (B_start > B_end) && (A_start < A_end)) { // A --> <-- B
                            pet_size = B_start - A_start;
                            trackInsert += pet_size;
                            if(pet_size >= low_iz && pet_size <= up_iz) {
                                ct_ok_contig++;
                                ct_ok_contig_hash += insert_size;
                            } else {
                                std::cout <<"Pairs unsatisfied in distance within a contig.  Pair (" << trackAll[matePairItr->first] << " - " << trackAll[mateListItr->first] << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << "\n";
                                ct_iz_issues++;
                                ct_iz_issues_hash += insert_size;
                            }
                        } else {
                            ct_illogical++;
                            ct_illogical_hash += insert_size;
                            std::cout << "Pairs unsatisfied in pairing logic within a contig.  Pair (" << trackAll[matePairItr->first] << " - " << trackAll[mateListItr->first] << ") on contig" << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << "\n";
                        }
                    }
                } else { // both pairs assembled
                    ct_single++;
                    ct_single_hash += insert_size;
                }
            } else { // if unseen
                if(matePair[matePairItr->first][mateListItr->first]->getBT() == 0) {
                    ct_multiple++;
                }
            }
        } // pairing read b
    } // pairing read a

// Summary of the contig pair issues
// std::cout << "------------- Putative issues with contig pairing - Summary  ----------------\n";
// for(errorPair = matePair.begin(); matePairItr != matePair.end(); matePairItr++) {
    uint64_t satisfied = ct_ok_pairs + ct_ok_contig;
    uint64_t unsatisfied = ct_problem_pairs + ct_iz_issues + ct_illogical;
    uint64_t ct_both_reads = ct_both * 2;

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
}

uint64_t getDistance(uint64_t insert_size, uint64_t length_i, uint64_t start_i, uint64_t start_j) {

    // L  ------  --------- R
    // i    ->        <-    j
    //      ....  ......    insert_span
    //      ============    insert_size

   uint64_t insert_span = (length_i - start_i) + start_j;
   uint64_t gap_or_overlap = insert_size - insert_span;

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








//**********SORTER*********
// bool cmp(pair<string, int>& a, 
//          pair<string, int>& b) 
// { 
//     return a.second < b.second; 
// } 
  
// // Function to sort the map according 
// // to value in a (key-value) pairs 
// void sort(map<string, int>& M) 
// { 
  
//     // Declare vector of pairs 
//     vector<pair<string, int> > A; 
  
//     // Copy key-value pair from Map 
//     // to vector of pairs 
//     for (auto& it : M) { 
//         A.push_back(it); 
//     } 
  
//     // Sort using comparator function 
//     sort(A.begin(), A.end(), cmp); 
  
//     // Print the sorted value 
//     for (auto& it : A) { 
  
//         cout << it.first << ' '
//              << it.second << endl; 
//     } 
// } 
