#ifndef LINKS_LINKS_HPP
#define LINKS_LINKS_HPP

#include "btllib/bloom_filter.hpp"
#include "btllib/nthash.hpp"
#include "btllib/order_queue.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/status.hpp"
#include "btllib/util.hpp"
#include "InputParser.hpp"

#include <unordered_map>
#include <unordered_set>

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <vector>
#include <mutex>
#include <shared_mutex>
#include <chrono>



class LINKS
{
  public:

    LINKS(InputParser* inputParser);
    void init_bloom_filter();
    void start_read_fasta();
    void start_read_contig();
    ~LINKS();

    InputParser* inputParser;

    std::string assemblyFile;
    std::string fofFile;
    std::vector<uint64_t> distances = {2000,4000};
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
    uint64_t bfOff = 0;

    static const size_t MAX_SIMULTANEOUS_INDEXLRS = 128;

    struct Read
    {
      Read() {}

      Read(size_t num, std::string id, std::string comment, std::string seq)
        : num(num)
        , id(std::move(id))
        , comment(std::move(comment))
        , seq(std::move(seq))
      {}

      size_t num = 0;
      std::string id;
      std::string comment;
      std::string seq;
    };

    struct MatePairInfo
    {
      MatePairInfo() = default;

      MatePairInfo(   bool seen,
                      int64_t insert_size)
        : seen(seen)
        , insert_size(insert_size)
      {}

      bool seen = false;
      int64_t insert_size = 0;
    };

    struct PairLinkInfo
    {
      PairLinkInfo() {}

      PairLinkInfo(   int64_t gaps,
                      uint64_t links,
                      std::string type)
        : gaps(gaps)
        , links(links)
        , type(type)
      {}

      int64_t gaps = 0;
      uint64_t links = 0;
      std::string type;
    };

    struct KmerInfo
    {
      KmerInfo() {}

      KmerInfo(   std::string tig,
                  uint64_t start,
                  uint64_t end,
                  uint64_t multiple,
                  bool orient)
                  : tig(tig)
                  , start(start)
                  , end(end)
                  , multiple(multiple)
                  , orient(orient)
      {}

      std::string tig = "";
      uint64_t start = 0;
      uint64_t end = 0;
      uint64_t multiple = 1;
      bool orient = 0;
    };
    
    typedef std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo> > mate_pair;
private:
  class Worker
    {
      public:
        void start() { t = std::thread(do_work, this); }
        void join() { t.join(); }

        virtual ~Worker() {}

        Worker& operator=(const Worker& worker) = delete;
        Worker& operator=(Worker&& worker) = delete;

      protected:
        LINKS& links;
        std::thread t;

        Worker(LINKS& links)
          : links(links)
        {}

        Worker(const Worker& worker)
          : Worker(worker.links)
        {}
        Worker(Worker&& worker) noexcept
          : Worker(worker.links)
        {}

        virtual void work() = 0;
        static void do_work(Worker* worker) { worker->work(); }
    };

    class InputWorker : public Worker
    {
      public:
        InputWorker(LINKS& links)
          : Worker(links)
        {}

        InputWorker(const InputWorker& worker)
          : InputWorker(worker.links)
        {}
        InputWorker(InputWorker&& worker) noexcept
          : InputWorker(worker.links)
        {}

        InputWorker& operator=(const InputWorker& worker) = delete;
        InputWorker& operator=(InputWorker&& worker) = delete;

        void work() override;
    };

    class ExtractMatePairWorker : public Worker
    {
      public:
        ExtractMatePairWorker(LINKS& links)
          : Worker(links)
        {}

        ExtractMatePairWorker(const ExtractMatePairWorker& worker)
          : ExtractMatePairWorker(worker.links)
        {}
        ExtractMatePairWorker(ExtractMatePairWorker&& worker) noexcept
          : ExtractMatePairWorker(worker.links)
        {}

        ExtractMatePairWorker& operator=(const ExtractMatePairWorker& worker) = delete;
        ExtractMatePairWorker& operator=(ExtractMatePairWorker&& worker) = delete;

        void work() override;
    };

    class PopulateMateInfoWorker : public Worker
    {
      public:
        PopulateMateInfoWorker(LINKS& links)
          : Worker(links)
        {}

        PopulateMateInfoWorker(const PopulateMateInfoWorker& worker)
          : PopulateMateInfoWorker(worker.links)
        {}
        PopulateMateInfoWorker(PopulateMateInfoWorker&& worker) noexcept
          : PopulateMateInfoWorker(worker.links)
        {}

        PopulateMateInfoWorker& operator=(const PopulateMateInfoWorker& worker) = delete;
        PopulateMateInfoWorker& operator=(PopulateMateInfoWorker&& worker) = delete;

        void work() override;
    };

    static const size_t BUFFER_SIZE = 4;
    static const size_t BLOCK_SIZE = 1;

    std::string seqfile;
    std::queue<std::string> longReads;
    unsigned threads;
    long id;
    size_t buffer_size = 32;
    size_t block_size = 8;

    btllib::KmerBloomFilter* make_bf(uint64_t bfElements, InputParser* linksArgParser);
    void extract_mate_pair(const std::string& seq);
    uint64_t get_file_size(std::string filename);
    bool does_file_exist(std::string fileName);
    std::atomic<bool> fasta{ false };

    btllib::KmerBloomFilter * bloom;
    btllib::SeqReader* reader;
    btllib::OrderQueueSPMC<Read> input_queue;
    InputWorker input_worker;
    std::vector<ExtractMatePairWorker> extract_mate_pair_workers;
    std::vector<PopulateMateInfoWorker> populate_mate_info_workers;
    std::shared_mutex mate_pair_mutex;
    std::mutex mates_mutex;
    

    std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo>> matePair;
    std::unordered_map<uint64_t, KmerInfo> trackAll;
    std::unordered_map<std::string, uint64_t> tigLength;
    // store second mates in a set
    std::unordered_set<uint64_t> mates;
};

inline btllib::KmerBloomFilter*
LINKS::make_bf(uint64_t bfElements, InputParser* linksArgParser){
  btllib::KmerBloomFilter * assemblyBF;
    if(linksArgParser->bfFile != "") {
        std::cout << "A Bloom filter was supplied (" << linksArgParser->bfFile << ") and will be used instead of building a new one from -f " << linksArgParser->assemblyFile << "\n";
        if(!does_file_exist(linksArgParser->bfFile)) {
            std::cout << "\nInvalid file: " << linksArgParser->bfFile <<  " -- fatal\n";
            exit(1);
        } else {
            std::cout << "Checking Bloom filter file " << linksArgParser->bfFile <<"...ok\n";
        }
        // std::cout << "Loading bloom filter of size " << getFileSize(linksArgParser->bfFile) << " from " << linksArgParser->bfFile << "\n";
        assemblyBF = new btllib::KmerBloomFilter(linksArgParser->bfFile);
    } else {
        uint64_t m = ceil((-1 * (double)bfElements * log(linksArgParser->fpr)) / (log(2) * log(2)));
        //uint64_t rem = 64 - (m % 64);
        m = ((uint64_t)(m / 8) + 1) * 8;
        std::cout << "HASHES CALC: " << std::to_string(((double)m / bfElements)) << " second: " << std::to_string(((double)m / bfElements) * log(2)) << "\n";
        unsigned hashFct = floor(((double)m / bfElements) * log(2));
        std::cout << "- Number of bfElements: " << bfElements << "\n"
                    << "- Input file path: " << linksArgParser->bfFile << "\n"
                    << "- Input file: " << linksArgParser->assemblyFile << "\n"
                    << "- kmersize: " << linksArgParser->k << "\n"
                    << "- m: " << m << "\n"
                    << "- fpr: " << linksArgParser->fpr << "\n"
                    << "- hashFct: " << hashFct << "\n";

        std::string reading_tigbloom_message = "\n\n=>Reading contig/sequence assembly file : " + std::to_string(time(0)) + "\n";
        // assemblyruninfo += reading_tigbloom_message;
        // std::cout << "- Filter output file : " << outFileBf << "\n";
        std::cout << "- Filter output file : " << linksArgParser->k << "\n";
        assemblyBF = new btllib::KmerBloomFilter(m/8, hashFct, linksArgParser->k);
        btllib::SeqReader assemblyReader(linksArgParser->assemblyFile, 8, 1);
        // int builder = 0;
        for (btllib::SeqReader::Record record; (record = assemblyReader.read());) {
            // if(builder % 100 == 0) {
            //     std::cout << "reading... builder: " << builder << "\n";
            // }
            // builder++;

            assemblyBF->insert(record.seq);
        }
        std::string bfmsg = "\n\nWriting Bloom filter to disk (" + linksArgParser->bfFile + ") : " + std::to_string(time(0)) + "\n";
        // assemblyruninfo += bfmsg;
        std::cout << bfmsg;
        assemblyBF->save("bftest.out");
        // std::cout << "Done mybf, printing stats...\n";
        
        //printBloomStats(*assemblyBF, std::cout);
        
    }
    return assemblyBF;
}

inline LINKS::LINKS(InputParser* inputParser)
        : inputParser(inputParser)
        , assemblyFile(inputParser->assemblyFile)
        , longReads(inputParser->longReads)
        //, seqfile(longReads.front())
        //, fofFile(inputParser->fofFile)
        , distances(inputParser->distances)
        , k(inputParser->k)
        , verbose(inputParser->verbose)
        , minSize(inputParser->minSize)
        , step(inputParser->step)
        , insertStdev(inputParser->insertStdev)
        , baseName(inputParser->baseName)
        , offset(inputParser->offset)
        , fpr(inputParser->fpr)
        , bfFile(inputParser->bfFile)
        , bfOff(inputParser->bfOff)
        , threads(inputParser->thread)
        , input_queue(buffer_size, block_size)
        , input_worker(*this)
        {}

inline void
LINKS::init_bloom_filter(){
  int64_t bf_elements = get_file_size(inputParser->assemblyFile);
  bloom = make_bf(bf_elements,inputParser);
}



inline void
LINKS::InputWorker::work()
{
  if (links.reader->get_format() == btllib::SeqReader::Format::FASTA) {
    links.fasta = true;
  } else {
    links.fasta = false;
  }

  decltype(links.input_queue)::Block block(links.block_size);
  size_t current_block_num = 0;
  Read read;

  uint counterx = 0;
  for(auto record : (*(links.reader))){
    block.data[block.count++] = Read(record.num,
                                    std::move(record.id),
                                    std::move(record.comment),
                                    std::move(record.seq));
    if (block.count == links.block_size) {
      block.num = current_block_num++;
      links.input_queue.write(block);
      block.count = 0;
    }
  }

  if (block.count > 0) {
    block.num = current_block_num++;
    links.input_queue.write(block);
  }
  for (unsigned i = 0; i < links.threads; i++) {
    block.num = current_block_num++;
    block.current = 0;
    block.count = 0;
    links.input_queue.write(block);
  }
}

inline void 
LINKS::start_read_fasta(){
  matePair.reserve(10000);
  mates.reserve(10000);
  
  reader = new btllib::SeqReader(std::move(longReads.front()), btllib::SeqReader::Flag::LONG_MODE);
  longReads.pop();
  extract_mate_pair_workers = std::vector<ExtractMatePairWorker>(threads, ExtractMatePairWorker(*this));
  input_worker.start();
  // start
  for (auto& worker : extract_mate_pair_workers) {
    worker.start();
  }
  // wait
  for (auto& worker : extract_mate_pair_workers) {
    worker.join();
  }
  input_worker.join();

  std::cout << "matepair size: " << matePair.size() << std::endl;
  std::cout << "mates size: " << mates.size() << std::endl;
}
inline void 
LINKS::start_read_contig(){
  trackAll.reserve(mates.size());

  delete(reader);
  reader = new btllib::SeqReader(assemblyFile, btllib::SeqReader::Flag::LONG_MODE);

  populate_mate_info_workers = std::vector<PopulateMateInfoWorker>(threads, PopulateMateInfoWorker(*this));

  for (auto& worker : populate_mate_info_workers) {
    worker.start();
  }
  // wait
  for (auto& worker : populate_mate_info_workers) {
    worker.join();
  }
}

inline void
LINKS::extract_mate_pair(const std::string& seq)
{

  for(auto distance : distances){
    //std::cout << "dist: " << distance << std::endl;
    uint64_t delta = distance - k;
    //uint64_t delta = distance;
    int breakFlag = 0;
    bool reverseExists = false;
    btllib::NtHash nthash(seq, bloom->get_hash_num(), k);
    btllib::NtHash nthashLead(seq, bloom->get_hash_num(), k, delta);

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
        mate_pair_mutex.lock_shared();
        mate_pair::iterator it = matePair.find(nthashLead.get_reverse_hash());
        mate_pair_mutex.unlock_shared();
        if( it != matePair.end() ) {
            std::unordered_map<uint64_t, MatePairInfo > &innerMap = it->second;
            mate_pair_mutex.lock_shared();
            std::unordered_map<uint64_t, MatePairInfo>::iterator innerit = innerMap.find(nthash.get_reverse_hash());
            mate_pair_mutex.unlock_shared();
            if( innerit != innerMap.end() ){
                mate_pair_mutex.lock();
                innerit->second.insert_size = distance;
                mate_pair_mutex.unlock();
                reverseExists = true;
            }
        }
        
        if(!reverseExists && bloom->contains(nthash.hashes()) && bloom->contains(nthashLead.hashes())) { // May need to change with forward reverse hashes
            mate_pair_mutex.lock_shared();
            mate_pair::iterator it = matePair.find(nthash.get_forward_hash());
            mate_pair_mutex.unlock_shared();
            if( it != matePair.end() ) {
                std::unordered_map<uint64_t, MatePairInfo> &innerMap = it->second;
                mate_pair_mutex.lock_shared();
                std::unordered_map<uint64_t, MatePairInfo>::iterator innerit = innerMap.find(nthashLead.get_forward_hash());
                mate_pair_mutex.unlock_shared();
                if( innerit != innerMap.end() ){
                    mate_pair_mutex.lock();
                    matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()].insert_size = distance;
                    mate_pair_mutex.unlock();
                }else{
                    mate_pair_mutex.lock();
                    matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()] = MatePairInfo(false, distance);
                    mate_pair_mutex.unlock();
                }
            }else{
                mate_pair_mutex.lock();
                matePair[nthash.get_forward_hash()][nthashLead.get_forward_hash()] = MatePairInfo(false, distance);
                mate_pair_mutex.unlock();
            }

            std::lock_guard<std::mutex> guard(mates_mutex);
            mates.insert(nthashLead.get_forward_hash());
        }
    }
  }
} 

inline void
LINKS::ExtractMatePairWorker::work()
{
  decltype(links.input_queue)::Block input_block(links.block_size);

  auto start = std::chrono::high_resolution_clock::now();
  for (;;) {
    if (input_block.current == input_block.count) {
      links.input_queue.read(input_block);
    }
    
    if (input_block.count == 0) {
      break;
    }
    
    Read& read = input_block.data[input_block.current++];

    if (links.k <= read.seq.size()) {
        links.extract_mate_pair(read.seq);
    } else {
      continue; // nothing
    }

  }
  auto finish = std::chrono::high_resolution_clock::now();
  auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
  std::cout << microseconds.count() << "µs\n";
}

inline void
LINKS::PopulateMateInfoWorker::work()
{
  decltype(links.input_queue)::Block input_block(links.block_size);

  auto start = std::chrono::high_resolution_clock::now();
  for (;;) {
    if (input_block.current == input_block.count) {
      links.input_queue.read(input_block);
    }
    
    if (input_block.count == 0) {
      break;
    }
    Read& read = input_block.data[input_block.current++];

    if (links.k <= read.seq.size()) {
        //links.extract_mate_pair(read.seq);
        continue;
    } else {
      continue; // nothing
    }

  }
  auto finish = std::chrono::high_resolution_clock::now();
  auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
  std::cout << microseconds.count() << "µs\n";
}

inline LINKS::~LINKS()
{
  
}

inline uint64_t 
LINKS::get_file_size(std::string filename)
{
    // This buffer is a stat struct that the information is placed concerning the file.
    struct stat stat_buf;
    // Stat method returns true if successfully completed
    int rc = stat(filename.c_str(), &stat_buf);
    // st_size holds the total size of the file in bytes
    return rc == 0 ? stat_buf.st_size : -1;
}

inline bool 
LINKS::does_file_exist(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}
#endif