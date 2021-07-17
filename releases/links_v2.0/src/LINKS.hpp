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



class LINKS
{
  public:

    LINKS(InputParser* inputParser);
    void init_bloom_filter();
    void start_read_fasta();
    void start_read_contig();
    void pair_contigs();
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
      /*
      MatePairInfo(MatePairInfo& copy)
        : seen(copy.seen)
        , insert_size(copy.insert_size)
      {}
*/
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
                      std::string type = "")
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
    size_t buffer_size = 4;
    size_t block_size = 1;

    btllib::KmerBloomFilter* make_bf(uint64_t bfElements, InputParser* linksArgParser);
    void extract_mate_pair(const std::string& seq, std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo>>& own_mate_pair,std::unordered_set<uint64_t>& own_mates);
    void populate_mate_info(const std::string& seq,const std::string contig_rank);
    
    int getDistanceBin(int distance);
    int getDistance(uint64_t insert_size, 
      uint64_t length_i, 
      uint64_t start_i, 
      uint64_t start_j);
    void addToPairMap(
      int isz,
      std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo>>>& pair,
      int distance,
      std::string kmer1_name,
      std::string kmer2_name,
      unsigned orient_enum);
    void merge_mate_pair_map(mate_pair& own_mate_pair);
    void merge_mates_set(std::unordered_set<uint64_t>& own_mates);
    // helper functions
    uint64_t get_file_size(std::string filename);
    bool does_file_exist(std::string fileName);
    std::atomic<bool> fasta{ false };

    //btllib::KmerBloomFilter * bloom;
    std::shared_ptr<btllib::KmerBloomFilter> bloom;
    //btllib::SeqReader* reader;
    std::shared_ptr<btllib::SeqReader> reader;
    std::shared_ptr<btllib::OrderQueueSPMC<Read>> input_queue;
    //InputWorker* input_worker;
    std::shared_ptr<InputWorker> input_worker;
    std::vector<ExtractMatePairWorker> extract_mate_pair_workers;
    std::vector<PopulateMateInfoWorker> populate_mate_info_workers;

    // commented out muro debug
    //std::shared_mutex mate_pair_mutex;
    //std::mutex mates_mutex;
    std::mutex mate_pair_mutex;
    std::mutex contig_rank_mutex;
    std::mutex tig_length_mutex; 
    std::shared_mutex track_all_mutex;
    

    std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo>> matePair;
    std::unordered_map<uint64_t, KmerInfo> trackAll;
    std::unordered_map<std::string, uint64_t> tigLength;
    // store second mates in a set
    std::unordered_set<uint64_t> mates;
    unsigned last_contig_rank = 0;
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
        uint builder = 0;
        for (btllib::SeqReader::Record record; (record = assemblyReader.read());) {
          if(builder % 10 == 0) {
                //std::cout << "reading... builder: " << builder << "\n";
          }
          builder++;
          //std::cout << "here 1\n";
          assemblyBF->insert(record.seq);
          //std::cout << "here 2\n";
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
        , input_queue(std::shared_ptr<btllib::OrderQueueSPMC<Read>>(new btllib::OrderQueueSPMC<Read>(buffer_size, block_size)))
        , input_worker(std::shared_ptr<InputWorker>(new InputWorker(*this)))
        {}

inline void
LINKS::init_bloom_filter(){
  int64_t bf_elements = get_file_size(inputParser->assemblyFile);
  bloom = std::shared_ptr<btllib::KmerBloomFilter>(make_bf(bf_elements,inputParser));
}



inline void
LINKS::InputWorker::work()
{
  //std::cout << "test 0\n";
  if (links.reader->get_format() == btllib::SeqReader::Format::FASTA) {
    links.fasta = true;
  } else {
    links.fasta = false;
  }
  //std::cout << "test 1\n";
  //decltype(*(links.input_queue))::Block block(links.block_size);
  //(std::remove_reference<decltype(*(links.input_queue))>)::Block input_block(links.block_size);
  btllib::OrderQueueSPMC<Read>::Block block(links.block_size); 

  size_t current_block_num = 0;
  size_t read_c = 0;
  Read read;
  //std::cout << "test 2\n";

  uint counterx = 0;
  for(auto record : (*(links.reader))){
    read_c++;
    if(read_c % 100 == 0){
      std::cout << "read_counter: " << read_c << std::endl;
    }

    block.data[block.count++] = Read(record.num,
                                    std::move(record.id),
                                    std::move(record.comment),
                                    std::move(record.seq));
    if (block.count == links.block_size) {
      block.num = current_block_num++;
      links.input_queue->write(block);
      block.count = 0;
    }
    //std::cout << "test 2.5\n";
  }
  //std::cout << "test 3\n";

  if (block.count > 0) {
    block.num = current_block_num++;
    links.input_queue->write(block);
  }
  for (unsigned i = 0; i < links.threads; i++) {
    block.num = current_block_num++;
    block.current = 0;
    block.count = 0;
    links.input_queue->write(block);
  }
}

inline void 
LINKS::start_read_fasta(){
  //matePair.reserve(10000);
  //mates.reserve(10000);
  
  reader = std::shared_ptr<btllib::SeqReader>(new btllib::SeqReader(std::move(longReads.front()), btllib::SeqReader::Flag::LONG_MODE));
  longReads.pop();
  extract_mate_pair_workers = std::vector<ExtractMatePairWorker>(threads, ExtractMatePairWorker(*this));
  input_worker->start();
  // start
  for (auto& worker : extract_mate_pair_workers) {
    worker.start();
  }
  // wait
  for (auto& worker : extract_mate_pair_workers) {
    worker.join();
    std::cout << "joined\n";
  }
  input_worker->join();

  std::cout << "after join matepair size: " << matePair.size() << std::endl;
  uint second_dim_size = 0;
  uint third_dim_size = 0;
  uint fourth_dim_size = 0;
  for (auto it = matePair.begin(); it != matePair.end(); it++) {
    //std::cout << *it << endl;
    //std::cout << "it->second.size(): " << it->second.size() << std::endl; 
    second_dim_size += it->second.size();
    for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      third_dim_size += it2->second.insert_size;
    }
  }
  std::cout << "second dim size: " << second_dim_size << std::endl;
  std::cout << "total insert size size: " << third_dim_size << std::endl; 
  std::cout << "mates size: " << mates.size() << std::endl;
  sleep(5);
}
inline void 
LINKS::start_read_contig(){
  trackAll.reserve(mates.size());
  //std::cout << "t test 0\n";
  //delete(reader);
  std::cout << "t test 0\n";
  reader.reset();
  reader = std::shared_ptr<btllib::SeqReader>(new btllib::SeqReader(assemblyFile, btllib::SeqReader::Flag::LONG_MODE));
  //delete input_worker;
  std::cout << "t test 1\n";
  input_worker.reset();
  std::cout << "t test 2\n";

  input_queue.reset();
  std::cout << "t test 3\n";
  input_queue = std::shared_ptr<btllib::OrderQueueSPMC<Read>>(new btllib::OrderQueueSPMC<Read>(buffer_size, block_size));
  //std::cout << "t test 1.2\n";

  std::cout << "t test 4\n";

  input_worker = std::shared_ptr<InputWorker>(new InputWorker(*this));
  input_worker->start();
  //std::cout << "t test 2\n";
  std::cout << "t test 5\n";

  populate_mate_info_workers = std::vector<PopulateMateInfoWorker>(threads, PopulateMateInfoWorker(*this));

  std::cout << "t test 6\n";  
  for (auto& worker : populate_mate_info_workers) {
    worker.start();
  }
  //std::cout << "t test 4\n";
  // wait
  std::cout << "t test 7\n";
  for (auto& worker : populate_mate_info_workers) {
    worker.join();
    std::cout << "t test 8\n";
  }
  //std::cout << "t test 5\n";
  input_worker->join();
  //std::cout << "t test 6\n";
  input_queue.reset();
  //std::cout << "t test 7\n";
  
  std::cout << "trackAll.size() " << trackAll.size() << std::endl;
}

inline void
LINKS::extract_mate_pair(const std::string& seq,
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo>>& own_mate_pair,
  std::unordered_set<uint64_t>& own_mates)
{
  //std::cout << "here 1" << std::endl;
  for(auto distance : distances){
    //std::cout << "dist: " << distance << std::endl;
    uint64_t delta = distance - k;
    //uint64_t delta = distance;
    int breakFlag = 0;
    bool reverseExists = false;
    btllib::NtHash nthash(seq, bloom->get_hash_num(), k);
    btllib::NtHash nthashLead(seq, bloom->get_hash_num(), k, delta);

    uint64_t hash_a, hash_b;    

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

        if(nthash.get_forward_hash() < nthashLead.get_reverse_hash()){
          hash_a = nthashLead.get_reverse_hash();
          hash_b = nthash.get_reverse_hash();
        }else{
          hash_a = nthash.get_forward_hash();
          hash_b = nthashLead.get_forward_hash();
        }

        /*
        mate_pair::iterator it = own_mate_pair.find(nthashLead.get_reverse_hash());

        if( it != own_mate_pair.end() ) {
            std::unordered_map<uint64_t, MatePairInfo > &innerMap = it->second;

            std::unordered_map<uint64_t, MatePairInfo>::iterator innerit = innerMap.find(nthash.get_reverse_hash());

            if( innerit != innerMap.end() ){

                innerit->second.insert_size = distance;

                reverseExists = true;
            }
        }
        */
        if(bloom->contains(nthash.hashes()) && bloom->contains(nthashLead.hashes())) { // May need to change with forward reverse hashes
            mate_pair::iterator it = own_mate_pair.find(hash_a);
            if( it != own_mate_pair.end() ) {
                std::unordered_map<uint64_t, MatePairInfo> &innerMap = it->second;
                std::unordered_map<uint64_t, MatePairInfo>::iterator innerit = innerMap.find(hash_b);
                if( innerit != innerMap.end() ){
                    own_mate_pair[hash_a][hash_b].insert_size = distance;
                }else{
                    own_mate_pair[hash_a][hash_b] = MatePairInfo(false, distance);
                }
            }else{
                own_mate_pair[hash_a][hash_b] = MatePairInfo(false, distance);
            }
            //std::lock_guard<std::mutex> guard(mates_mutex);
            own_mates.insert(hash_b);
        }
    }
  }
} 

inline void
LINKS::populate_mate_info(const std::string& seq,
                          const std::string contig_rank){
  
  btllib::NtHash ntHashContig(seq, bloom->get_hash_num(), k); // hashFunc can be 1 after first step
  //std::cout << "t test 20\n";
  //std::cout << "matepair size:" << matePair.size() << " t test 21\n";

  uint second_dim_size = 0;
    uint third_dim_size = 0;
    uint fourth_dim_size = 0;
    for (auto it = matePair.begin(); it != matePair.end(); it++) {
    //std::cout << *it << endl;
    //std::cout << "it->second.size(): " << it->second.size() << std::endl; 
    second_dim_size += it->second.size();
    for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      third_dim_size += it2->second.insert_size;
    }
  }
  //std::cout << "second dim size: " << second_dim_size << std::endl;
  //std::cout << "total insert size size: " << third_dim_size << std::endl; 

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
      //std::cout << "t test 21\n";
        // Forward part
	    if(matePair.find(ntHashContig.get_forward_hash()) != matePair.end() || 
            mates.find(ntHashContig.get_forward_hash()) != mates.end()) {
            track_all_mutex.lock_shared();
            
            //std::cout << "matepair size:" << matePair.size() << " test xx\n";
            //std::cout << "t test yy\n";
            if(trackAll.find(ntHashContig.get_forward_hash()) == trackAll.end()) {
              //std::cout << "t test zz\n";
              track_all_mutex.unlock_shared();
              track_all_mutex.lock();
              trackAll[ntHashContig.get_forward_hash()] = KmerInfo(contig_rank, i, i + k, 1, false);
              //std::cout << "trackAll not found\n";
              track_all_mutex.unlock();
            } else {
              track_all_mutex.unlock_shared();
              track_all_mutex.lock();
              trackAll[ntHashContig.get_forward_hash()].multiple +=1;
              track_all_mutex.unlock();
            }
        }

        // Reverse part
        if(matePair.find(ntHashContig.get_reverse_hash()) != matePair.end() || 
            mates.find(ntHashContig.get_reverse_hash()) != mates.end()) {
            track_all_mutex.lock_shared();
            if(trackAll.find(ntHashContig.get_reverse_hash()) == trackAll.end()) {
              track_all_mutex.unlock_shared();
              track_all_mutex.lock();
              trackAll[ntHashContig.get_reverse_hash()] = KmerInfo(contig_rank, i, i + k, 1, true);
              track_all_mutex.unlock();
            } else {
              track_all_mutex.unlock_shared();
              track_all_mutex.lock();
              trackAll[ntHashContig.get_reverse_hash()].multiple +=1;
              track_all_mutex.unlock();
            }
        }
    }
}
inline void 
LINKS::merge_mates_set(std::unordered_set<uint64_t>& own_mates){
  //main_mates.insert(own_mates.begin(), own_mates.end());
}
inline void
LINKS::merge_mate_pair_map(mate_pair& own_mate_pair){

  //std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo> > 

  for (auto own_first_mate = own_mate_pair.begin(); 
    own_first_mate != own_mate_pair.end(); own_first_mate++) {
    //std::cout << "own_first_mate: " << own_first_mate->first << std::endl;
    auto main_first_mate = matePair.find(own_first_mate->first);
    if(main_first_mate != matePair.end()){
      //std::cout << "yes-found" << std::endl;
      for (auto own_second_mate = own_first_mate->second.begin(); 
        own_second_mate != own_first_mate->second.end(); own_second_mate++) {

        auto main_second_mate = main_first_mate->second.find(own_second_mate->first);
        if(main_second_mate == main_first_mate->second.end()){
          matePair[own_first_mate->first][own_second_mate->first] = MatePairInfo(false, own_second_mate->second.insert_size);
        }
      }
    } else{
      //std::cout << "no-found" << std::endl;
      for (auto own_second_mate = own_first_mate->second.begin(); 
        own_second_mate != own_first_mate->second.end(); own_second_mate++) {
          matePair[own_first_mate->first][own_second_mate->first] = MatePairInfo(false, own_second_mate->second.insert_size);
      }
    }
  }
} 

inline void
LINKS::ExtractMatePairWorker::work()
{
  //std::cout << "y test 1\n";
  //decltype(*(links.input_queue))::Block input_block(links.block_size);
  //(std::remove_reference<decltype(*(links.input_queue))>)::Block input_block(links.block_size);
  btllib::OrderQueueSPMC<Read>::Block input_block(links.block_size); 

  std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo>> own_mate_pair;
  std::unordered_set<uint64_t> own_mates;

  //std::cout << "y test 2\n";
  for (;;) {
    if (input_block.current == input_block.count) {
      links.input_queue->read(input_block);
    }
    
    if (input_block.count == 0) {
      break;
    }
    
    Read& read = input_block.data[input_block.current++];

    if (links.k <= read.seq.size()) {
        links.extract_mate_pair(read.seq,own_mate_pair,own_mates);
    } else {
      continue; // nothing
    }
  }
    links.mate_pair_mutex.lock();
    links.merge_mate_pair_map(own_mate_pair);
    own_mate_pair.clear();
    links.merge_mates_set(own_mates);
    std::cout << "main mate pair first dim size: " << links.matePair.size() << std::endl;
    uint second_dim_size = 0;
    uint third_dim_size = 0;
    uint fourth_dim_size = 0;
    for (auto it = links.matePair.begin(); it != links.matePair.end(); it++) {
    //std::cout << *it << endl;
    //std::cout << "it->second.size(): " << it->second.size() << std::endl; 
    second_dim_size += it->second.size();
    for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      third_dim_size += it2->second.insert_size;

    }
  }
    std::cout << "second dim size: " << second_dim_size << std::endl;
    std::cout << "total insert size size: " << third_dim_size << std::endl; 
    links.mate_pair_mutex.unlock();
  //std::cout << "y test 3\n";
}

inline void
LINKS::PopulateMateInfoWorker::work()
{
  //std::cout << "q test 1\n";
  //decltype(*(links.input_queue))::Block input_block(links.block_size);
  //(std::remove_reference<decltype(*(links.input_queue))>)::Block input_block(links.block_size);
  btllib::OrderQueueSPMC<Read>::Block input_block(links.block_size);
  unsigned cur_contig_rank;

  std::cout << "t test 10\n";
  //std::cout << "q test 2\n";
  for (;;) {
    //std::cout << "q test 3\n";
    if (input_block.current == input_block.count) {
      //std::cout << "q test 3.5\n";
      links.input_queue->read(input_block);
    }
    //std::cout << "q test 4\n";
    if (input_block.count == 0) {
      break;
    }
    Read& read = input_block.data[input_block.current++];
    //std::cout << "q test 5\n";

    links.contig_rank_mutex.lock();
    cur_contig_rank = links.last_contig_rank;
    ++links.last_contig_rank;
    links.contig_rank_mutex.unlock();

    if (links.k <= read.seq.size()) {
        //std::cout << "t test 11\n";
        //std::cout << "read size: " << read.seq.size() << std::endl;
        //std::cout << "std::to_string(cur_contig_rank) " << std::to_string(cur_contig_rank) << std::endl;
        //std::cout << "tigLength.size() " << links.tigLength.size() << std::endl;
        links.tig_length_mutex.lock(); 
        links.tigLength[std::to_string(cur_contig_rank)] = read.seq.size();
        links.tig_length_mutex.unlock();
        //std::cout << "t test 12\n";
        links.populate_mate_info(read.seq,std::to_string(cur_contig_rank));
        //std::cout << "q test 6\n";
        continue;
    } else {
      //std::cout << "q test 7\n";
      continue; // nothing
    }
  }
}

inline LINKS::~LINKS()
{
 std::cout << "here deleting LINKS" << std::endl;
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

inline void 
LINKS::addToPairMap(
    int isz,
    std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo>>>& pair,
    int distance,
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

        //std::cout << "orient enum: " << orient_enum << std::endl;

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
        //std::cout << "!!!!!!!" << std::endl;
        //std::cout << std::get<0>(first_pair) << " " << std::get<1>(first_pair) << std::endl;
        //std::cout << std::get<0>(second_pair) << " " << std::get<1>(second_pair) << std::endl;

        if(pair.find(std::get<0>(first_pair)) == pair.end() 
            || pair[std::get<0>(first_pair)].find(isz) == pair[std::get<0>(first_pair)].end() 
            || pair[std::get<0>(first_pair)][isz].find(std::get<1>(first_pair)) == pair[std::get<0>(first_pair)][isz].end()) {
            // std::cout << "Checkpoint 7.1 adding to pair new GAPSLINKS\n";
            pair[std::get<0>(first_pair)][isz][std::get<1>(first_pair)] = PairLinkInfo(distance,1);
        } else {
            // std::cout << "Checkpoint 7.2 adding to pair existing gapslings\n";
            pair[std::get<0>(first_pair)][isz][std::get<1>(first_pair)].gaps += distance;
            pair[std::get<0>(first_pair)][isz][std::get<1>(first_pair)].links += 1;
        }
        if(pair.find(std::get<0>(second_pair)) == pair.end() 
            || pair[std::get<0>(second_pair)].find(isz) == pair[std::get<0>(second_pair)].end() 
            || pair[std::get<0>(second_pair)][isz].find(std::get<1>(second_pair)) == pair[std::get<0>(second_pair)][isz].end()) {
            // std::cout << "Checkpoint 7.3 adding to pair new GAPSLINKSs\n";
            pair[std::get<0>(second_pair)][isz][std::get<1>(second_pair)] = PairLinkInfo(distance,1);
        } else {
            pair[std::get<0>(second_pair)][isz][std::get<1>(second_pair)].gaps += distance;
            pair[std::get<0>(second_pair)][isz][std::get<1>(second_pair)].links += 1;
        }
}

inline int 
LINKS::getDistanceBin(int distance)
{
    return distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 10000;
}

inline int 
LINKS::getDistance(uint64_t insert_size, uint64_t length_i, uint64_t start_i, uint64_t start_j) {

   int insert_span = (length_i - start_i) + start_j;
   int gap_or_overlap = insert_size - insert_span;

   return gap_or_overlap;
}
inline void
LINKS::pair_contigs()
{
  /*
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
*/
  std::string issues = inputParser->baseName + ".pairing_issues";
  std::string distribution = inputParser->baseName + ".pairing_distribution.csv";
  std::string tigpair_checkpoint = inputParser->baseName + ".tigpair_checkpoint.tsv"; // add a checkpoint file, prevent re-running LINKS from scratch if crash
  std::string simplepair_checkpoint = inputParser->baseName + ".simplepair_checkpoint.tsv"; // add a checkpoint file, prevent re-running LINKS from scratch if cras
  uint64_t totalPairs = 0;
  uint64_t ct_illogical = 0, ct_ok_contig = 0, ct_ok_pairs = 0, ct_problem_pairs = 0, ct_iz_issues = 0, ct_single = 0, ct_multiple = 0, ct_both = 0, trackInsert = 0;
  std::unordered_map<uint64_t, uint64_t> ct_single_hash, ct_both_hash, ct_illogical_hash, ct_ok_contig_hash, ct_ok_pairs_hash, ct_problem_pairs_hash, ct_iz_issues_hash;
  // Mapping of tiga_head -> insertSize -> tigb_head -> links & gaps
  std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo> > > pair;
  std::unordered_map<std::string, std::unordered_map<std::string, PairLinkInfo> >simplePair;
  std::unordered_map<std::string, PairLinkInfo> err;
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
  std::unordered_map<uint64_t, MatePairInfo>::iterator mateListItr;
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo> >::iterator matePairItr;
  
  uint counter_1 = 0, counter_2 = 0, counter_3 = 0, counter_4 = 0, counter_5 = 0, counter_6 = 0, counter_7 = 0;  
  for(matePairItr = matePair.begin(); matePairItr != matePair.end(); matePairItr++) {
        for(mateListItr = matePairItr->second.begin(); mateListItr != matePairItr->second.end(); mateListItr++) {
        ++counter_1;    
        //std::cout << "trackAll[matePairItr->first].multiple " << trackAll[matePairItr->first].multiple << std::endl;
        //std::cout << "trackAll[mateListItr->first].multiple " << trackAll[mateListItr->first].multiple << std::endl;
	if( mateListItr->second.seen == false &&                 //matepair is not seen
                trackAll.find(matePairItr->first) != trackAll.end() &&  //first mate is tracked
                trackAll[matePairItr->first].multiple == 1 &&      //first mate seen once
                trackAll.find(mateListItr->first) != trackAll.end() &&  //second mate is tracked
                trackAll[mateListItr->first].multiple == 1) {      //second mate is seen once
                ++counter_2;
		mateListItr->second.seen = true;
                insert_size = matePair[matePairItr->first][mateListItr->first].insert_size;
                min_allowed = -1 * (insertStdev * insert_size); // check int
                low_iz = insert_size + min_allowed; // check int
                up_iz = insert_size - min_allowed; // check int
                if(verbose) std::cout << "Pair read1Hash=" << matePairItr->first << " read2Hash=" << mateListItr->first << "\n";
                if(trackAll[matePairItr->first].tig != "" && trackAll[mateListItr->first].tig != "") { //double check if tig names not null
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

                    //std::cout << "pre kmer1: " << kmer1.tig << " kmer2: " << kmer2.tig << std::endl;

                    if(kmer1.tig != kmer2.tig) { // paired reads located on <> contigs
                    //std::cout << "aft kmer1: " << kmer1.tig << " kmer2: " << kmer2.tig << std::endl;
                    //std::cout << "kmer1 ori: " << kmer1.orient << " kmer2 ori: " << kmer2.orient << std::endl;
                    //std::cout << "tigLength[kmer1.tig] " << tigLength[kmer1.tig] 
                    //<< " kmer1.start " << kmer1.start << std::endl;
                    //std::cout << "tigLength[kmer2.tig] " << tigLength[kmer2.tig] 
                    //<< " kmer2.start " << kmer2.start << std::endl;
                        // MURATHAN DEBUG 11.5.21
                        if(!kmer1.orient){             // if kmer1 is forward
                            if(!kmer2.orient){         // if kmer2 is forward
                                ++counter_4;
                                //std::cout << "here\n";
        //std::cout << "tigLength[kmer1.tig] " << tigLength[kmer1.tig] 
        //<< " kmer1.start " << kmer1.start << " kmer2.start " << kmer2.start << std::endl;
				distance = getDistance(insert_size, tigLength[kmer1.tig], kmer1.start, kmer2.start);
        //std::cout << "distance: " << distance << std::endl;
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.tig, kmer2.tig, 0);
                                }
                            }else{ 
				++counter_5;  
                             distance = getDistance(insert_size, tigLength[kmer1.tig], kmer1.start, tigLength[kmer2.tig] - kmer2.end);
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.tig, kmer2.tig, 1);
                                }
                            }
                        }else{                              // if kmer1 is reverse
				if(!kmer2.orient){ 
				++counter_6;
                                distance = getDistance(insert_size, tigLength[kmer1.tig], tigLength[kmer1.tig] - kmer1.end, kmer2.start);
                                if(distance > min_allowed && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.tig, kmer2.tig, 2);
                                }
                            }else{                          // if kmer2 is reverse
                                ++counter_7;
				distance = getDistance(insert_size, tigLength[kmer2.tig], kmer2.end, kmer1.end);
                                if(distance > min_allowed  && distance < insert_size){
                                    isz = getDistanceBin(distance);
                                    addToPairMap( isz, pair, distance, kmer1.tig, kmer2.tig, 3);
                                }
                            }
                        }
                        //std::cout << "distance: " << distance << std::endl;
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
                if(matePair[matePairItr->first][mateListItr->first].seen == false) {
                    //std::cout << "UNSEEN getBT() increment ct\n";
                    ct_multiple++;
                }

            }
        } // pairing read b
    } // pairing read a

/*
	std::cout << "counter_1: " << counter_1 << std::endl;
	std::cout << "counter_2: " << counter_2 << std::endl;
	std::cout << "counter_3: " << counter_3 << std::endl;
	std::cout << "counter_4: " << counter_4 << std::endl;
	std::cout << "counter_5: " << counter_5 << std::endl;
	std::cout << "counter_6: " << counter_6 << std::endl;
	std::cout << "counter_7: " << counter_7 << std::endl;
	std::cout << "pair size: " << pair.size() << std::endl;
*/

  // std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo>>
  // (12, (15, (false, 1000)))
  uint second_dim_size = 0;
  uint third_dim_size = 0;
  uint fourth_dim_size = 0;
  for (auto it = matePair.begin(); it != matePair.end(); it++) {
    //std::cout << *it << endl;
    //std::cout << "it->second.size(): " << it->second.size() << std::endl; 
    second_dim_size += it->second.size();
    for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      //third_dim_size += it2->second.links;
      /*for (auto it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
        fourth_dim_size += it3->second.size(); 
      }*/
    }
  } 
  /*
  std::cout << "second dim size: " << second_dim_size << std::endl;
  std::cout << "third dim size: " << third_dim_size << std::endl;
  std::cout << "fourth dim size: " << fourth_dim_size << std::endl;
  */

  // std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo> > >:
  //("f13",(1000, ("r12",(12000,12))))
  uint pair_second_dim_size = 0;
  uint pair_third_dim_size = 0;
  uint pair_fourth_dim_size = 0;
  for (auto it = pair.begin(); it != pair.end(); it++) {
    //std::cout << *it << endl;
    //std::cout << "it->second.size(): " << it->second.size() << std::endl;
    //std::cout << "*********************" << std::endl; 
    pair_second_dim_size += it->second.size();
    for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      pair_third_dim_size += it2->second.size();
      //std::cout << "---------------------" << std::endl;
      for (auto it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
        //std::cout << "map: " << it3->first << std::endl;
        pair_fourth_dim_size += it3->second.links;
      }
    }
  } 
  
  std::cout << "second dim size: " << pair_second_dim_size << std::endl;
  std::cout << "third dim size: " << pair_third_dim_size << std::endl;
  std::cout << "fourth dim size: " << pair_fourth_dim_size << std::endl;
  

	//uint second_dim_size = 0;
	//std;;unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo>>>::iterator pair_iterator;
	//for(itr = pair_iterator.begin(); itr != pair_iterator.end(); itr++){
	//	second_dim_size += itr->second.size();
	
	//}
	//std::cout << "second dim size: " << second_dim_size << std::endl;

    // Summary of the contig pair issues
    //std::cout << "------------- Putative issues with contig pairing - Summary  ----------------\n";
    //std::cout << "err map size: " << err.size() << "\n"; 
    
    //sortErr(err); TODO: uncomment
    std::unordered_map<std::string, PairLinkInfo>::iterator errItr;
    for(errItr = err.begin(); errItr != err.end(); errItr++) {
        double mean_iz = 0;
        if(errItr->second.links) {
            mean_iz = errItr->second.gaps / errItr->second.links;
        }
        // std::cout << "Pair " << errItr->first << " has " << errItr->second.getLinks() << " Links and mean distance of = " << mean_iz << "\n";
    }

    uint64_t satisfied = ct_ok_pairs + ct_ok_contig;
    uint64_t unsatisfied = ct_problem_pairs + ct_iz_issues + ct_illogical;
    uint64_t ct_both_reads = ct_both * 2;

    
    //std::cout << "THESE ARE THE FILTERINGS:\n"<< "filter 1: "<< std::to_string(filter1) << "\n" << "filter 2: "<< std::to_string(filter2) << "\n" << "filter 3: "<< std::to_string(filter3) << "\n" << "filter 4: "<< std::to_string(filter4) << "\n";
    //std::cout << "THESE ARE THE COUNTERS:\n" << std::to_string(CheckCounterBase) << "\n0 " << std::to_string(Check0Counter) << "\n1 " <<  std::to_string(Check1Counter) << "\n2 " <<  std::to_string(Check2Counter) << "\n3 " <<  std::to_string(Check3Counter) << "\n4 " <<  std::to_string(Check4Counter) << "\n5 " <<  std::to_string(Check5Counter) << "\n6 " <<  std::to_string(Check6Counter) << "\n7 " <<  std::to_string(Check7Counter) << "\n8 " <<  std::to_string(Check8Counter) << "\n9 " <<  std::to_string(Check9Counter) << "\n10 " << std::to_string(Check10Counter) << "\n11 " << std::to_string(Check11Counter) << "\n12 " << std::to_string(Check12Counter) << "\n13 " << std::to_string(Check13Counter) << "\n14 " << std::to_string(Check14Counter) << "\n15 " << std::to_string(Check15Counter) << "\n16 " << std::to_string(Check16Counter) << "\n17 " << std::to_string(Check17Counter) << "\n18 " << std::to_string(Check18Counter) << "\n19" << std::to_string(Check19Counter) << "\n20 " << std::to_string(Check20Counter) << "\n21 " << std::to_string(Check21Counter) << "\n22 " << std::to_string(Check22Counter) << "\n23 " << std::to_string(Check23Counter) << "\n24 " << std::to_string(Check24Counter) << "\n25 " << std::to_string(Check25Counter) << "\n26 " << std::to_string(Check26Counter) << "\n";
    std::cout << "\n===========PAIRED K-MER STATS===========\n";
    //std::cout << "Total number of pairs extracted from -s " << longReadsFile << " " << totalPairs << "\n";
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
    //sortInsertSize(ct_both_hash); TODO: uncomment
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
    std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo> > >::iterator pairItr;
    std::cout << "size of TIGPAIR: " << pair.size() << "\n";
    for(pairItr = pair.begin(); pairItr != pair.end(); pairItr++) {
        std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo> >::iterator insertSizes;
        for(insertSizes = pairItr->second.begin(); insertSizes != pairItr->second.end(); insertSizes++) {
            std::unordered_map<std::string, PairLinkInfo>::iterator secPairItr;
            for(secPairItr = insertSizes->second.begin(); secPairItr != insertSizes->second.end(); secPairItr++) {
                tigpairCheckpointFile << insertSizes->first /*distance*/ << "\t" << pairItr->first << "\t" << secPairItr->first  << "\t" << secPairItr->second.links  << "\t" << secPairItr->second.gaps << "\n";
            }
        }
    }
    // if (distFile.is_open()) {} else {}
    // foreach my $is (sort {$a<=>$b} keys %$track_insert){
    //     print CSV "$is,$track_insert->{$is}\n";
    // }
    tigpairCheckpointFile.close();

}

#endif
