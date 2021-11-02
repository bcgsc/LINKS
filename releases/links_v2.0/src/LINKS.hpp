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
//#include <shared_mutex>
#include <chrono>

class LINKS
{
  public:

    LINKS(InputParser* input_parser);
    void init_bloom_filter();
    void start_read_fasta();
    void start_read_contig();
    void pair_contigs();
    ~LINKS();

    InputParser* input_parser;

    std::string assembly_file;
    std::string fof_file;
    std::vector<uint32_t> distances = {2000,4000};
    std::vector<uint16_t> step_sizes = {2};
    uint64_t k = 15;
    bool verbose = false;
    uint64_t min_links = 5;
    uint64_t min_size = 500;
    float max_link_ratio = 0.3;
    //uint64_t step = 2;
    // Added for MPET
    uint64_t read_length;         // MPET
    float insert_stdev = 0.1;      // MPET (Need to adjust to a wider-range of distances when dealing with MPET) 
    std::string base_name;   // When set, this will override the MPET-induced changes on -e
    uint64_t offset = 0;
    std::string bf_file;
    float fpr = 0.001;
    bool bf_off = false;
    unsigned threads = 3;

    static const size_t MAX_SIMULTANEOUS_INDEXLRS = 128;

    struct BufferMatePairData
    {
      BufferMatePairData() {}

      BufferMatePairData( uint64_t kmer_1_hash, 
                          uint64_t kmer_2_hash, 
                          uint32_t distance)
        : kmer_1_hash(kmer_1_hash)
        , kmer_2_hash(kmer_2_hash)
        , distance(distance)
      {}

      uint64_t kmer_1_hash; 
      uint64_t kmer_2_hash;
      uint32_t distance;
    };

    struct BufferMateData
    {
      BufferMateData() {}

      BufferMateData( uint64_t hash,
                      std::string tig,
                      uint64_t start,
                      uint64_t end,
                      uint64_t multiple,
                      bool orient)
                      : hash(hash) 
                      , tig(tig)
                      , start(start)
                      , end(end)
                      , multiple(multiple)
                      , orient(orient)
      {}

      uint64_t hash;
      std::string tig = "";
      uint64_t start = 0;
      uint64_t end = 0;
      uint64_t multiple = 1;
      bool orient = 0;
    };

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


    struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
      return p.first ^ p.second;
      }
    };

    typedef std::unordered_map<std::pair<uint64_t,uint64_t>, MatePairInfo, pair_hash> mate_pair_type;
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


    std::string seq_file;
    std::queue<std::string> long_reads;
    //long id;

    size_t read_buffer_size = 16;
    size_t read_block_size = 4;

    size_t mate_pair_buffer_size = 6;
    size_t mate_pair_block_size = 1000000;

    btllib::KmerBloomFilter* make_bf(uint64_t bf_elements, InputParser* links_arg_parser);
    void extract_mate_pair( const std::string& seq,
                            btllib::OrderQueueSPMC<BufferMatePairData>::Block& mate_pair_block);
    void populate_mate_info(const std::string& seq,
                            const std::string contig_rank,
                            btllib::OrderQueueSPMC<BufferMateData>::Block& mate_block);
    
    int get_distance_bin(int distance);
    int get_distance(uint64_t insert_size, 
      uint64_t length_i, 
      uint64_t start_i, 
      uint64_t start_j);
    void add_to_pair_map(
      int isz,
      std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo>>>& pair,
      int distance,
      std::string kmer1_name,
      std::string kmer2_name,
      unsigned orient_enum);
    void merge_mate_pair_map(mate_pair_type& own_new_mate_pair);
    void merge_mates_set(std::unordered_set<uint64_t>& own_mates);
    void merge_track_all(std::unordered_map<uint64_t, KmerInfo>& own_track_all);
    void write_from_block_to_map();
    void write_from_block_to_set();
    // helper functions
    uint64_t get_file_size(std::string file_name);
    bool does_file_exist(std::string file_name);
    std::atomic<bool> fasta{ false };

    std::shared_ptr<btllib::KmerBloomFilter> bloom;
    std::shared_ptr<btllib::SeqReader> reader;
    std::shared_ptr<btllib::OrderQueueSPMC<Read>> input_queue;
    
    btllib::OrderQueueSPMC<BufferMatePairData> mate_pair_input_queue;
    btllib::OrderQueueSPMC<BufferMateData> mate_input_queue;
    std::shared_ptr<InputWorker> input_worker;
    std::vector<ExtractMatePairWorker> extract_mate_pair_workers;
    std::vector<PopulateMateInfoWorker> populate_mate_info_workers;

    std::mutex tig_length_mutex;
    
    //mate_pair_type new_mate_pair;
    mate_pair_type test_mate_pair;

    std::unordered_map<uint64_t, KmerInfo> track_all_test;
    std::unordered_map<std::string, uint64_t> tig_length;
    // store second mates in a set
    std::unordered_set<uint64_t> mates;

    std::atomic<std::size_t> mate_pair_current_block_num = { 0 };
    std::atomic<std::size_t> mate_pair_threads_done_writing = { 0 };
    std::atomic<std::size_t> mate_current_block_num = { 0 };
    std::atomic<std::size_t> mate_threads_done_writing = { 0 };
};

inline btllib::KmerBloomFilter*
LINKS::make_bf(uint64_t bf_elements, InputParser* links_arg_parser){
  btllib::KmerBloomFilter * assembly_BF;
    if(links_arg_parser->bf_file != "") {
        std::cout << "A Bloom filter was supplied (" << links_arg_parser->bf_file << ") and will be used instead of building a new one from -f " << links_arg_parser->assembly_file << "\n";
        if(!does_file_exist(links_arg_parser->bf_file)) {
            std::cout << "\nInvalid file: " << links_arg_parser->bf_file <<  " -- fatal\n";
            exit(1);
        } else {
            std::cout << "Checking Bloom filter file " << links_arg_parser->bf_file <<"...ok\n";
        }
        // std::cout << "Loading bloom filter of size " << getFileSize(linksArgParser->bfFile) << " from " << linksArgParser->bfFile << "\n";
        assembly_BF = new btllib::KmerBloomFilter(links_arg_parser->bf_file);
    } else {
        uint64_t m = ceil((-1 * (double)bf_elements * log(links_arg_parser->fpr)) / (log(2) * log(2)));
        //uint64_t rem = 64 - (m % 64);
        m = ((uint64_t)(m / 8) + 1) * 8;
        std::cout << "HASHES CALC: " << std::to_string(((double)m / bf_elements)) << " second: " << std::to_string(((double)m / bf_elements) * log(2)) << "\n";
        unsigned hash_fct = floor(((double)m / bf_elements) * log(2));
        std::cout << "- Number of bfElements: " << bf_elements << "\n"
                    << "- Input file path: " << links_arg_parser->bf_file << "\n"
                    << "- Input file: " << links_arg_parser->assembly_file << "\n"
                    << "- kmersize: " << links_arg_parser->k << "\n"
                    << "- m: " << m << "\n"
                    << "- fpr: " << links_arg_parser->fpr << "\n"
                    << "- hashFct: " << hash_fct << "\n";

        std::string reading_tigbloom_message = "\n\n=>Reading contig/sequence assembly file : " + std::to_string(time(0)) + "\n";
        // assemblyruninfo += reading_tigbloom_message;
        // std::cout << "- Filter output file : " << outFileBf << "\n";
        std::cout << "- Filter output file : " << links_arg_parser->k << "\n";
        assembly_BF = new btllib::KmerBloomFilter(m/8, hash_fct, links_arg_parser->k);
        btllib::SeqReader assembly_reader(links_arg_parser->assembly_file, 8, 1);
        size_t builder = 0;
        for (btllib::SeqReader::Record record; (record = assembly_reader.read());) {
          if(builder % 10 == 0) {
                //std::cout << "reading... builder: " << builder << "\n";
          }
          builder++;
          //std::cout << "here 1\n";
          assembly_BF->insert(record.seq);
          //std::cout << "here 2\n";
        }
        std::string bfmsg = "\n\nWriting Bloom filter to disk (" + links_arg_parser->bf_file + ") : " + std::to_string(time(0)) + "\n";
        // assemblyruninfo += bfmsg;
        std::cout << bfmsg;
        assembly_BF->save(".bloom");
        // std::cout << "Done mybf, printing stats...\n";
        //printBloomStats(*assemblyBF, std::cout);
        
    }
    return assembly_BF;
}

inline LINKS::LINKS(InputParser* input_parser)
        : input_parser(input_parser)
        , assembly_file(input_parser->assembly_file)
        , long_reads(input_parser->long_reads)
        , distances(input_parser->distances)
        , k(input_parser->k)
        , verbose(input_parser->verbose)
        , min_size(input_parser->min_size)
        , step_sizes(input_parser->step_sizes)
        , insert_stdev(input_parser->insert_stdev)
        , base_name(input_parser->base_name)
        , offset(input_parser->offset)
        , fpr(input_parser->fpr)
        , bf_file(input_parser->bf_file)
        , bf_off(input_parser->bf_off)
        , threads(input_parser->thread)
        , input_queue(std::shared_ptr<btllib::OrderQueueSPMC<Read>>(new btllib::OrderQueueSPMC<Read>(read_buffer_size, read_block_size)))
        , mate_pair_input_queue(mate_pair_buffer_size, mate_pair_block_size)
        , mate_input_queue(mate_pair_buffer_size, mate_pair_block_size)
        , input_worker(std::shared_ptr<InputWorker>(new InputWorker(*this)))
        {}

inline void
LINKS::init_bloom_filter(){
  int64_t bf_elements = get_file_size(input_parser->assembly_file);
  bloom = std::shared_ptr<btllib::KmerBloomFilter>(make_bf(bf_elements,input_parser));
}

inline void
LINKS::write_from_block_to_map(){
  btllib::OrderQueueSPMC<BufferMatePairData>::Block mate_pair_block(mate_pair_block_size);

  //test_mate_pair.reserve(1000000000);
  //mates.reserve(1000000000); 

  for (;;) {
    if (mate_pair_block.current == mate_pair_block.count) {
      mate_pair_input_queue.read(mate_pair_block);
      /*
      if (mate_pair_block.count != 0) {
        std::cout << "block num after read: " << mate_pair_block.num << std::endl;
      }else{
        std::cout << "block num after read2: " << mate_pair_block.num << std::endl;
      }
      */
    }
    
    if (mate_pair_block.count == 0) {
      break;
    }
    
    BufferMatePairData& mate_data = mate_pair_block.data[mate_pair_block.current++];
    

    if(test_mate_pair.find(std::make_pair(mate_data.kmer_1_hash,mate_data.kmer_2_hash)) == test_mate_pair.end()){
      test_mate_pair[std::make_pair(mate_data.kmer_1_hash,mate_data.kmer_2_hash)] = MatePairInfo(false, mate_data.distance);
      mates.insert(mate_data.kmer_1_hash); // with new data structure have to insert both
      mates.insert(mate_data.kmer_2_hash);
    }
  }
  //std::cout << "writer leaves\n";
  std::cout << "test_mate_pair.size() " << test_mate_pair.size() << std::endl;
  std::cout << "mates.size() " << mates.size() << std::endl;

}
inline void
LINKS::write_from_block_to_set(){
  btllib::OrderQueueSPMC<BufferMateData>::Block mate_block(mate_pair_block_size);

  // for fast testing
  // track_all_test.reserve(10000000);

  for (;;) {
    if (mate_block.current == mate_block.count) {
      mate_input_queue.read(mate_block);
      /*
      if (mate_block.count != 0) {
        std::cout << "block num after read: " << mate_block.num << std::endl;
      }else{
        std::cout << "block num after read2: " << mate_block.num << std::endl;
      }
      */
    }
    
    if (mate_block.count == 0) {
      break;
    }
    
    BufferMateData& mate_data = mate_block.data[mate_block.current++];

    if(track_all_test.find(mate_data.hash) == track_all_test.end()) {
      track_all_test[mate_data.hash] = KmerInfo(mate_data.tig, mate_data.start, mate_data.end, mate_data.multiple, mate_data.orient);
    } else {
      track_all_test[mate_data.hash].multiple +=1;
    }
  }
  std::cout << "writer leaves\n";
  std::cout << "test_mate_pair.size() " << test_mate_pair.size() << std::endl;
  std::cout << "mates.size() " << mates.size() << std::endl;
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
  btllib::OrderQueueSPMC<Read>::Block block(links.read_block_size); 

  size_t current_block_num = 0;
  size_t read_c = 0;
  Read read;

  size_t counterx = 0;
  for(auto record : (*(links.reader))){
    read_c++;
    if(read_c % 10000 == 0){
      std::cout << "read_counter: " << read_c << std::endl;
    }

    block.data[block.count++] = Read(record.num, //read_c instead of record.num
                                    std::move(record.id),
                                    std::move(record.comment),
                                    std::move(record.seq));
    if (block.count == links.read_block_size) {
      block.num = current_block_num++;
      links.input_queue->write(block);
      block.count = 0;
    }
  }

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

  reader = std::shared_ptr<btllib::SeqReader>(new btllib::SeqReader(std::move(long_reads.front()), btllib::SeqReader::Flag::LONG_MODE));
  long_reads.pop();
  
  extract_mate_pair_workers = std::vector<ExtractMatePairWorker>(threads, ExtractMatePairWorker(*this));

  input_worker->start();
  std::thread writer_thread(&LINKS::write_from_block_to_map, this);

  for (auto& worker : extract_mate_pair_workers) {
    std::cout << "Here1\n";
    worker.start();
  }
  // wait
  for (auto& worker : extract_mate_pair_workers) {
    std::cout << "Here2\n";
    worker.join();
  }
  std::cout << "Here3\n";
  writer_thread.join();
  std::cout << "Here4\n";
  input_worker->join();
  std::cout << "Here5\n";
}
inline void 
LINKS::start_read_contig(){
  //trackAll.reserve(mates.size());
  reader.reset();
  reader = std::shared_ptr<btllib::SeqReader>(new btllib::SeqReader(assembly_file, btllib::SeqReader::Flag::LONG_MODE));
  
  input_worker.reset();
  input_worker = std::shared_ptr<InputWorker>(new InputWorker(*this));

  input_queue.reset();
  input_queue = std::shared_ptr<btllib::OrderQueueSPMC<Read>>(new btllib::OrderQueueSPMC<Read>(read_buffer_size, read_block_size));

  input_worker->start();

  std::thread writer_thread(&LINKS::write_from_block_to_set, this);

  populate_mate_info_workers = std::vector<PopulateMateInfoWorker>(threads, PopulateMateInfoWorker(*this));

  for (auto& worker : populate_mate_info_workers) {
    worker.start();
  }
  for (auto& worker : populate_mate_info_workers) {
    worker.join();
  }
  std::cout << "all workers joined" << std::endl;
  writer_thread.join();
  std::cout << "writer joined" << std::endl;
  input_worker->join();
  std::cout << "input worker joined" << std::endl;
  input_queue.reset();
  
  std::cout << "track_all_test.size() " << track_all_test.size() << std::endl;
}

inline void
LINKS::extract_mate_pair(const std::string& seq,
  btllib::OrderQueueSPMC<BufferMatePairData>::Block& mate_pair_block)
{
  uint step_index = 0;
  for(auto distance : distances){
    uint cur_step_size = step_sizes[step_index];
    step_index++;
    //std::cout << "distance: " << distance << " step size: " << cur_step_size << std::endl;
    uint64_t delta = distance - k;
    int break_flag = 0;
    bool reverse_exists = false;
    btllib::NtHash ntHash(seq, bloom->get_hash_num(), k, offset);
    btllib::NtHash ntHash_lead(seq, bloom->get_hash_num(), k, delta + offset);

    uint64_t hash_a, hash_b;    

    for (size_t i = 0; ntHash.roll() && ntHash_lead.roll(); i+=cur_step_size) {
        // roll for the number of steps
        break_flag = 0;
        reverse_exists = false;
        // for step ----
        for(uint j = 1; j < cur_step_size; j++) {
            if(!ntHash_lead.roll() || !ntHash.roll()) {
                break_flag = 1;
            }
        }  
        if(break_flag){break;}
        // for step ----

        if(ntHash.get_forward_hash() < ntHash_lead.get_reverse_hash()){
          hash_a = ntHash_lead.get_reverse_hash();
          hash_b = ntHash.get_reverse_hash();
        }else{
          hash_a = ntHash.get_forward_hash();
          hash_b = ntHash_lead.get_forward_hash();
        }

        if(bloom->contains(ntHash.hashes()) && bloom->contains(ntHash_lead.hashes())) { // May need to change with forward reverse hashes
          mate_pair_block.data[mate_pair_block.count++] = BufferMatePairData(hash_a,
                                  hash_b,
                                  distance);
          if (mate_pair_block.count == mate_pair_block_size) {
            mate_pair_block.num = mate_pair_current_block_num++;
            mate_pair_input_queue.write(mate_pair_block);
            mate_pair_block.count = 0;
          } 
        }
    }
  }
} 

inline void
LINKS::populate_mate_info(const std::string& seq,
                          const std::string contig_rank,
                          btllib::OrderQueueSPMC<BufferMateData>::Block& mate_info_block){
  
  btllib::NtHash ntHash_contig(seq, bloom->get_hash_num(), k); // hashFunc can be 1 after first step

  int breakFlag = 0;
  for (size_t i = 0; ntHash_contig.roll(); i+=1) {

      i = ntHash_contig.get_pos();
      if(mates.find(ntHash_contig.get_forward_hash()) != mates.end()){
          mate_info_block.data[mate_info_block.count++] = BufferMateData(
                      ntHash_contig.get_forward_hash(),
                      contig_rank, 
                      i, 
                      i + k, 
                      1, 
                      false);
          if (mate_info_block.count == mate_pair_block_size) {
            mate_info_block.num = mate_current_block_num++;
            mate_input_queue.write(mate_info_block);
            mate_info_block.count = 0;
          }
        }

        // Reverse part
        if(mates.find(ntHash_contig.get_reverse_hash()) != mates.end()){
          mate_info_block.data[mate_info_block.count++] = BufferMateData(
              ntHash_contig.get_reverse_hash(),
              contig_rank, 
              i, 
              i + k, 
              1, 
              true);
          if (mate_info_block.count == mate_pair_block_size) {
            mate_info_block.num = mate_current_block_num++;
            mate_input_queue.write(mate_info_block);
            mate_info_block.count = 0;
          }
        }
    }
}
inline void
LINKS::ExtractMatePairWorker::work()
{
  btllib::OrderQueueSPMC<Read>::Block input_block(links.read_block_size);

  btllib::OrderQueueSPMC<BufferMatePairData>::Block mate_pair_block(links.mate_pair_block_size); 

  for (;;) {
  //for (int ll=0; ll<100000; ll++) {
    if (input_block.current == input_block.count) {
      links.input_queue->read(input_block);
    }
    
    if (input_block.count == 0) {
      break;
    }
    
    Read& read = input_block.data[input_block.current++];

    if (links.k <= read.seq.size()) {
        links.extract_mate_pair(read.seq,mate_pair_block);
    } else {
      continue; // nothing
    }
  }
  if (mate_pair_block.count > 0) {
    mate_pair_block.num = links.mate_pair_current_block_num++;
    //std::cout << "block num after write2: " << mate_pair_block.num << std::endl;
    links.mate_pair_input_queue.write(mate_pair_block);
  }
  links.mate_pair_threads_done_writing++;
  if(links.mate_pair_threads_done_writing == links.threads){
    for (size_t i = 0; i < links.threads; i++)
    {
      mate_pair_block.num = links.mate_pair_current_block_num++;
      mate_pair_block.current = 0;
      mate_pair_block.count = 0;
      links.mate_pair_input_queue.write(mate_pair_block);
    }
  }
}

inline void
LINKS::PopulateMateInfoWorker::work()
{
  btllib::OrderQueueSPMC<Read>::Block input_block(links.read_block_size);
  btllib::OrderQueueSPMC<BufferMateData>::Block mate_block(links.mate_pair_block_size); 
  std::unordered_map<uint64_t, KmerInfo> own_track_all;

  std::cout << "t test 10\n";
  for (;;) {
    if (input_block.current == input_block.count) {
      links.input_queue->read(input_block);
    }
    if (input_block.count == 0) {
      break;
    }
    Read& read = input_block.data[input_block.current++];

    if (links.k <= read.seq.size()) {
        links.tig_length_mutex.lock(); 
        links.tig_length[std::to_string(read.num+1)] = read.seq.size();
        links.tig_length_mutex.unlock();
        links.populate_mate_info(read.seq,std::to_string(read.num+1), mate_block);
        continue;
    } else {
      continue; // nothing
    }
  }
  if (mate_block.count > 0) {
    mate_block.num = links.mate_current_block_num++;
    links.mate_input_queue.write(mate_block);
  }
  links.mate_threads_done_writing++;
  if(links.mate_threads_done_writing == links.threads){
  
  for (size_t i = 0; i < links.threads; i++)
  {
    mate_block.num = links.mate_current_block_num++;
    mate_block.current = 0;
    mate_block.count = 0;
    links.mate_input_queue.write(mate_block);
  }
  
  }
  std::cout << "populate mate info thread done\n";
}

inline LINKS::~LINKS()
{
}

inline uint64_t 
LINKS::get_file_size(std::string file_name)
{
    // This buffer is a stat struct that the information is placed concerning the file.
    struct stat stat_buf;
    // Stat method returns true if successfully completed
    int rc = stat(file_name.c_str(), &stat_buf);
    // st_size holds the total size of the file in bytes
    return rc == 0 ? stat_buf.st_size : -1;
}

inline bool 
LINKS::does_file_exist(std::string file_name)
{
    std::ifstream infile(file_name);
    return infile.good();
}

inline void 
LINKS::add_to_pair_map(
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
LINKS::get_distance_bin(int distance)
{
    return distance < 0 ? -1 : distance == 10 ? 10 : distance < 500 ? 500 : distance < 5000 ? 5000 : 10000;
}

inline int 
LINKS::get_distance(uint64_t insert_size, uint64_t length_i, uint64_t start_i, uint64_t start_j) {

   int insert_span = (length_i - start_i) + start_j;
   int gap_or_overlap = insert_size - insert_span;

   return gap_or_overlap;
}
inline void
LINKS::pair_contigs()
{
  std::cout << "in real pair contigs\n";
  std::string issues = input_parser->base_name + ".pairing_issues";
  std::string distribution = input_parser->base_name + ".pairing_distribution.csv";
  std::string tigpair_checkpoint = input_parser->base_name + ".tigpair_checkpoint.tsv"; // add a checkpoint file, prevent re-running LINKS from scratch if crash
  std::string simplepair_checkpoint = input_parser->base_name + ".simplepair_checkpoint.tsv"; // add a checkpoint file, prevent re-running LINKS from scratch if cras
  uint64_t totalPairs = 0;
  uint64_t ct_illogical = 0, ct_ok_contig = 0, ct_ok_pairs = 0, ct_problem_pairs = 0, ct_iz_issues = 0, ct_single = 0, ct_multiple = 0, ct_both = 0, track_insert = 0;
  std::unordered_map<uint64_t, uint64_t> ct_single_hash, ct_both_hash, ct_illogical_hash, ct_ok_contig_hash, ct_ok_pairs_hash, ct_problem_pairs_hash, ct_iz_issues_hash;
  // Mapping of tiga_head -> insertSize -> tigb_head -> links & gaps
  std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo> > > pair;
  //std::unordered_map<std::string, std::unordered_map<std::string, PairLinkInfo> >simplePair;
  std::unordered_map<std::string, PairLinkInfo> err;
  //std::string order1;
  //std::string order2;
  if(verbose) std::cout << "Pairing contigs...\n";
  std::ofstream issues_file;
  issues_file.open (issues);
  int64_t insert_size = 0;
  int min_allowed = 0;
  uint low_iz = 0;
  uint up_iz = 0;
  int distance;
  int isz;

  std::string tig_a, tig_b, ftig_a, ftig_b, rtig_a, rtig_b;
  uint64_t A_length = 0, A_start = 0, A_end = 0, B_length = 0, B_start = 0, B_end = 0;

  KmerInfo kmer1, kmer2;

  size_t counter = 0;
  size_t percent_size = test_mate_pair.size() / 100;
  mate_pair_type::iterator mate_pair_iterator;
  for (mate_pair_iterator = test_mate_pair.begin(); mate_pair_iterator != test_mate_pair.end(); mate_pair_iterator++)
  {
    if(counter % percent_size == 0){
      std::cout << "Done: %" << uint(counter / percent_size) << std::endl;
    }
    ++counter;
    //++counter_1;
    if (mate_pair_iterator->second.seen == false &&                                      //matepair is not seen
        track_all_test.find(mate_pair_iterator->first.first) != track_all_test.end() &&  //first mate is tracked
        track_all_test[mate_pair_iterator->first.first].multiple == 1 &&                 //first mate seen once
        track_all_test.find(mate_pair_iterator->first.second) != track_all_test.end() && //second mate is tracked
        track_all_test[mate_pair_iterator->first.second].multiple == 1)
    {

      //++counter_2;

      mate_pair_iterator->second.seen = true;

      insert_size = test_mate_pair[std::make_pair(mate_pair_iterator->first.first, mate_pair_iterator->first.second)].insert_size;

      min_allowed = -1 * (insert_stdev * insert_size); // check int
      low_iz = insert_size + min_allowed;             // check int
      up_iz = insert_size - min_allowed;              // check int

      //if(verbose) std::cout << "Pair read1Hash=" << matePairItr->first << " read2Hash=" << mateListItr->first << "\n";

      if (track_all_test[mate_pair_iterator->first.first].tig != "" && track_all_test[mate_pair_iterator->first.second].tig != "")
      { //double check if tig names not null
        //++counter_3;
        ct_both++;
        if (ct_both_hash.find(insert_size) == ct_both_hash.end())
        {
          ct_both_hash[insert_size] = 1;
        }
        else
        {
          ct_both_hash[insert_size] = ct_both_hash[insert_size] + 1;
        }

        kmer1 = track_all_test[mate_pair_iterator->first.first];
        kmer2 = track_all_test[mate_pair_iterator->first.second];

        if (kmer1.tig != kmer2.tig)
        { // paired reads located on <> contigs

          // MURATHAN DEBUG 11.5.21
          if (!kmer1.orient)
          { // if kmer1 is forward
            if (!kmer2.orient)
            { // if kmer2 is forward
              //++counter_4;
              distance = get_distance(insert_size, tig_length[kmer1.tig], kmer1.start, kmer2.start);
              if (distance > min_allowed && distance < insert_size)
              {
                isz = get_distance_bin(distance);
                add_to_pair_map(isz, pair, distance, kmer1.tig, kmer2.tig, 0);
              }
            }
            else
            {
              //++counter_5;
              distance = get_distance(insert_size, tig_length[kmer1.tig], kmer1.start, tig_length[kmer2.tig] - kmer2.end);
              if (distance > min_allowed && distance < insert_size)
              {
                isz = get_distance_bin(distance);
                add_to_pair_map(isz, pair, distance, kmer1.tig, kmer2.tig, 1);
              }
            }
          }
          else
          { // if kmer1 is reverse
            if (!kmer2.orient)
            {
              //++counter_6;
              distance = get_distance(insert_size, tig_length[kmer1.tig], tig_length[kmer1.tig] - kmer1.end, kmer2.start);
              if (distance > min_allowed && distance < insert_size)
              {
                isz = get_distance_bin(distance);
                add_to_pair_map(isz, pair, distance, kmer1.tig, kmer2.tig, 2);
              }
            }
            else
            { // if kmer2 is reverse
              //++counter_7;
              distance = get_distance(insert_size, tig_length[kmer2.tig], kmer2.end, kmer1.end);
              if (distance > min_allowed && distance < insert_size)
              {
                isz = get_distance_bin(distance);
                add_to_pair_map(isz, pair, distance, kmer1.tig, kmer2.tig, 3);
              }
            }
          }
        }
        else
        { // Clone, paired reads located on the same contig -- could be used to investigate misassemblies
          if (verbose)
            std::cout << "Pair (" << mate_pair_iterator->first.first << " and " << mate_pair_iterator->first.second << ") located on same contig " << tig_a << " (" << A_length << " nt)\n";
            //std::cout << "Pair (" << matePairItr->first << " and " << mateListItr->first << ") located on same contig " << tig_a << " (" << A_length << " nt)\n";
          uint64_t pet_size = 0;
          if (A_start > B_start && (B_start < B_end) && (A_start > A_end))
          { // B --> <-- A
            pet_size = A_start - B_start;
            track_insert += pet_size;
            if (pet_size >= low_iz && pet_size <= up_iz)
            {
              ct_ok_contig++;
              if (ct_ok_contig_hash.find(insert_size) == ct_ok_contig_hash.end())
              {
                ct_ok_contig_hash[insert_size] = 1;
              }
              else
              {
                ct_ok_contig_hash[insert_size] = ct_ok_contig_hash[insert_size] + 1;
              }
            }
            else
            {
              issues_file << "Pairs unsatisfied in distance within a contig.  Pair (" << mate_pair_iterator->first.first << " and " << mate_pair_iterator->first.second << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << " CALCULATED DISTANCE APART: " << pet_size << "\n";
              //issues_file << "Pairs unsatisfied in distance within a contig.  Pair (" << matePairItr->first << " - " << mateListItr->first << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << " CALCULATED DISTANCE APART: " << pet_size << "\n";
              ct_iz_issues++;
              if (ct_iz_issues_hash.find(insert_size) == ct_iz_issues_hash.end())
              {
                ct_iz_issues_hash[insert_size] = 1;
              }
              else
              {
                ct_iz_issues_hash[insert_size] = ct_iz_issues_hash[insert_size] + 1;
              }
            }
          }
          else if (B_start > A_start && (B_start > B_end) && (A_start < A_end))
          { // A --> <-- B
            pet_size = B_start - A_start;
            track_insert += pet_size;
            if (pet_size >= low_iz && pet_size <= up_iz)
            {
              ct_ok_contig++;
              if (ct_ok_contig_hash.find(insert_size) == ct_ok_contig_hash.end())
              {
                ct_ok_contig_hash[insert_size] = 1;
              }
              else
              {
                ct_ok_contig_hash[insert_size] = ct_ok_contig_hash[insert_size] + 1;
              }
            }
            else
            {
              issues_file << "Pairs unsatisfied in distance within a contig.  Pair (" << mate_pair_iterator->first.first << " and " << mate_pair_iterator->first.second << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << "\n";
              //issues_file << "Pairs unsatisfied in distance within a contig.  Pair (" << matePairItr->first << " - " << mateListItr->first << ") on contig " << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << "\n";
              ct_iz_issues++;
              if (ct_iz_issues_hash.find(insert_size) == ct_iz_issues_hash.end())
              {
                ct_iz_issues_hash[insert_size] = 1;
              }
              else
              {
                ct_iz_issues_hash[insert_size] = ct_iz_issues_hash[insert_size] + 1;
              }
            }
          }
          else
          {
            ct_illogical++;
            if (ct_illogical_hash.find(insert_size) == ct_illogical_hash.end())
            {
              ct_illogical_hash[insert_size] = 1;
            }
            else
            {
              ct_illogical_hash[insert_size] = ct_illogical_hash[insert_size] + 1;
            }
            // FOLLOWING IS NOT A DEBUGGING PRINT
            //issues_file << "Pairs unsatisfied in pairing logic within a contig.  Pair (" << matePairItr->first << " - " << mateListItr->first << ") on contig" << tig_a << " (" << A_length << " nt) Astart:" << A_start << " Aend:" << A_end << " Bstart:" << B_start << " Bend:" << B_end << "\n";
          }
        }
      }
      else
      { // both pairs assembled
        ct_single++;
        if (ct_single_hash.find(insert_size) == ct_single_hash.end())
        {
          ct_single_hash[insert_size] = 1;
        }
        else
        {
          ct_single_hash[insert_size] = ct_single_hash[insert_size] + 1;
        }
      }
    }
    else
    { // if unseen
      // std::cout << "UNSEEN\n";
      if (test_mate_pair[std::make_pair(mate_pair_iterator->first.first, mate_pair_iterator->first.second)].seen == false)
      {
        ct_multiple++;
      }
    }
  } // pairing read b
    // pairing read a

    // Summary of the contig pair issues
    //std::cout << "------------- Putative issues with contig pairing - Summary  ----------------\n";
    //std::cout << "err map size: " << err.size() << "\n"; 
    
    //sortErr(err); TODO: uncomment
    std::unordered_map<std::string, PairLinkInfo>::iterator err_itr;
    for(err_itr = err.begin(); err_itr != err.end(); err_itr++) {
        double mean_iz = 0;
        if(err_itr->second.links) {
            mean_iz = err_itr->second.gaps / err_itr->second.links;
        }
        // std::cout << "Pair " << err_itr->first << " has " << err_itr->second.getLinks() << " Links and mean distance of = " << mean_iz << "\n";
    }

    uint64_t satisfied = ct_ok_pairs + ct_ok_contig;
    uint64_t unsatisfied = ct_problem_pairs + ct_iz_issues + ct_illogical;
    uint64_t ct_both_reads = ct_both * 2;

    
    //std::cout << "THESE ARE THE FILTERINGS:\n"<< "filter 1: "<< std::to_string(filter1) << "\n" << "filter 2: "<< std::to_string(filter2) << "\n" << "filter 3: "<< std::to_string(filter3) << "\n" << "filter 4: "<< std::to_string(filter4) << "\n";
    //std::cout << "THESE ARE THE COUNTERS:\n" << std::to_string(CheckCounterBase) << "\n0 " << std::to_string(Check0Counter) << "\n1 " <<  std::to_string(Check1Counter) << "\n2 " <<  std::to_string(Check2Counter) << "\n3 " <<  std::to_string(Check3Counter) << "\n4 " <<  std::to_string(Check4Counter) << "\n5 " <<  std::to_string(Check5Counter) << "\n6 " <<  std::to_string(Check6Counter) << "\n7 " <<  std::to_string(Check7Counter) << "\n8 " <<  std::to_string(Check8Counter) << "\n9 " <<  std::to_string(Check9Counter) << "\n10 " << std::to_string(Check10Counter) << "\n11 " << std::to_string(Check11Counter) << "\n12 " << std::to_string(Check12Counter) << "\n13 " << std::to_string(Check13Counter) << "\n14 " << std::to_string(Check14Counter) << "\n15 " << std::to_string(Check15Counter) << "\n16 " << std::to_string(Check16Counter) << "\n17 " << std::to_string(Check17Counter) << "\n18 " << std::to_string(Check18Counter) << "\n19" << std::to_string(Check19Counter) << "\n20 " << std::to_string(Check20Counter) << "\n21 " << std::to_string(Check21Counter) << "\n22 " << std::to_string(Check22Counter) << "\n23 " << std::to_string(Check23Counter) << "\n24 " << std::to_string(Check24Counter) << "\n25 " << std::to_string(Check25Counter) << "\n26 " << std::to_string(Check26Counter) << "\n";
    std::cout << "\n===========PAIRED K-MER STATS===========\n";
    //std::cout << "Total number of pairs extracted from -s " << long_readsFile << " " << totalPairs << "\n";
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
    std::unordered_map<uint64_t, uint64_t>::iterator itr_IS;
    std::cout << "ct_both_hash map size: " << err.size() << "\n"; 
    for(itr_IS = ct_both_hash.begin(); itr_IS != ct_both_hash.end(); itr_IS++) {
        std::cout <<  "--------k-mers separated by "<< itr_IS->first << " bp (outer distance)--------\n";
        //int64_t maopt = -1 * (insert_stdev * itr_IS->first);
        //int64_t low_izopt = itr_IS->first + maopt;
        //int64_t up_izopt = itr_IS->first - maopt;
        /*
        std::cout <<  "MIN: " << low_izopt << " MAX: " << up_izopt << "  as defined by  " << itr_IS->first << "  *  " << insert_stdev << " \n";
        std::cout <<  "At least one sequence/pair missing:  " << ct_single_hash[itr_IS->first] << " \n";
        std::cout <<  "Assembled pairs:  " << itr_IS->second << " \n";
        std::cout <<  "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target:  " << ct_ok_contig_hash[itr_IS->first] << " \n";
        std::cout <<  "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds):  " << ct_iz_issues_hash[itr_IS->first] << " \n";
        std::cout <<  "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->):  " << ct_illogical_hash[itr_IS->first] << " \n";
        std::cout <<  "\t---\n";
        std::cout <<  "\tSatisfied in distance/logic within a given contig pair (pre-scaffold):  " << ct_ok_pairs_hash[itr_IS->first] << " \n";
        std::cout <<  "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds):  " << ct_problem_pairs_hash[itr_IS->first] << " \n";
        */
    }
    std::cout << "============================================\n";
    std::ofstream dist_file;
    dist_file.open ("distfile.txt");
    // if (dist_file.is_open()) {} else {}
    // foreach my $is (sort {$a<=>$b} keys %$track_insert){
    //     print CSV "$is,$track_insert->{$is}\n";
    // }
    dist_file.close();

    // TIGPAIR CHECKPOINT
    std::ofstream tigpair_checkpoint_file;
    tigpair_checkpoint_file.open (tigpair_checkpoint);
    std::unordered_map<std::string, std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo> > >::iterator pair_itr;
    std::cout << "size of TIGPAIR: " << pair.size() << "\n";
    for(pair_itr = pair.begin(); pair_itr != pair.end(); pair_itr++) {
        std::unordered_map<int64_t, std::unordered_map<std::string, PairLinkInfo> >::iterator insert_sizes;
        for(insert_sizes = pair_itr->second.begin(); insert_sizes != pair_itr->second.end(); insert_sizes++) {
            std::unordered_map<std::string, PairLinkInfo>::iterator sec_pair_itr;
            for(sec_pair_itr = insert_sizes->second.begin(); sec_pair_itr != insert_sizes->second.end(); sec_pair_itr++) {
                tigpair_checkpoint_file << insert_sizes->first /*distance*/ << "\t" << pair_itr->first << "\t" << sec_pair_itr->first  << "\t" << sec_pair_itr->second.links  << "\t" << sec_pair_itr->second.gaps << "\n";
            }
        }
    }
    // if (dist_file.is_open()) {} else {}
    // foreach my $is (sort {$a<=>$b} keys %$track_insert){
    //     print CSV "$is,$track_insert->{$is}\n";
    // }
    tigpair_checkpoint_file.close();

}

#endif
