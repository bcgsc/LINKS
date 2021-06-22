//#ifndef LINKS_LINKS_HPP
//#define LINKS_LINKS_HPP
#pragma once

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

class LINKS
{

  public:

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

  //Record get_minimizers();
  void extract_mate_pairs();

  LINKS(InputParser* inputParser);

  ~LINKS();

  InputParser inputParser;


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

private:

    std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo>> matePair;

    std::unordered_map<uint64_t, KmerInfo> trackAll;
    std::unordered_map<std::string, uint64_t> tigLength;

    // store second mates in a set
    std::unordered_set<uint64_t> mates;
   

    // Stage 1 -- populate bloom filter with contigs
    btllib::KmerBloomFilter * myFilter;
    
    std::atomic<bool> fasta{ false };

    static const size_t BUFFER_SIZE = 4;
    static const size_t BLOCK_SIZE = 1;

    std::string seqfile;
    std::queue<std::string> longReads;
    unsigned threads;
    long id;
    size_t buffer_size = 32;
    size_t block_size = 8;
    /*
    
    static const size_t SHORT_MODE_BUFFER_SIZE = 32;
    static const size_t SHORT_MODE_BLOCK_SIZE = 32;

    static const size_t LONG_MODE_BUFFER_SIZE = 4;
    static const size_t LONG_MODE_BLOCK_SIZE = 1;
    */

    btllib::OrderQueueSPMC<Read> input_queue;

    void extract_mate_pair(const std::string& seq);

    static long* ready_blocks_owners()
    {
      thread_local static long var[MAX_SIMULTANEOUS_INDEXLRS];
      return var;
    }

    static std::atomic<long>& last_id()
    {
      static std::atomic<long> var(0);
      return var;
    }

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


  //std::atomic<bool> fasta{ false };
  //OrderQueueSPMC<Read> input_queue;
  //OrderQueueMPSC<Record> output_queue;

  /*using OutputQueueType = decltype(output_queue);
  static std::unique_ptr<OutputQueueType::Block>* ready_blocks_array()
  {
    thread_local static std::unique_ptr<decltype(output_queue)::Block>
      var[MAX_SIMULTANEOUS_INDEXLRS];
    return var;
  }*/

  /*
  static long* ready_blocks_owners()
  {
    thread_local static long var[MAX_SIMULTANEOUS_INDEXLRS];
    return var;
  }

  static std::atomic<long>& last_id()
  {
    static std::atomic<long> var(0);
    return var;
  }
  */

  btllib::SeqReader reader;
  InputWorker input_worker;
  std::vector<ExtractMatePairWorker> extract_mate_pair_workers;
};

inline void
LINKS::InputWorker::work()
{
  std::cout << "here for fun" << std::endl;
  if (links.reader.get_format() == btllib::SeqReader::Format::FASTA) {
    links.fasta = true;
  } else {
    links.fasta = false;
  }

  decltype(links.input_queue)::Block block(links.block_size);
  size_t current_block_num = 0;
  Read read;

  uint counterx = 0;
  for (btllib::SeqReader::Record record; (record = links.reader.read());) {
   //for(auto record : links.reader){
    std::cout << "counterx: " << counterx << std::endl;
    ++counterx;
    std::cout << "counter num: " << record.num << std::endl;
    std::cout << "counter id: " << record.id << std::endl;
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
LINKS::extract_mate_pair(const std::string& seq)
{

  std::cout << "seq size: " << seq.size() << std::endl;
  std::cout << "id: " << std::this_thread::get_id() << std::endl;
} 

inline void
LINKS::ExtractMatePairWorker::work()
{
  std::cout << "here for a intellectual talk" << std::endl;
  decltype(links.input_queue)::Block input_block(links.block_size);
  //decltype(indexlr.output_queue)::Block output_block(indexlr.block_size);
  
  for (;;) {
    if (input_block.current == input_block.count) {
      /*if (output_block.count > 0) {
        output_block.num = input_block.num;
        indexlr.output_queue.write(output_block);
        output_block.current = 0;
        output_block.count = 0;
      }*/
      links.input_queue.read(input_block);
    }
    
    if (input_block.count == 0) {
      /*
      output_block.num = input_block.num;
      output_block.current = 0;
      output_block.count = 0;
      indexlr.output_queue.write(output_block);
      */
      break;
    }
    
    Read& read = input_block.data[input_block.current++];
    //Record record;
    //record.num = read.num;
    //if (links.output_id()) {
    //  record.id = std::move(read.id);
    //}
    /*
    if (links.output_bx()) {
      record.barcode = links.extract_barcode(record.id, read.comment);
    }
    */
    //record.readlen = read.seq.size();

    /*
    check_info(indexlr.verbose && indexlr.k > read.seq.size(),
               "Indexlr: skipped seq " + std::to_string(read.num) +
                 " on line " +
                 std::to_string(read.num * (indexlr.fasta ? 2 : 4) + 2) +
                 "; k (" + std::to_string(indexlr.k) + ") > seq length (" +
                 std::to_string(read.seq.size()) + ")");

    check_info(indexlr.verbose && indexlr.w > read.seq.size() - indexlr.k + 1,
               "Indexlr: skipped seq " + std::to_string(read.num) +
                 " on line " +
                 std::to_string(read.num * (indexlr.fasta ? 2 : 4) + 2) +
                 "; w (" + std::to_string(indexlr.w) + ") > # of hashes (" +
                 std::to_string(read.seq.size() - indexlr.k + 1) + ")");
    */
    if (links.k <= read.seq.size()) {
        links.extract_mate_pair(read.seq);
    } else {
      continue; // nothing
    }

    //output_block.data[output_block.count++] = std::move(record);
  }
  std::cout << "here!" << std::endl;
}

inline LINKS::LINKS(InputParser* inputParser)
        : assemblyFile(inputParser->assemblyFile)
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
        , threads(1)
        , input_queue(buffer_size, block_size)
        , input_worker(*this)
        , extract_mate_pair_workers(
             std::vector<ExtractMatePairWorker>(threads, ExtractMatePairWorker(*this)))
        , reader(std::move(longReads.front()), btllib::SeqReader::Flag::LONG_MODE)
        {
          longReads.pop();
          input_worker.start();
          for (auto& worker : extract_mate_pair_workers) {
              worker.start();
          }
        }

inline LINKS::~LINKS()
{
  //reader.close();
  for (auto& worker : extract_mate_pair_workers) {
    worker.join();
  }
  input_worker.join();
  
}


//#endif