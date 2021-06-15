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
    size_t buffer_size = 8;
    size_t block_size = 2;
    /*
    
    static const size_t SHORT_MODE_BUFFER_SIZE = 32;
    static const size_t SHORT_MODE_BLOCK_SIZE = 32;

    static const size_t LONG_MODE_BUFFER_SIZE = 4;
    static const size_t LONG_MODE_BLOCK_SIZE = 1;
    */

    btllib::OrderQueueSPMC<Read> input_queue;

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

/*
  std::atomic<bool> fasta{ false };
  OrderQueueSPMC<Read> input_queue;
  OrderQueueMPSC<Record> output_queue;

  using OutputQueueType = decltype(output_queue);
  static std::unique_ptr<OutputQueueType::Block>* ready_blocks_array()
  {
    thread_local static std::unique_ptr<decltype(output_queue)::Block>
      var[MAX_SIMULTANEOUS_INDEXLRS];
    return var;
  }

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
  if (links.reader.get_format() == btllib::SeqReader::Format::FASTA) {
    links.fasta = true;
  } else {
    links.fasta = false;
  }

  decltype(links.input_queue)::Block block(links.block_size);
  size_t current_block_num = 0;
  Read read;

  for (btllib::SeqReader::Record record; (record = links.reader.read());) {
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
LINKS::ExtractMatePairWorker::work()
{
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
/*
inline LINKS::LINKS(std::string seqfile,
                        const size_t k,
                        const size_t w,
                        const unsigned flags,
                        const unsigned threads,
                        const bool verbose,
                        const BloomFilter& bf1,
                        const BloomFilter& bf2)
  : seqfile(std::move(seqfile))
  , k(k)
  , w(w)
  , flags(flags)
  , threads(threads)
  , verbose(verbose)
  , id(++last_id())
  , buffer_size(short_mode() ? SHORT_MODE_BUFFER_SIZE : LONG_MODE_BUFFER_SIZE)
  , block_size(short_mode() ? SHORT_MODE_BLOCK_SIZE : LONG_MODE_BLOCK_SIZE)
  , filter_in_bf(filter_in() ? bf1 : Indexlr::dummy_bf())
  , filter_out_bf(filter_out() ? filter_in() ? bf2 : bf1 : Indexlr::dummy_bf())
  , filter_in_enabled(filter_in())
  , filter_out_enabled(filter_out())
  , input_queue(buffer_size, block_size)
  , output_queue(buffer_size, block_size)
  , reader(this->seqfile,
           short_mode() ? SeqReader::Flag::SHORT_MODE
                        : SeqReader::Flag::LONG_MODE)
  , input_worker(*this)
  , minimize_workers(
      std::vector<MinimizeWorker>(threads, MinimizeWorker(*this)))
{
  check_error(!short_mode() && !long_mode(),
              "Indexlr: no mode selected, either short or long mode flag must "
              "be provided.");
  check_error(threads == 0,
              "Indexlr: Number of processing threads cannot be 0.");
  input_worker.start();
  for (auto& worker : minimize_workers) {
    worker.start();
  }
}
*/

inline LINKS::~LINKS()
{
  //reader.close();
  for (auto& worker : extract_mate_pair_workers) {
    worker.join();
  }
  input_worker.join();
  
}

/*
inline static void
filter_hashed_kmer(LINKS::HashedKmer& hk,
                   bool filter_in,
                   bool filter_out,
                   const BloomFilter& filter_in_bf,
                   const BloomFilter& filter_out_bf)
{
  if (filter_in && filter_out) {
    std::vector<uint64_t> tmp;
    tmp = { hk.min_hash };
    if (!filter_in_bf.contains(tmp) || filter_out_bf.contains(tmp)) {
      hk.min_hash = std::numeric_limits<uint64_t>::max();
    }
  } else if (filter_in) {
    if (!filter_in_bf.contains({ hk.min_hash })) {
      hk.min_hash = std::numeric_limits<uint64_t>::max();
    }
  } else if (filter_out) {
    if (filter_out_bf.contains({ hk.min_hash })) {
      hk.min_hash = std::numeric_limits<uint64_t>::max();
    }
  }
}

inline static void
calc_minimizer(const std::vector<LINKS::HashedKmer>& hashed_kmers_buffer,
               const LINKS::Minimizer*& min_current,
               const size_t idx,
               ssize_t& min_idx_left,
               ssize_t& min_idx_right,
               ssize_t& min_pos_prev,
               const size_t w,
               std::vector<LINKS::Minimizer>& minimizers)
{
  min_idx_left = ssize_t(idx + 1 - w);
  min_idx_right = ssize_t(idx + 1);
  const auto& min_left =
    hashed_kmers_buffer[min_idx_left % hashed_kmers_buffer.size()];
  const auto& min_right =
    hashed_kmers_buffer[(min_idx_right - 1) % hashed_kmers_buffer.size()];

  if (min_current == nullptr || min_current->pos < min_left.pos) {
    min_current = &min_left;
    // Use of operator '<=' returns the minimum that is furthest from left.
    for (ssize_t i = min_idx_left; i < min_idx_right; i++) {
      const auto& min_i = hashed_kmers_buffer[i % hashed_kmers_buffer.size()];
      if (min_i.min_hash <= min_current->min_hash) {
        min_current = &min_i;
      }
    }
  } else if (min_right.min_hash <= min_current->min_hash) {
    min_current = &min_right;
  }
  if (ssize_t(min_current->pos) > min_pos_prev &&
      min_current->min_hash != std::numeric_limits<uint64_t>::max()) {
    min_pos_prev = ssize_t(min_current->pos);
    minimizers.push_back(*min_current);
  }
}

inline std::vector<LINKS::Minimizer>
LINKS::minimize(const std::string& seq) const
{
  if ((k > seq.size()) || (w > seq.size() - k + 1)) {
    return {};
  }
  std::vector<Minimizer> minimizers;
  minimizers.reserve(2 * (seq.size() - k + 1) / w);
  std::vector<HashedKmer> hashed_kmers_buffer(w + 1);
  ssize_t min_idx_left, min_idx_right, min_pos_prev = -1;
  const Minimizer* min_current = nullptr;
  size_t idx = 0;
  for (NtHash nh(seq, 2, k); nh.roll(); ++idx) {
    auto& hk = hashed_kmers_buffer[idx % hashed_kmers_buffer.size()];

    hk = HashedKmer(nh.hashes()[0],
                    nh.hashes()[1],
                    nh.get_pos(),
                    nh.forward(),
                    output_seq() ? seq.substr(nh.get_pos(), k) : "");

    filter_hashed_kmer(
      hk, filter_in(), filter_out(), filter_in_bf.get(), filter_out_bf.get());

    if (idx + 1 >= w) {
      calc_minimizer(hashed_kmers_buffer,
                     min_current,
                     idx,
                     min_idx_left,
                     min_idx_right,
                     min_pos_prev,
                     w,
                     minimizers);
    }
  }
  return minimizers;
}
*/
inline void
LINKS::extract_mate_pairs()
{
  /*
  if (ready_blocks_owners()[id % MAX_SIMULTANEOUS_INDEXLRS] != id) {
    ready_blocks_array()[id % MAX_SIMULTANEOUS_INDEXLRS] =
      std::unique_ptr<decltype(output_queue)::Block>(
        new decltype(output_queue)::Block(block_size));
    ready_blocks_owners()[id % MAX_SIMULTANEOUS_INDEXLRS] = id;
  }
  auto& block = *(ready_blocks_array()[id % MAX_SIMULTANEOUS_INDEXLRS]);
  if (block.count <= block.current) {
    output_queue.read(block);
    if (block.count <= block.current) {
      output_queue.close();
      block = decltype(output_queue)::Block(block_size);
      return Record();
    }
  }
  return std::move(block.data[block.current++]);
  */
  return;
}

//#endif