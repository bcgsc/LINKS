#include "btllib/bloom_filter.hpp"
#include "btllib/nthash.hpp"
#include "btllib/order_queue.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/status.hpp"
#include "btllib/util.hpp"
#include "InputParser.hpp"
#include "Worker.hpp"

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

    /**
     * Construct a SeqReader to read sequences from a given path.
     *
     * @param seqfile Filepath to read sequences from. Pass "-" to read from
     * stdin.
     * @param k k-mer size for the minimizer.
     * @param w window size when selecting minimizers.
     * @param flags Modifier flags. Specifiying either short or long mode flag is
     * mandatory; other flags are optional.
     * @param threads Maximum number of processing threads to use. Must be at
     * least 1.
     * @param verbose Whether to output informational messages during processing.
     * @param bf1 A Bloom filter to use for either filtering minimizers in or out,
     * if one of the flags is enabled. If both are enabled, bf1 is used for
     * filtering in.
     * @param bf2 A second Bloom filter to use when both filtering minimizers in
     * and out is enabled. bf2 is used for filtering out.
     */

    LINKS(InputParser::InputParser inputParser);
        : assemblyFile(inputParser.assemblyFile)
        , fofFile(inputParser.fofFile)
        , distances(inputParser.distances)
        , k(inputParser.k)
        , verbose(inputParser.verbose)
        , minSize(inputParser.minSize)
        , step(inputParser.step)
        , insertStdev(inputParser.insertStdev)
        , baseName(inputParser.baseName)
        , offset(inputParser.offset)
        , fpr(inputParser.fpr)
        , bfFile(inputParser.bfFile)
        , bfOff(inputParser.bfOff)

    /*
    LINKS(  std::string seq,
            size_t k;
            unsigned threads = 5,
            bool verbose = false);
    */

    ~LINKS();

    InputParser::InputParser inputParser;
    uint64_t k;


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

private:
    std::unordered_map<uint64_t, std::unordered_map<uint64_t, MatePairInfo>> matePair;

    std::unordered_map<uint64_t, KmerInfo> trackAll;
    std::unordered_map<std::string, uint64_t> tigLength;

    // store second mates in a set
    std::unordered_set<uint64_t> mates;
   

    // Stage 1 -- populate bloom filter with contigs
    btllib::KmerBloomFilter * myFilter;
    
    std::atomic<bool> fasta{ false };
    /*
    
    static const size_t SHORT_MODE_BUFFER_SIZE = 32;
    static const size_t SHORT_MODE_BLOCK_SIZE = 32;

    static const size_t LONG_MODE_BUFFER_SIZE = 4;
    static const size_t LONG_MODE_BLOCK_SIZE = 1;
    */

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
  InputWorker::InputWorker input_worker;
  std::vector<ExtractMatePairWorker::ExtractMatePairWorker> extract_mate_pair_workers;
};

inline LINKS::LINKS(InputParser::InputParser inputParser)
        : assemblyFile(inputParser.assemblyFile)
        , fofFile(inputParser.fofFile)
        , distances(inputParser.distances)
        , k(inputParser.k)
        , verbose(inputParser.verbose)
        , minSize(inputParser.minSize)
        , step(inputParser.step)
        , insertStdev(inputParser.insertStdev)
        , baseName(inputParser.baseName)
        , offset(inputParser.offset)
        , fpr(inputParser.fpr)
        , bfFile(inputParser.bfFile)
        , bfOff(inputParser.bfOff){
            
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

inline LINKS::~Indexlr()
{
  reader.close();
  for (auto& worker : extract_mate_pair_workers) {
    worker.join();
  }
  input_worker.join();
}

// Minimerize a sequence: Find the minimizers of a vector of hash values
// representing a sequence.
/* Algorithm
v is a vector of non-negative integers
w is the window size
Invariants
    0 <  w <= v.size() - 1
    0 <= l <= r <= v.size() - 1
Initial conditions
    M    = NIL       Final set of minimizers, empty initially
    min  = -1        Minimum element
    i    = -1        Index of minimum element
    prev = -1        Index of previous minimum element
    l    = 0         Index of left end of window
    r    = l + w - 1 Index of right end of window
Computation
At each window, if the previous minimum is out of scope, find the new,
right-most, minimum or else, check with only the right-most element to determine
if that is the new minimum. A minimizer is added to the final vector only if
it's index has changed. for each window of v bounded by [l, r] if (i < l) i =
index of minimum element in [l, r], furthest from l. else if (v[r] <= v[i]) i =
r min = v[i] if (i != prev) { prev = i M <- M + m
    }
    l = l + 1        Move window's left bound by one element
    r = l + w - 1    Set window's right bound
}*/

inline std::string
Indexlr::extract_barcode(const std::string& id, const std::string& comment)
{
  const static std::string BARCODE_PREFIX = "BX:Z:";
  if (starts_with(comment, BARCODE_PREFIX)) {
    const auto space_pos = comment.find(' ');
    if (space_pos != std::string::npos) {
      return comment.substr(BARCODE_PREFIX.size(),
                            space_pos - BARCODE_PREFIX.size());
    }
    return comment.substr(BARCODE_PREFIX.size());
  }
  const auto pound_pos = id.find('#');
  if (pound_pos != std::string::npos) {
    const auto slash_pos = id.find('/');
    if (slash_pos > pound_pos) {
      return id.substr(pound_pos + 1, slash_pos - (pound_pos + 1));
    }
  }
  return "NA";
}

inline static void
filter_hashed_kmer(Indexlr::HashedKmer& hk,
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
calc_minimizer(const std::vector<Indexlr::HashedKmer>& hashed_kmers_buffer,
               const Indexlr::Minimizer*& min_current,
               const size_t idx,
               ssize_t& min_idx_left,
               ssize_t& min_idx_right,
               ssize_t& min_pos_prev,
               const size_t w,
               std::vector<Indexlr::Minimizer>& minimizers)
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

inline std::vector<Indexlr::Minimizer>
Indexlr::minimize(const std::string& seq) const
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

inline Indexlr::Record
Indexlr::get_minimizers()
{
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
}

inline void
Indexlr::InputWorker::work()
{
  if (indexlr.reader.get_format() == SeqReader::Format::FASTA) {
    indexlr.fasta = true;
  } else {
    indexlr.fasta = false;
  }

  decltype(indexlr.input_queue)::Block block(indexlr.block_size);
  size_t current_block_num = 0;
  SeqReader::Record record;
  Read read;
  while ((record = indexlr.reader.read())) {
    block.data[block.count++] = Read(record.num,
                                     std::move(record.id),
                                     std::move(record.comment),
                                     std::move(record.seq));
    if (block.count == indexlr.block_size) {
      block.num = current_block_num++;
      indexlr.input_queue.write(block);
      block.count = 0;
    }
  }
  if (block.count > 0) {
    block.num = current_block_num++;
    indexlr.input_queue.write(block);
  }
  for (unsigned i = 0; i < indexlr.threads; i++) {
    block.num = current_block_num++;
    block.current = 0;
    block.count = 0;
    indexlr.input_queue.write(block);
  }
}

inline void
Indexlr::MinimizeWorker::work()
{
  decltype(indexlr.input_queue)::Block input_block(indexlr.block_size);
  decltype(indexlr.output_queue)::Block output_block(indexlr.block_size);

  for (;;) {
    if (input_block.current == input_block.count) {
      if (output_block.count > 0) {
        output_block.num = input_block.num;
        indexlr.output_queue.write(output_block);
        output_block.current = 0;
        output_block.count = 0;
      }
      indexlr.input_queue.read(input_block);
    }
    if (input_block.count == 0) {
      output_block.num = input_block.num;
      output_block.current = 0;
      output_block.count = 0;
      indexlr.output_queue.write(output_block);
      break;
    }
    Read& read = input_block.data[input_block.current++];
    Record record;
    record.num = read.num;
    if (indexlr.output_id()) {
      record.id = std::move(read.id);
    }
    if (indexlr.output_bx()) {
      record.barcode = indexlr.extract_barcode(record.id, read.comment);
    }
    record.readlen = read.seq.size();

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

    if (indexlr.k <= read.seq.size() &&
        indexlr.w <= read.seq.size() - indexlr.k + 1) {
      record.minimizers = indexlr.minimize(read.seq);
    } else {
      record.minimizers = {};
    }

    output_block.data[output_block.count++] = std::move(record);
  }
}