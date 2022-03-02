#ifndef LINKS_LINKS_HPP
#define LINKS_LINKS_HPP

#include "InputParser.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/nthash.hpp"
#include "btllib/order_queue.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/status.hpp"
#include "btllib/util.hpp"

#include <unordered_map>
#include <unordered_set>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

class LINKS {
public:
  /* ------------------- main functions ------------------- */
  LINKS(InputParser input_parser);
  void init_bloom_filter();
  void start_read_fasta();
  void start_read_contig();
  void pair_contigs();
  ~LINKS();
  /* ------------------- main functions ------------------- */
private:
  /* ------------------- structs ------------------- */
  struct BufferMatePairData {
    BufferMatePairData() {}

    BufferMatePairData(uint64_t kmer_1_hash, uint64_t kmer_2_hash,
                       uint32_t distance)
        : kmer_1_hash(kmer_1_hash), kmer_2_hash(kmer_2_hash),
          distance(distance) {}

    uint64_t kmer_1_hash;
    uint64_t kmer_2_hash;
    uint32_t distance;
  };

  struct BufferMateData {
    BufferMateData() {}

    BufferMateData(uint64_t hash, std::string tig, uint64_t start, uint64_t end,
                   uint64_t multiple, bool orient)
        : hash(hash), tig(tig), start(start), end(end), multiple(multiple),
          orient(orient) {}

    uint64_t hash;
    std::string tig = "";
    uint64_t start = 0;
    uint64_t end = 0;
    uint64_t multiple = 1;
    bool orient = 0;
  };

  struct Read {
    Read() {}

    Read(size_t num, std::string id, std::string comment, std::string seq)
        : num(num), id(std::move(id)), comment(std::move(comment)),
          seq(std::move(seq)) {}

    size_t num = 0;
    std::string id;
    std::string comment;
    std::string seq;
  };

  struct MatePairInfo {
    MatePairInfo() = default;

    MatePairInfo(bool seen, int64_t insert_size)
        : seen(seen), insert_size(insert_size) {}

    bool seen = false;
    int64_t insert_size = 0;
  };

  struct PairLinkInfo {
    PairLinkInfo() {}

    PairLinkInfo(int64_t gaps, uint64_t links, std::string type = "")
        : gaps(gaps), links(links), type(type) {}

    int64_t gaps = 0;
    uint64_t links = 0;
    std::string type;
  };

  struct KmerInfo {
    KmerInfo() {}

    KmerInfo(std::string tig, uint64_t start, uint64_t end, uint64_t multiple,
             bool orient)
        : tig(tig), start(start), end(end), multiple(multiple), orient(orient) {
    }

    std::string tig = "";
    uint64_t start = 0;
    uint64_t end = 0;
    uint64_t multiple = 1;
    bool orient = 0;
  };
  /* ------------------- structs ------------------- */

  /* xor operation overrides default hash function of map */
  struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
      return p.first ^ p.second;
    }
  };

  /* ------------------- worker classes ------------------- */
  class Worker {
  public:
    void start() { t = std::thread(do_work, this); }
    void join() { t.join(); }

    virtual ~Worker() {}

    Worker &operator=(const Worker &worker) = delete;
    Worker &operator=(Worker &&worker) = delete;

  protected:
    LINKS &links;
    std::thread t;

    Worker(LINKS &links) : links(links) {}

    Worker(const Worker &worker) : Worker(worker.links) {}
    Worker(Worker &&worker) noexcept : Worker(worker.links) {}

    virtual void work() = 0;
    static void do_work(Worker *worker) { worker->work(); }
  };

  class InputWorker : public Worker {
  public:
    InputWorker(LINKS &links) : Worker(links) {}

    InputWorker(const InputWorker &worker) : InputWorker(worker.links) {}
    InputWorker(InputWorker &&worker) noexcept : InputWorker(worker.links) {}

    InputWorker &operator=(const InputWorker &worker) = delete;
    InputWorker &operator=(InputWorker &&worker) = delete;

    void work() override;
  };

  class ExtractMatePairWorker : public Worker {
  public:
    ExtractMatePairWorker(LINKS &links) : Worker(links) {}

    ExtractMatePairWorker(const ExtractMatePairWorker &worker)
        : ExtractMatePairWorker(worker.links) {}
    ExtractMatePairWorker(ExtractMatePairWorker &&worker) noexcept
        : ExtractMatePairWorker(worker.links) {}

    ExtractMatePairWorker &
    operator=(const ExtractMatePairWorker &worker) = delete;
    ExtractMatePairWorker &operator=(ExtractMatePairWorker &&worker) = delete;

    void work() override;
  };

  class PopulateMateInfoWorker : public Worker {
  public:
    PopulateMateInfoWorker(LINKS &links) : Worker(links) {}

    PopulateMateInfoWorker(const PopulateMateInfoWorker &worker)
        : PopulateMateInfoWorker(worker.links) {}
    PopulateMateInfoWorker(PopulateMateInfoWorker &&worker) noexcept
        : PopulateMateInfoWorker(worker.links) {}

    PopulateMateInfoWorker &
    operator=(const PopulateMateInfoWorker &worker) = delete;
    PopulateMateInfoWorker &operator=(PopulateMateInfoWorker &&worker) = delete;

    void work() override;
  };
  /* ------------------- worker classes ------------------- */

  /* ------------------- helper functions ------------------- */
  btllib::KmerBloomFilter *make_bf(uint64_t bf_elements);
  void extract_mate_pair(
      const std::string &seq,
      btllib::OrderQueueSPMC<BufferMatePairData>::Block &mate_pair_block);
  void
  populate_mate_info(const std::string &seq, const std::string contig_rank,
                     btllib::OrderQueueSPMC<BufferMateData>::Block &mate_block);
  void add_to_pair_map(
      int isz,
      std::unordered_map<
          std::string,
          std::unordered_map<
              int64_t, std::unordered_map<std::string, PairLinkInfo>>> &pair,
      int distance, std::string kmer1_name, std::string kmer2_name,
      unsigned orient_enum);
  void write_from_block_to_map();
  void write_from_block_to_set();

  int get_distance_bin(int distance);
  int get_distance(uint64_t insert_size, uint64_t length_i, uint64_t start_i,
                   uint64_t start_j);
  uint64_t get_file_size(std::string file_name);
  bool does_file_exist(std::string file_name);
  /* ------------------- helper functions ------------------- */

  /* ------------------- constants ------------------- */
  static const size_t PRINT_READ_COUNT_PERIOD = 100000;
  static const size_t READ_BUFFER_SIZE = 16;
  static const size_t READ_BLOCK_SIZE = 4;
  static const size_t MATE_PAIR_BUFFER_SIZE = 6;
  static const size_t MATE_PAIR_BLOCK_SIZE = 1000000;
  /* ------------------- constants ------------------- */

  /* ------------------- main fields called from various functions ------------------- */
  InputParser input_parser;
  
  std::shared_ptr<btllib::KmerBloomFilter> bloom;
  std::shared_ptr<btllib::SeqReader> reader;

  std::shared_ptr<InputWorker> input_worker;
  std::shared_ptr<btllib::OrderQueueSPMC<Read>> input_queue;
  btllib::OrderQueueSPMC<BufferMatePairData> mate_pair_input_queue;
  btllib::OrderQueueSPMC<BufferMateData> mate_input_queue;
  
  std::vector<ExtractMatePairWorker> extract_mate_pair_workers;
  std::vector<PopulateMateInfoWorker> populate_mate_info_workers;

  typedef std::unordered_map<std::pair<uint64_t, uint64_t>, MatePairInfo,
                          pair_hash>
  mate_pair_type;
  mate_pair_type mate_pair;
  // stage coordination variable. false when reading reads(stage 1), true when reading contigs(stage 2)
  bool reading_contig = false;
  std::atomic<bool> fasta{false};

  // track information of kmers
  std::unordered_map<uint64_t, KmerInfo> track_all;
  // contig length map
  std::unordered_map<std::string, uint64_t> tig_length;
  // store second mates in a set
  std::unordered_set<uint64_t> mates;
  /* ------------------- main fields called from various functions ------------------- */

  /* ------------------- thread coordination variables ------------------- */
  std::mutex tig_length_mutex;
  std::atomic<std::size_t> mate_pair_current_block_num = {0};
  std::atomic<std::size_t> mate_pair_threads_done_writing = {0};
  std::atomic<std::size_t> mate_current_block_num = {0};
  std::atomic<std::size_t> mate_threads_done_writing = {0};
  /* ------------------- thread coordination variables ------------------- */

  /* ------------------- verbose counters ------------------- */
  uint ct_ok_pairs = 0;      // last stage verbose
  uint ct_problem_pairs = 0; // last stage verbose
  /* ------------------- verbose counters ------------------- */
};

inline btllib::KmerBloomFilter *LINKS::make_bf(uint64_t bf_elements) {
  btllib::KmerBloomFilter *assembly_BF;
  if (input_parser.bf_file != "") {
    std::cout << "A Bloom filter was supplied (" << input_parser.bf_file
              << ") and will be used instead of building a new one\n";
    if (!does_file_exist(input_parser.bf_file)) {
      std::cout << "\nInvalid file: " << input_parser.bf_file
                << " -- fatal\n";
      exit(1);
    }

    assembly_BF = new btllib::KmerBloomFilter(input_parser.bf_file);
  } else {
    uint64_t m = ceil((-1 * (double)bf_elements * log(input_parser.fpr)) /
                      (log(2) * log(2)));

    m = ((uint64_t)(m / 8) + 1) * 8;

    unsigned hash_fct = floor(((double)m / bf_elements) * log(2));
    std::cout << "- Number of BF Elements: " << bf_elements << "\n"
              << "- Input file path: " << input_parser.base_name << "\n"
              << "- Input file: " << input_parser.assembly_file << "\n"
              << "- K: " << input_parser.k << "\n"
              << "- Fpr: " << input_parser.fpr << "\n"
              << "- #Hash functions: " << hash_fct << "\n";

    std::string reading_tigbloom_message =
        "\n\n=>Reading contig/sequence assembly file : " +
        std::to_string(time(0)) + "\n";

    std::cout << "- Filter output file : " << input_parser.k << "\n";
    assembly_BF =
        new btllib::KmerBloomFilter(m / 8, hash_fct, input_parser.k);
    btllib::SeqReader assembly_reader(input_parser.assembly_file, 8, 1);
    size_t builder = 0;
    for (btllib::SeqReader::Record record; (record = assembly_reader.read());) {

      builder++;
      assembly_BF->insert(record.seq);
    }
    std::string bfmsg = "\n\nWriting Bloom filter to disk (" +
                        input_parser.base_name + ".bloom" + ")\n";
    std::cout << bfmsg;
    assembly_BF->save(input_parser.base_name + ".bloom");
  }
  return assembly_BF;
}

inline LINKS::LINKS(InputParser input_parser)
    : input_parser(input_parser),
      input_queue(std::shared_ptr<btllib::OrderQueueSPMC<Read>>(
          new btllib::OrderQueueSPMC<Read>(READ_BUFFER_SIZE, READ_BLOCK_SIZE))),
      mate_pair_input_queue(MATE_PAIR_BUFFER_SIZE, MATE_PAIR_BLOCK_SIZE),
      mate_input_queue(MATE_PAIR_BUFFER_SIZE, MATE_PAIR_BLOCK_SIZE),
      input_worker(std::shared_ptr<InputWorker>(new InputWorker(*this))) {}

inline void LINKS::init_bloom_filter() {
  int64_t bf_elements = get_file_size(input_parser.assembly_file);
  bloom = std::shared_ptr<btllib::KmerBloomFilter>(
      make_bf(bf_elements));
}

inline void LINKS::write_from_block_to_map() {
  btllib::OrderQueueSPMC<BufferMatePairData>::Block mate_pair_block(
      MATE_PAIR_BLOCK_SIZE);

  for (;;) {
    if (mate_pair_block.current == mate_pair_block.count) {
      mate_pair_input_queue.read(mate_pair_block);
    }

    if (mate_pair_block.count == 0) {
      break;
    }

    BufferMatePairData &mate_data =
        mate_pair_block.data[mate_pair_block.current++];

    if (mate_pair.find(std::make_pair(
            mate_data.kmer_1_hash, mate_data.kmer_2_hash)) == mate_pair.end()) {
      mate_pair[std::make_pair(mate_data.kmer_1_hash, mate_data.kmer_2_hash)] =
          MatePairInfo(false, mate_data.distance);
      mates.insert(mate_data.kmer_1_hash);
      mates.insert(mate_data.kmer_2_hash);
    }
  }
}
inline void LINKS::write_from_block_to_set() {
  btllib::OrderQueueSPMC<BufferMateData>::Block mate_block(
      MATE_PAIR_BLOCK_SIZE);

  for (;;) {
    if (mate_block.current == mate_block.count) {
      mate_input_queue.read(mate_block);
    }

    if (mate_block.count == 0) {
      break;
    }

    BufferMateData &mate_data = mate_block.data[mate_block.current++];

    if (track_all.find(mate_data.hash) == track_all.end()) {
      track_all[mate_data.hash] =
          KmerInfo(mate_data.tig, mate_data.start, mate_data.end,
                   mate_data.multiple, mate_data.orient);
    } else {
      track_all[mate_data.hash].multiple += 1;
    }
  }
}

inline void LINKS::InputWorker::work() {
  btllib::OrderQueueSPMC<Read>::Block block(links.READ_BLOCK_SIZE);
  size_t current_block_num = 0;
  size_t read_counter = 0;
  std::string read_file;
  Read read;

  while(links.reading_contig || links.input_parser.long_reads.size() > 0) {
    if(!links.reading_contig){
      read_file = links.input_parser.long_reads.front();
      links.input_parser.long_reads.pop();
      links.reader = std::shared_ptr<btllib::SeqReader>(
        new btllib::SeqReader(read_file, btllib::SeqReader::Flag::LONG_MODE));
    } else {
        links.reader = std::shared_ptr<btllib::SeqReader>(
            new btllib::SeqReader(links.input_parser.assembly_file, btllib::SeqReader::Flag::LONG_MODE));
    }

    if (links.reader->get_format() == btllib::SeqReader::Format::FASTA) {
      links.fasta = true;
    } else {
      links.fasta = false;
    }

    std::cout << "Reading: " << read_file << std::endl;

    for (auto record : (*(links.reader))) {
      if(links.reading_contig && record.seq.length() < links.input_parser.min_size){
        continue;
      }
      block.data[block.count++] =
          Read(record.num, std::move(record.id), std::move(record.comment),
               std::move(record.seq));
      if (block.count == links.READ_BLOCK_SIZE) {
        block.num = current_block_num++;
        links.input_queue->write(block);
        block.count = 0;
      }
      read_counter++;
      if (read_counter % PRINT_READ_COUNT_PERIOD == 0) {
        if (read_counter == PRINT_READ_COUNT_PERIOD) {
          std::cout << "Processed read count: " << read_counter;
        } else {
          std::cout << " - " << read_counter << std::flush; // avoid buffering of stdout
        }
      }
    }
    if(links.reading_contig){
      links.reading_contig = false;
    }
  }

  if (block.count > 0) {
    block.num = current_block_num++;
    links.input_queue->write(block);
  }
  for (unsigned i = 0; i < links.input_parser.threads; i++) {
    block.num = current_block_num++;
    block.current = 0;
    block.count = 0;
    links.input_queue->write(block);
  }
  if (read_counter < PRINT_READ_COUNT_PERIOD) {
    std::cout << "Processed read count: " << read_counter << std::endl;
  } else {
    std::cout << " - " << read_counter << std::endl; // final 'processed reads' print
  }
}

inline void LINKS::start_read_fasta() {
  extract_mate_pair_workers =
      std::vector<ExtractMatePairWorker>(input_parser.threads, ExtractMatePairWorker(*this));

  input_worker->start();
  std::thread writer_thread(&LINKS::write_from_block_to_map, this);

  for (auto &worker : extract_mate_pair_workers) {
    worker.start();
  }
  // wait
  for (auto &worker : extract_mate_pair_workers) {
    worker.join();
  }
  writer_thread.join();
  input_worker->join();
}
inline void LINKS::start_read_contig() {
  reading_contig = true;

  input_worker.reset();
  input_worker = std::shared_ptr<InputWorker>(new InputWorker(*this));

  input_queue.reset();
  input_queue = std::shared_ptr<btllib::OrderQueueSPMC<Read>>(
      new btllib::OrderQueueSPMC<Read>(READ_BUFFER_SIZE, READ_BLOCK_SIZE));

  input_worker->start();

  std::thread writer_thread(&LINKS::write_from_block_to_set, this);

  populate_mate_info_workers = std::vector<PopulateMateInfoWorker>(
      input_parser.threads, PopulateMateInfoWorker(*this));

  for (auto &worker : populate_mate_info_workers) {
    worker.start();
  }
  for (auto &worker : populate_mate_info_workers) {
    worker.join();
  }
  writer_thread.join();
  input_worker->join();
  input_queue.reset();
}

inline void LINKS::extract_mate_pair(
    const std::string &seq,
    btllib::OrderQueueSPMC<BufferMatePairData>::Block &mate_pair_block) {
  uint step_index = 0;
  for (auto distance : input_parser.distances) {
    uint cur_step_size = input_parser.step_sizes[step_index];
    step_index++;
    uint64_t delta = distance;
    int break_flag = 0;
    bool reverse_exists = false;
    btllib::NtHash ntHash(seq, bloom->get_hash_num(), input_parser.k, input_parser.offset);
    btllib::NtHash ntHash_lead(seq, bloom->get_hash_num(), input_parser.k, delta + input_parser.offset);

    uint64_t hash_a, hash_b;

    for (size_t i = 0; ntHash.roll() && ntHash_lead.roll();
         i += cur_step_size) {
      // roll for the number of steps
      break_flag = 0;
      reverse_exists = false;
      // for step ----
      for (uint j = 1; j < cur_step_size; j++) {
        if (!ntHash_lead.roll() || !ntHash.roll()) {
          break_flag = 1;
        }
      }
      if (break_flag) {
        break;
      }
      // for step ----

      if (ntHash.get_forward_hash() < ntHash_lead.get_reverse_hash()) {
        hash_a = ntHash_lead.get_reverse_hash();
        hash_b = ntHash.get_reverse_hash();
      } else {
        hash_a = ntHash.get_forward_hash();
        hash_b = ntHash_lead.get_forward_hash();
      }

      if (bloom->contains(ntHash.hashes()) &&
          bloom->contains(ntHash_lead.hashes())) { // May need to change with
                                                   // forward reverse hashes
        mate_pair_block.data[mate_pair_block.count++] =
            BufferMatePairData(hash_a, hash_b, delta);
        if (mate_pair_block.count == MATE_PAIR_BLOCK_SIZE) {
          mate_pair_block.num = mate_pair_current_block_num++;
          mate_pair_input_queue.write(mate_pair_block);
          mate_pair_block.count = 0;
        }
      }
    }
  }
}

inline void LINKS::populate_mate_info(
    const std::string &seq, const std::string contig_rank,
    btllib::OrderQueueSPMC<BufferMateData>::Block &mate_info_block) {

  btllib::NtHash ntHash_contig(seq, bloom->get_hash_num(),
                               input_parser.k); // hashFunc can be 1 after first step

  int breakFlag = 0;
  for (size_t i = 0; ntHash_contig.roll(); i += 1) {

    i = ntHash_contig.get_pos();
    if (mates.find(ntHash_contig.get_forward_hash()) != mates.end()) {
      mate_info_block.data[mate_info_block.count++] = BufferMateData(
          ntHash_contig.get_forward_hash(), contig_rank, i, i + input_parser.k, 1, false);
      if (mate_info_block.count == MATE_PAIR_BLOCK_SIZE) {
        mate_info_block.num = mate_current_block_num++;
        mate_input_queue.write(mate_info_block);
        mate_info_block.count = 0;
      }
    }

    // Reverse part
    if (mates.find(ntHash_contig.get_reverse_hash()) != mates.end()) {
      mate_info_block.data[mate_info_block.count++] = BufferMateData(
          ntHash_contig.get_reverse_hash(), contig_rank, i, i + input_parser.k, 1, true);
      if (mate_info_block.count == MATE_PAIR_BLOCK_SIZE) {
        mate_info_block.num = mate_current_block_num++;
        mate_input_queue.write(mate_info_block);
        mate_info_block.count = 0;
      }
    }
  }
}
inline void LINKS::ExtractMatePairWorker::work() {
  btllib::OrderQueueSPMC<Read>::Block input_block(links.READ_BLOCK_SIZE);

  btllib::OrderQueueSPMC<BufferMatePairData>::Block mate_pair_block(
      links.MATE_PAIR_BLOCK_SIZE);

  for (;;) {
    if (input_block.current == input_block.count) {
      links.input_queue->read(input_block);
    }

    if (input_block.count == 0) {
      break;
    }

    Read &read = input_block.data[input_block.current++];

    if (links.input_parser.k <= read.seq.size()) {
      links.extract_mate_pair(read.seq, mate_pair_block);
    } else {
      continue; // nothing
    }
  }
  if (mate_pair_block.count > 0) {
    mate_pair_block.num = links.mate_pair_current_block_num++;
    links.mate_pair_input_queue.write(mate_pair_block);
  }
  links.mate_pair_threads_done_writing++;
  if (links.mate_pair_threads_done_writing == links.input_parser.threads) {
    for (size_t i = 0; i < links.input_parser.threads; i++) {
      mate_pair_block.num = links.mate_pair_current_block_num++;
      mate_pair_block.current = 0;
      mate_pair_block.count = 0;
      links.mate_pair_input_queue.write(mate_pair_block);
    }
  }
}

inline void LINKS::PopulateMateInfoWorker::work() {
  btllib::OrderQueueSPMC<Read>::Block input_block(links.READ_BLOCK_SIZE);
  btllib::OrderQueueSPMC<BufferMateData>::Block mate_block(
      links.MATE_PAIR_BLOCK_SIZE);
  std::unordered_map<uint64_t, KmerInfo> own_track_all;

  for (;;) {
    if (input_block.current == input_block.count) {
      links.input_queue->read(input_block);
    }
    if (input_block.count == 0) {
      break;
    }
    Read &read = input_block.data[input_block.current++];

    if (links.input_parser.k <= read.seq.size()) {
      links.tig_length_mutex.lock();
      links.tig_length[std::to_string(read.num + 1)] = read.seq.size();
      links.tig_length_mutex.unlock();
      links.populate_mate_info(read.seq, std::to_string(read.num + 1),
                               mate_block);
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
  if (links.mate_threads_done_writing == links.input_parser.threads) {

    for (size_t i = 0; i < links.input_parser.threads; i++) {
      mate_block.num = links.mate_current_block_num++;
      mate_block.current = 0;
      mate_block.count = 0;
      links.mate_input_queue.write(mate_block);
    }
  }
}

inline LINKS::~LINKS() {}

inline uint64_t LINKS::get_file_size(std::string file_name) {
  // This buffer is a stat struct that the information is placed concerning the
  // file.
  struct stat stat_buf;
  // Stat method returns true if successfully completed
  int rc = stat(file_name.c_str(), &stat_buf);
  // st_size holds the total size of the file in bytes
  return rc == 0 ? stat_buf.st_size : -1;
}

inline bool LINKS::does_file_exist(std::string file_name) {
  std::ifstream infile(file_name);
  return infile.good();
}

inline void LINKS::add_to_pair_map(
    int isz,
    std::unordered_map<
        std::string,
        std::unordered_map<
            int64_t, std::unordered_map<std::string, PairLinkInfo>>> &pair,
    int distance, std::string kmer1_name, std::string kmer2_name,
    unsigned orient_enum) {
  ct_ok_pairs++;
  std::tuple<std::string, std::string> first_pair;
  std::tuple<std::string, std::string> second_pair;

  std::string ftig_a = "f" + kmer1_name;
  std::string ftig_b = "f" + kmer2_name;

  std::string rtig_a = "r" + kmer1_name;
  std::string rtig_b = "r" + kmer2_name;

  switch (orient_enum) {
  case 0: // A B or rB rA
    first_pair = std::make_tuple(ftig_a, ftig_b);
    second_pair = std::make_tuple(rtig_b, rtig_a);
    break;
  case 1: // A B or A rB
    first_pair = std::make_tuple(ftig_a, rtig_b);
    second_pair = std::make_tuple(ftig_b, rtig_a);
    break;
  case 2: // rA B or rB A
    first_pair = std::make_tuple(rtig_a, ftig_b);
    second_pair = std::make_tuple(rtig_b, ftig_a);
    break;
  case 3: // rA rB or B A
    first_pair = std::make_tuple(rtig_a, rtig_b);
    second_pair = std::make_tuple(ftig_b, ftig_a);
    break;
  default:
    break;
  }

  if (pair.find(std::get<0>(first_pair)) == pair.end() ||
      pair[std::get<0>(first_pair)].find(isz) ==
          pair[std::get<0>(first_pair)].end() ||
      pair[std::get<0>(first_pair)][isz].find(std::get<1>(first_pair)) ==
          pair[std::get<0>(first_pair)][isz].end()) {
    pair[std::get<0>(first_pair)][isz][std::get<1>(first_pair)] =
        PairLinkInfo(distance, 1);
  } else {
    pair[std::get<0>(first_pair)][isz][std::get<1>(first_pair)].gaps +=
        distance;
    pair[std::get<0>(first_pair)][isz][std::get<1>(first_pair)].links += 1;
  }
  if (pair.find(std::get<0>(second_pair)) == pair.end() ||
      pair[std::get<0>(second_pair)].find(isz) ==
          pair[std::get<0>(second_pair)].end() ||
      pair[std::get<0>(second_pair)][isz].find(std::get<1>(second_pair)) ==
          pair[std::get<0>(second_pair)][isz].end()) {
    pair[std::get<0>(second_pair)][isz][std::get<1>(second_pair)] =
        PairLinkInfo(distance, 1);
  } else {
    pair[std::get<0>(second_pair)][isz][std::get<1>(second_pair)].gaps +=
        distance;
    pair[std::get<0>(second_pair)][isz][std::get<1>(second_pair)].links += 1;
  }
}

inline int LINKS::get_distance_bin(int distance) {
  return distance < 0      ? -1
         : distance == 10  ? 10
         : distance < 500  ? 500
         : distance < 5000 ? 5000
                           : 10000;
}

inline int LINKS::get_distance(uint64_t insert_size, uint64_t length_i,
                               uint64_t start_i, uint64_t start_j) {

  int insert_span = (length_i - start_i) + start_j;
  int gap_or_overlap = insert_size - insert_span;

  return gap_or_overlap;
}
inline void LINKS::pair_contigs() {
  std::string issues = input_parser.base_name + ".pairing_issues";
  std::string distribution =
      input_parser.base_name + ".pairing_distribution.csv";
  std::string tigpair_checkpoint =
      input_parser.base_name +
      ".tigpair_checkpoint.tsv"; // add a checkpoint file, prevent re-running
                                 // LINKS from scratch if crash
  std::string simplepair_checkpoint =
      input_parser.base_name +
      ".simplepair_checkpoint.tsv"; // add a checkpoint file, prevent re-running
                                    // LINKS from scratch if cras
  uint64_t totalPairs = 0;
  uint64_t ct_illogical = 0, ct_ok_contig = 0, ct_problem_pairs = 0,
           ct_iz_issues = 0, ct_single = 0, ct_multiple = 0, ct_both = 0,
           track_insert = 0;
  std::unordered_map<uint64_t, uint64_t> ct_both_hash, ct_illogical_hash,
      ct_ok_contig_hash, ct_ok_pairs_hash, ct_problem_pairs_hash,
      ct_iz_issues_hash;
  // Mapping of tiga_head -> insertSize -> tigb_head -> links & gaps
  std::unordered_map<
      std::string, std::unordered_map<
                       int64_t, std::unordered_map<std::string, PairLinkInfo>>>
      pair;
  std::unordered_map<std::string, PairLinkInfo> err;

  if (input_parser.verbose)
    std::cout << "Pairing contigs...\n";
  std::ofstream issues_file;
  issues_file.open(issues);
  int64_t insert_size = 0;
  int min_allowed = 0;
  uint low_iz = 0;
  uint up_iz = 0;
  int distance;
  int isz;

  std::string tig_a, tig_b, ftig_a, ftig_b, rtig_a, rtig_b;

  KmerInfo kmer1, kmer2;

  size_t counter = 0;
  size_t percent_size = mate_pair.size() / 25;
  mate_pair_type::iterator mate_pair_iterator;
  for (mate_pair_iterator = mate_pair.begin();
       mate_pair_iterator != mate_pair.end(); mate_pair_iterator++) {
    ++counter;
    if (mate_pair_iterator->second.seen == false && // matepair is not seen
        track_all[mate_pair_iterator->first.first].multiple ==
            1 && // first mate seen once
        track_all[mate_pair_iterator->first.second].multiple == 1) {

      mate_pair_iterator->second.seen = true;

      insert_size = mate_pair[std::make_pair(mate_pair_iterator->first.first,
                                             mate_pair_iterator->first.second)]
                        .insert_size;

      min_allowed = -1 * (input_parser.insert_stdev * insert_size); // check int
      low_iz = insert_size + min_allowed;              // check int
      up_iz = insert_size - min_allowed;               // check int

      if (track_all.find(mate_pair_iterator->first.first) !=
              track_all.end() && // first mate is tracked
          track_all.find(mate_pair_iterator->first.second) !=
              track_all.end() // second mate is tracked
      ) {

        kmer1 = track_all[mate_pair_iterator->first.first];
        kmer2 = track_all[mate_pair_iterator->first.second];

        if (kmer1.tig != kmer2.tig) { // paired reads located on <> contigs

          if (!kmer1.orient) {   // if kmer1 is forward
            if (!kmer2.orient) { // if kmer1 is forward
              // if kmer2 is forward
              distance = get_distance(insert_size, tig_length[kmer1.tig],
                                      kmer1.start, kmer2.start);
              if (distance > min_allowed && distance < insert_size) {
                isz = get_distance_bin(distance);
                add_to_pair_map(isz, pair, distance, kmer1.tig, kmer2.tig, 0);
              } else {
                ct_problem_pairs++;
              }
            } else { // if kmer1 is forward
              // if kmer2 is reverse
              distance =
                  get_distance(insert_size, tig_length[kmer1.tig], kmer1.start,
                               tig_length[kmer2.tig] - kmer2.end);
              if (distance > min_allowed && distance < insert_size) {
                isz = get_distance_bin(distance);
                add_to_pair_map(isz, pair, distance, kmer1.tig, kmer2.tig, 1);
              } else {
                ct_problem_pairs++;
              }
            }
          } else { // if kmer1 is reverse
            // if kmer2 is forward
            if (!kmer2.orient) {
              distance =
                  get_distance(insert_size, tig_length[kmer1.tig],
                               tig_length[kmer1.tig] - kmer1.end, kmer2.start);
              if (distance > min_allowed && distance < insert_size) {
                isz = get_distance_bin(distance);
                add_to_pair_map(isz, pair, distance, kmer1.tig, kmer2.tig, 2);
              } else {
                ct_problem_pairs++;
              }
            } else { // if kmer1 is reverse
              // if kmer2 is reverse
              distance = get_distance(insert_size, tig_length[kmer2.tig],
                                      kmer2.end, kmer1.end);
              if (distance > min_allowed && distance < insert_size) {
                isz = get_distance_bin(distance);
                add_to_pair_map(isz, pair, distance, kmer1.tig, kmer2.tig, 3);
              } else {
                ct_problem_pairs++;
              }
            }
          }
        } else { // Clone, paired reads located on the same contig -- could be
                 // used to investigate misassemblies
          if (input_parser.verbose)
            std::cout << "Pair (" << mate_pair_iterator->first.first << " and "
                      << mate_pair_iterator->first.second
                      << ") located on same contig " << tig_a << " ("
                      << tig_length[kmer1.tig] << " nt)\n";
          uint64_t pet_size = 0;
          if (kmer1.orient == kmer2.orient) { // if kmer1 is forward
            if (kmer1.start > kmer2.start) {  // if kmer2 is forward
                                              // A-> B->
              pet_size = kmer1.start - kmer2.start;
              if (pet_size >= low_iz && pet_size <= up_iz) {
                ct_ok_contig++;
              } else {
                issues_file
                    << "Pairs unsatisfied in distance within a contig.  Pairs "
                       "on "
                    << "contig id: " << kmer1.tig
                    << " (length: " << tig_length[kmer1.tig] << ")\n"
                    << "kmer start: " << kmer1.start
                    << " kmer end: " << kmer1.end
                    << " orientation: " << kmer1.orient << "\n"
                    << "kmer start: " << kmer2.start
                    << " kmer end: " << kmer2.end
                    << " orientation: " << kmer2.orient << "\n"
                    << "Calculated distance apart on contig: " << pet_size
                    << " - Distance on read: " << insert_size << "\n";
                ct_iz_issues++;
              }
            } else {
              pet_size = kmer2.start - kmer1.start;
              if (pet_size >= low_iz && pet_size <= up_iz) {
                ct_ok_contig++;
              } else {
                issues_file
                    << "Pairs unsatisfied in distance within a contig.  Pairs "
                       "on "
                    << "contig id: " << kmer1.tig
                    << " (length: " << tig_length[kmer1.tig] << ")\n"
                    << "kmer start: " << kmer1.start
                    << " kmer end: " << kmer1.end
                    << " orientation: " << kmer1.orient << "\n"
                    << "kmer start: " << kmer2.start
                    << " kmer end: " << kmer2.end
                    << " orientation: " << kmer2.orient << "\n"
                    << "Calculated distance apart on contig: " << pet_size
                    << " - Distance on read: " << insert_size << "\n";
                ct_iz_issues++;
              }
            }
          } else {
            ct_illogical++;
          }
        }
      }
    } else if (track_all.find(mate_pair_iterator->first.first) ==
                   track_all.end() ^
               track_all.find(mate_pair_iterator->first.second) ==
                   track_all.end() // one of the kmers is not assembled
    ) { // either one of the kmers is not assembled
      ct_single++;
    } else if (track_all.find(mate_pair_iterator->first.first) !=
                   track_all.end() &&
               track_all.find(mate_pair_iterator->first.second) !=
                   track_all.end()) {
      ct_multiple++;
    }
  } // pairing read b
  // pairing read a
  std::cout << "All contig pairs traversed\n";
  issues_file.close();

  std::unordered_map<std::string, PairLinkInfo>::iterator err_itr;
  for (err_itr = err.begin(); err_itr != err.end(); err_itr++) {
    double mean_iz = 0;
    if (err_itr->second.links) {
      mean_iz = err_itr->second.gaps / err_itr->second.links;
    }
  }

  uint64_t satisfied = ct_ok_pairs + ct_ok_contig;
  uint64_t unsatisfied = ct_problem_pairs + ct_iz_issues + ct_illogical;
  uint64_t ct_both_reads = ct_both * 2;

  std::cout << "\n===========PAIRED K-MER STATS===========\n";
  std::cout << "At least one sequence/pair missing from contigs: " << ct_single
            << "\n";
  std::cout << "Ambiguous kmer pairs (both kmers are ambiguous): "
            << ct_multiple << "\n";
  std::cout << "Assembled pairs: " << ct_both << " (" << ct_both_reads
            << " sequences)\n";
  std::cout << "\tSatisfied in distance/logic within contigs (i.e. -> <-, "
               "distance on target): "
            << ct_ok_contig << "\n";
  std::cout << "\tUnsatisfied in distance within contigs (i.e. distance "
               "out-of-bounds): "
            << ct_iz_issues << "\n";
  std::cout << "\tUnsatisfied pairing logic within contigs (i.e. illogical "
               "pairing ->->, <-<- or <-->): "
            << ct_illogical << "\n";
  std::cout << "\t---\n";
  std::cout << "\tSatisfied in distance/logic within a given contig pair "
               "(pre-scaffold): "
            << ct_ok_pairs << "\n";
  std::cout << "\tUnsatisfied in distance within a given contig pair (i.e. "
               "calculated distances out-of-bounds): "
            << ct_problem_pairs << "\n";
  std::cout << "\t---\n";
  std::cout << "Total satisfied: " << satisfied
            << "\tunsatisfied: " << unsatisfied;

  std::cout << "\n============================================\n\n";

  // TIGPAIR CHECKPOINT
  std::ofstream tigpair_checkpoint_file;
  tigpair_checkpoint_file.open(tigpair_checkpoint);
  std::unordered_map<
      std::string,
      std::unordered_map<
          int64_t, std::unordered_map<std::string, PairLinkInfo>>>::iterator
      pair_itr;

  for (pair_itr = pair.begin(); pair_itr != pair.end(); pair_itr++) {
    std::unordered_map<int64_t,
                       std::unordered_map<std::string, PairLinkInfo>>::iterator
        insert_sizes;
    for (insert_sizes = pair_itr->second.begin();
         insert_sizes != pair_itr->second.end(); insert_sizes++) {
      std::unordered_map<std::string, PairLinkInfo>::iterator sec_pair_itr;
      for (sec_pair_itr = insert_sizes->second.begin();
           sec_pair_itr != insert_sizes->second.end(); sec_pair_itr++) {
        tigpair_checkpoint_file
            << insert_sizes->first /*distance*/ << "\t" << pair_itr->first
            << "\t" << sec_pair_itr->first << "\t" << sec_pair_itr->second.links
            << "\t" << sec_pair_itr->second.gaps << "\n";
      }
    }
  }
  tigpair_checkpoint_file.close();
}

#endif
