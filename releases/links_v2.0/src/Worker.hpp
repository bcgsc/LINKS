
#ifndef LINKS_WORKER_HPP
#define LINKS_WORKER_HPP

/*
#include "bloom_filter.hpp"
#include "nthash.hpp"
#include "order_queue.hpp"
#include "seq_reader.hpp"
#include "status.hpp"
#include "util.hpp"
*/

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

#include "LINKS.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/seq_reader.hpp"

class Worker
  {
  public:
    void start() { t = std::thread(do_work, this); }
    void join() { t.join(); }

    virtual ~Worker() {}

    Worker& operator=(const Worker& worker) = delete;
    Worker& operator=(Worker&& worker) = delete;

  protected:
    Worker(LINKS& links)
      : links(links)
    {}

    Worker(const Worker& worker)
      : Worker(worker.links)
    {}
    Worker(Worker&& worker) noexcept
      : Worker(worker.links)
    {}

    LINKS& links;

    virtual void work() = 0;
    static void do_work(Worker* worker) { worker->work(); }

    std::thread t;
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

  btllib::SeqReader reader;
  InputWorker input_worker;
  std::vector<ExtractMatePairWorker> extract_mate_pair_workers;
};

inline void
InputWorker::work()
{
  if (links.reader.get_format() == btllib::SeqReader::Format::FASTA) {
    links.fasta = true;
  } else {
    links.fasta = false;
  }

  decltype(links.input_queue)::Block block(links.block_size);
  size_t current_block_num = 0;
  SeqReader::Record record;
  Read read;
  while ((record = links.reader.read())) {
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
ExtractMatePairWorker::work()
{
  decltype(links.input_queue)::Block input_block(links.block_size);
  decltype(links.output_queue)::Block output_block(links.block_size);

  for (;;) {
    if (input_block.current == input_block.count) {
      if (output_block.count > 0) {
        output_block.num = input_block.num;
        links.output_queue.write(output_block);
        output_block.current = 0;
        output_block.count = 0;
      }
      links.input_queue.read(input_block);
    }
    if (input_block.count == 0) {
      output_block.num = input_block.num;
      output_block.current = 0;
      output_block.count = 0;
      links.output_queue.write(output_block);
      break;
    }
    Read& read = input_block.data[input_block.current++];
    Record record;
    record.num = read.num;
    if (links.output_id()) {
      record.id = std::move(read.id);
    }
    if (links.output_bx()) {
      record.barcode = links.extract_barcode(record.id, read.comment);
    }
    record.readlen = read.seq.size();

    check_info(links.verbose && links.k > read.seq.size(),
               "Indexlr: skipped seq " + std::to_string(read.num) +
                 " on line " +
                 std::to_string(read.num * (links.fasta ? 2 : 4) + 2) +
                 "; k (" + std::to_string(links.k) + ") > seq length (" +
                 std::to_string(read.seq.size()) + ")");

    check_info(links.verbose && links.w > read.seq.size() - links.k + 1,
               "Indexlr: skipped seq " + std::to_string(read.num) +
                 " on line " +
                 std::to_string(read.num * (links.fasta ? 2 : 4) + 2) +
                 "; w (" + std::to_string(links.w) + ") > # of hashes (" +
                 std::to_string(read.seq.size() - links.k + 1) + ")");

    if (links.k <= read.seq.size() &&
        links.w <= read.seq.size() - links.k + 1) {
      record.minimizers = links.minimize(read.seq);
    } else {
      record.minimizers = {};
    }

    output_block.data[output_block.count++] = std::move(record);
  }
}

#endif