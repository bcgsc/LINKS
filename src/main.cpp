#include "InputParser.hpp"
#include "LINKS.hpp"
#include <chrono>
#include <stdlib.h> // for sleep

std::string format_duration(std::chrono::milliseconds ms) {
  using namespace std::chrono;
  auto secs = duration_cast<seconds>(ms);
  ms -= duration_cast<milliseconds>(secs);
  auto mins = duration_cast<minutes>(secs);
  secs -= duration_cast<seconds>(mins);
  auto hour = duration_cast<hours>(mins);
  mins -= duration_cast<minutes>(hour);

  std::stringstream ss;
  ss << hour.count() << "h " << mins.count() << "m " << secs.count() << "s "
     << ms.count() << "ms";
  return ss.str();
}

int main(int argc, char **argv) {
  InputParser linksArgParser(argc, argv);
  if (!linksArgParser.arguments_satisfied) {
    return -1;
  }
  linksArgParser.print_opts();
  std::cout << std::endl;
  LINKS links(linksArgParser);

  std::cout << "STAGE: Bloom filter constructing/reading" << std::flush;
  auto start = std::chrono::high_resolution_clock::now();
  links.init_bloom_filter();
  auto finish = std::chrono::high_resolution_clock::now();
  auto milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
  std::cout << "Bloom filter generated or loaded in: "
            << format_duration(milliseconds) << "\n\n";

  std::cout << "STAGE: Reading sequence reads" << std::flush;;
  start = std::chrono::high_resolution_clock::now();
  links.start_read_fasta();
  finish = std::chrono::high_resolution_clock::now();
  milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
  std::cout << "Read fasta stage completed in: "
            << format_duration(milliseconds) << "\n\n";

  std::cout << "STAGE: Reading contigs" << std::flush;;
  start = std::chrono::high_resolution_clock::now();
  links.start_read_contig();
  finish = std::chrono::high_resolution_clock::now();
  milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
  std::cout << "Read contig stage completed in: "
            << format_duration(milliseconds) << "\n\n";

  std::cout << "STAGE: Pairing contigs" << std::flush;;
  start = std::chrono::high_resolution_clock::now();
  links.pair_contigs();
  finish = std::chrono::high_resolution_clock::now();
  milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
  std::cout << "Pair contig stage completed in: "
            << format_duration(milliseconds) << "\n\n";

  return 0;
}
