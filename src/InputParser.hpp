#ifndef LINKS_INPUTPARSER_HPP
#define LINKS_INPUTPARSER_HPP

#include <fstream>
#include <getopt.h>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

// Globals
#define BASE_TEN 10
static const std::string progname = "LINKS";
static const std::string version = "2.0.1";

class InputParser {
private:
  std::vector<std::string> tokens;
  int argc;
  int c;
  int optindex = 0;
  int help = 0;
  char *end = nullptr;

  std::vector<uint32_t> split_distance_input(std::string input) {
    std::vector<uint32_t> distances;

    std::stringstream ss(input);

    while (ss.good()) {
      std::string substr;
      getline(ss, substr, ',');
      distances.push_back(static_cast<uint32_t>(std::stoul(substr)));
    }
    return distances;
  }

  std::vector<uint16_t> split_step_sizes_input(std::string input) {
    std::vector<uint16_t> step_sizes;

    std::stringstream ss(input);

    while (ss.good()) {
      std::string substr;
      getline(ss, substr, ',');
      step_sizes.push_back(static_cast<uint16_t>(std::stoul(substr)));
    }
    return step_sizes;
  }

  std::queue<std::string> parse_fof_input(std::string fof_file) {
    std::queue<std::string> long_reads;

    std::ifstream infile(fof_file);
    std::string read_file;
    while (infile >> read_file) {
      long_reads.push(read_file);
    }

    return long_reads;
  }

public:
  std::string assembly_file = "";
  std::string fof_file = "";
  bool arguments_satisfied = true;
  std::queue<std::string> long_reads;
  std::string distances_text = "";
  std::vector<uint32_t> distances = {4000};
  std::string step_sizes_text = "2";
  std::vector<uint16_t> step_sizes = {2};
  uint64_t k = 15;
  bool verbose = false;
  uint64_t min_links = 5;
  uint64_t min_size = 500;
  float max_link_ratio = 0.3;
  float insert_stdev = 0.1;
  std::string base_name = "";
  uint64_t offset = 0;
  std::string bf_file = "";
  float fpr = 0.001;
  bool bf_off = false;
  uint thread = 3;

  static void print_usage() {
    std::cout
        << "Usage:  " << progname << " " << version << "\n"
        << "  -f  sequences to scaffold (Multi-FASTA format, required)\n"
           "  -s  file-of-filenames, full path to long sequence reads"
           "pairs [see below] (Multi-FASTA/fastq format, required)\n"
           "  -d  distance between k-mer pairs (ie. target distances to "
           "re-scaffold on. default -d 4000, optional)\n"
           "  \tMultiple distances are separated by comma. eg. -d "
           "500,1000,2000,3000\n"
           "  -j  threads  (default -j 3, optional) "
           "  -k  k-mer value (default -k 15, optional)\n"
           "  -t  step of sliding window when extracting k-mer pairs from long "
           "reads (default -t 2, optional)\n"
           "  \tMultiple steps are separated by comma. eg. -t 10,5\n"
           "  -o  offset position for extracting k-mer pairs (default -o 0, "
           "optional)\n"
           "  -e  error (%) allowed on -d distance   e.g. -e 0.1  == distance "
           "+/- 10% (default -e 0.1, optional)\n"
           "  -l  minimum number of links (k-mer pairs) to compute scaffold "
           "(default -l 5, optional)\n"
           "  -a  maximum link ratio between two best contig pairs (default -a "
           "0.3, optional)\n"
           "  \t *higher values lead to least accurate scaffolding*\n"
           "  -z  minimum contig length to consider for scaffolding (default "
           "-z 500, optional)\n"
           "  -b  base name for your output files (optional)\n"
           "  -r  Bloom filter input file for sequences supplied in -s "
           "(optional, if none provided will output to .bloom)\n"
           "  \t NOTE: BLOOM FILTER MUST BE DERIVED FROM THE SAME FILE "
           "SUPPLIED IN -f WITH SAME -k VALUE\n"
           "  \t IF YOU DO NOT SUPPLY A BLOOM FILTER, ONE WILL BE CREATED "
           "(.bloom)\n"
           "  -p  Bloom filter false positive rate (default -p 0.001, "
           "optional; increase to prevent memory allocation errors)\n"
           "  -x  Turn off Bloom filter functionality (-x 1 = yes, default = "
           "no, optional)\n"
           "  -v  Runs in verbose mode (-v 1 = yes, default = no, optional)\n";
  }
  void print_opts() {
    std::cout << "\nLINKS running options:"
              << "\n";
    std::cout << "  -f " << assembly_file << "\n"
              << "  -s " << fof_file << "\n"
              << "  -d " << distances_text << "\n"
              << "  -k " << k << "\n"
              << "  -j " << thread << "\n"
              << "  -t " << step_sizes_text << "\n"
              << "  -o " << offset << "\n"
              << "  -e " << insert_stdev << "\n"
              << "  -l " << min_links << "\n"
              << "  -a " << max_link_ratio << "\n"
              << "  -z " << min_size << "\n"
              << "  -b " << base_name << "\n"
              << "  -r " << bf_file << "\n"
              << "  -p " << fpr << "\n"
              << "  -x " << bf_off << "\n"
              << "  -v " << verbose << "\n";
  }

  InputParser() = default;

  InputParser(int &argc, char **argv) {
    this->argc = argc;
    static const struct option longopts[] = {{"help", no_argument, &help, 1},
                                             {nullptr, 0, nullptr, 0}};
    while ((c = getopt_long(argc, argv, "f:s:d:k:t:j:o:e:l:a:z:b:r:p:x:v:",
                            longopts, &optindex)) != -1) {
      switch (c) {
      case 0:
        break;
      case 'f':
        assembly_file.assign(optarg);
        break;
      case 's':
        // full path for fof
        fof_file.assign(optarg);
        long_reads = parse_fof_input(fof_file);
        break;
      case 'd':
        distances_text.assign(optarg);
        distances = split_distance_input(distances_text);
        break;
      case 'k':
        k = strtoul(optarg, &end, BASE_TEN);
        break;
      case 't':
        step_sizes_text.assign(optarg);
        step_sizes = split_step_sizes_input(step_sizes_text);
        break;
      case 'j':
        thread = strtoul(optarg, &end, BASE_TEN);
        break;
      case 'o':
        offset = strtoul(optarg, &end, BASE_TEN);
        break;
      case 'e':
        insert_stdev = strtof(optarg, &end);
        break;
      case 'l':
        min_links = strtoul(optarg, &end, BASE_TEN);
        break;
      case 'a':
        max_link_ratio = strtof(optarg, &end);
        break;
      case 'z':
        min_size = strtoul(optarg, &end, BASE_TEN);
        break;
      case 'b':
        base_name.assign(optarg);
        break;
      case 'r':
        bf_file.assign(optarg);
        break;
      case 'p':
        fpr = strtof(optarg, &end);
        break;
      case 'x':
        if (strtoul(optarg, &end, BASE_TEN) == 1) {
          bf_off = true;
        }
        break;
      case 'v':
        if (strtoul(optarg, &end, BASE_TEN) == 1) {
          verbose = true;
        }
        break;
      default:
        exit(EXIT_FAILURE);
      }
    }
    if (distances.size() > step_sizes.size()) {
      uint last_step = step_sizes.back();
      for (int i = step_sizes.size(); i < distances.size(); i++) {
        step_sizes.push_back(last_step);
      }
    }
    std::ifstream fof_file_file(fof_file);
    if (fof_file_file.peek() == std::ifstream::traits_type::eof()) {
      std::cerr << "\n File of files for reads cannot be empty (-s)\n";
      arguments_satisfied = false;
    } else {
      std::cout << "here with: " << fof_file << std::endl;
    }
    if (distances.size() < step_sizes.size()) {
      std::cerr << "\n Number of provided distances can't be lower than number "
                   "of step sizes provided.\n";
      arguments_satisfied = false;
    }
    if (base_name == "") {
      // copied from v1.8.7 -- add pid to basename
      base_name =
          assembly_file + "scaff_s-" + fof_file + "_d" + distances_text + "_k" +
          std::to_string(k) + "_e" + std::to_string(insert_stdev) + "_l" +
          std::to_string(min_links) + "_a" + std::to_string(max_link_ratio) +
          "_z" + std::to_string(min_size) + "_t" + step_sizes_text + "_o" +
          std::to_string(offset) + "_r-" + (bf_off ? "-" : bf_file) + "_p" +
          std::to_string(fpr) + "_x" + std::to_string(bf_off);
    }
    if (assembly_file == "" || fof_file == "") {
      print_usage();
      arguments_satisfied = false;
    }
  }
};
#endif
