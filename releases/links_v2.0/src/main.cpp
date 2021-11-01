#include "LINKS.hpp"
#include "InputParser.hpp"
#include <chrono>
#include <stdlib.h> // for sleep

int main(int argc, char** argv){
    InputParser* linksArgParser = new InputParser(argc, argv);
    if(!linksArgParser->arguments_satisfied){return -1;}
    linksArgParser->print_opts();
    std::cout << "fasta file: " << linksArgParser->assembly_file << std::endl;
    LINKS* links = new LINKS(linksArgParser);
    std::cout << "back from LINKS create bf\n";

    auto start = std::chrono::high_resolution_clock::now();
    links->init_bloom_filter();
    std::cout << "back from init bf\n";
    auto finish = std::chrono::high_resolution_clock::now();
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << milliseconds.count() << " ms\n";

    std::cout << "Read fasta stage\n";
    start = std::chrono::high_resolution_clock::now();
    links->start_read_fasta();
    std::cout << "Read fasta stage finish\n";
    finish = std::chrono::high_resolution_clock::now();
    milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << "Read fasta stage took: " << milliseconds.count() << " ms\n";

    std::cout << "Read contig stage\n";
    start = std::chrono::high_resolution_clock::now();
    links->start_read_contig();
    finish = std::chrono::high_resolution_clock::now();
    milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << "Read contig stage took: " << milliseconds.count() << " ms\n";

    std::cout << "Pair contigs stage\n";
    start = std::chrono::high_resolution_clock::now();
    links->pair_contigs();
    finish = std::chrono::high_resolution_clock::now();
    milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << "Pair contigs stage took: " << milliseconds.count() << " ms\n";

    delete(linksArgParser);
    delete(links);
    return 0;
}