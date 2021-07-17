#include "LINKS.hpp"
#include "InputParser.hpp"
#include <chrono>
#include <stdlib.h> // for sleep

int main(int argc, char** argv){
    InputParser* linksArgParser = new InputParser(argc, argv);
    std::cout << "fasta file: " << linksArgParser->assemblyFile << std::endl;
    LINKS* links = new LINKS(linksArgParser);

    auto start = std::chrono::high_resolution_clock::now();
    links->init_bloom_filter();
    std::cout << "back from init bf\n";
    auto finish = std::chrono::high_resolution_clock::now();
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << milliseconds.count() << "ms\n";

    std::cout << "Read fasta stage\n";
    start = std::chrono::high_resolution_clock::now();
    links->start_read_fasta();
    std::cout << "Read fasta stage finish\n";
    finish = std::chrono::high_resolution_clock::now();
    milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << milliseconds.count() << "ms\n";

    std::cout << "Read contig stage\n";
    start = std::chrono::high_resolution_clock::now();
    links->start_read_contig();
    finish = std::chrono::high_resolution_clock::now();
    milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << milliseconds.count() << "ms\n";

    std::cout << "Pair contigs stage\n";
    start = std::chrono::high_resolution_clock::now();
    links->pair_contigs();
    finish = std::chrono::high_resolution_clock::now();
    milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << milliseconds.count() << "ms\n";

    delete(linksArgParser);
    std::cout << "test 4: " << std::endl;
    delete(links);
    std::cout << "test 6: " << std::endl;
    return 0;
}