#include "LINKS.hpp"
#include "InputParser.hpp"
#include <chrono>

int main(int argc, char** argv){
    InputParser* linksArgParser = new InputParser(argc, argv);
    std::cout << "fasta file: " << linksArgParser->assemblyFile << std::endl;
    LINKS* links = new LINKS(linksArgParser);

    auto start = std::chrono::high_resolution_clock::now();
    links->init_bloom_filter();
    auto finish = std::chrono::high_resolution_clock::now();
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << milliseconds.count() << "ms\n";

    start = std::chrono::high_resolution_clock::now();
    links->start_read_fasta();
    finish = std::chrono::high_resolution_clock::now();
    milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
    std::cout << milliseconds.count() << "ms\n";

    delete(linksArgParser);
    std::cout << "test 4: " << std::endl;
    delete(links);
    std::cout << "test 6: " << std::endl;
    return 0;
}