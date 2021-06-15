#include "LINKS.hpp"
#include "InputParser.hpp"

int main(int argc, char** argv){
    InputParser* linksArgParser = new InputParser(argc, argv);
    std::cout << "fasta file: " << linksArgParser->assemblyFile << std::endl;
    LINKS* links = new LINKS(linksArgParser);
    std::cout << "test 2:" << std::endl;
    delete(linksArgParser);
    std::cout << "test 4: " << std::endl;
    delete(links);
    std::cout << "test 6: " << std::endl;
    return 0;
}