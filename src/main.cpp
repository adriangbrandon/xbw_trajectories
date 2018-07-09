#include <iostream>
#include <vector>
#include <trie.hpp>
#include <xbw.hpp>


int main(int argc, char **argv) {

    std::vector<uint32_t> w1 = {7, 3, 5, 2, 1};
    std::vector<uint32_t> w2 = {7, 3, 2, 1};
    std::vector<uint32_t> w3 = {7, 3, 6, 8, 1};
    std::vector<uint32_t> w4 = {7, 3, 5, 9, 1};
    std::vector<uint32_t> w5 = {7, 4, 8, 1};
    std::vector<uint32_t> w6 = {7, 4, 5, 9, 1};
    std::vector<uint32_t> w7 = {7, 10, 5, 8, 1};
    std::vector<uint32_t> w8 = {5, 10, 5, 8, 1};

    xbw_trie::trie<> trie;
    trie.insert(w1);
    trie.insert(w2);
    trie.insert(w3);
    trie.insert(w4);
    trie.insert(w5);
    trie.insert(w6);
    trie.insert(w7);
    trie.insert(w8);

    //trie.print();
    std::cout << std::endl;
    std::vector<xbw_trie::xbw_node<>> nodes;
    xbw_trie::xbw<> xbwa(trie, 1, nodes);

    for(uint64_t i = 0; i < 25; ++i){
        size_t ini, end;
        xbwa.get_children(i, ini, end);
        std::cout << "n: " << i << " Children ini: " << ini << " end: " << end << std::endl;
        //std::cout << "n: " << i << " parent: " << xbwa.get_parent(i) << std::endl;
    }

    std::vector<uint32_t> pattern = {2, 3};
   //std::vector<uint32_t> pattern = {9,1};

    auto solution = xbwa.sub_path_search(pattern);
    std::cout << "Solution" << std::endl;
    for(uint64_t i = 0; i < solution.size(); ++i){
        std::cout << solution[i] << std::endl;
    }


    std::vector<uint64_t> solutions;
    for(uint64_t i = 0; i < nodes.size(); ++i){
        if(nodes[i].contains(pattern)) solutions.push_back(i);
    }

    std::cout << "Solutions" << std::endl;
    for(uint64_t i = 0; i < solutions.size(); ++i){
        std::cout << solutions[i] << std::endl;
    }

}