//
// Created by adrian on 25/06/18.
//

#include <iostream>
#include <vector>
#include <trie.hpp>
#include <xbw.hpp>
#include <xbw_cnt.hpp>

template<class t_array>
void print_array(const t_array &array, const std::string name){
    std::cout << name << ": ";
    for(size_t i = 0; i < array.size(); ++i){
        std::cout << array[i]  << ", ";
    }
    std::cout << std::endl;
}

int main(int argc, char **argv) {

    size_t traj_size = (size_t) atoi(argv[1]);
    std::mt19937 rng;
    std::vector<std::vector<uint32_t >> words(traj_size, std::vector<uint32_t>());
    {
        std::uniform_int_distribution<uint64_t> distribution(2,4);
        auto dice = bind(distribution, rng);
        for(size_t i = 0; i < words.size(); ++i){
            for (size_t j=0; j<10; ++j) {
                words[i].push_back(dice());
            }
            words[i].push_back(1);
        }
    }

    std::ofstream out("proba.txt");
    for(uint64_t i = 0; i < words.size(); ++i){
        out.write((char*) words[i].data(), words[i].size() * sizeof(uint32_t));
    }
    out.close();

    std::vector<xbw_trie::xbw_node<>> nodes;
    xbw_trie::trie<> trie;
    for(size_t i = 0; i < words.size(); ++i){
        trie.insert(words[i]);
    }
    std::cout << "Trie created" << std::endl;
    xbw_trie::xbw<> xbwa(trie, 1, nodes);
    sdsl::util::clear(trie);
    xbw_trie::xbw_cnt<> xbwCnt("proba.txt");

    {
        std::uniform_int_distribution<uint64_t> distribution(2, 4);
        auto dice = bind(distribution, rng);
        for(size_t length = 1; length < 10; ++length){
            for(size_t i = 0; i < 1000; ++i){
                std::vector<uint32_t> pattern;

                for(size_t j = 0; j < length; ++j){
                    pattern.push_back(dice());
                }
                std::vector<uint64_t> solutions;
                for(uint64_t j = 0; j < nodes.size(); ++j){
                    if(nodes[j].contains(pattern)) solutions.push_back(j);
                }

                print_array(pattern, "Pattern");
                auto result = xbwa.sub_path_search(pattern);
                std::cout << "Total results (" << i << "): " << result.size() << std::endl;
                std::cout << std::endl;
                if(result.size() != solutions.size()){
                    std::cout << "Pattern: " << std::endl;
                    std::cout << std::endl;
                    std::cout << "Error size" << std::endl;
                    std::cout << "Expected: " << solutions.size() << std::endl;
                    std::cout << "Obtained: " << result.size() << std::endl;
                    print_array(result, "result");
                    exit(0);
                }
                for(uint64_t j = 0; j < result.size(); ++j){
                    if(result[j] != solutions[j]){
                        std::cout << "Error solution" << std::endl;
                        std::cout << "Expected: " << solutions.size() << std::endl;
                        std::cout << "Obtained: " << result.size() << std::endl;
                        print_array(result, "result");
                        exit(0);
                    }
                }

                auto result_cnt = xbwCnt.sub_path_search(pattern);
                if(result.size() != result_cnt.size()){
                    std::cout << "Pattern CNT: " << std::endl;
                    std::cout << std::endl;
                    std::cout << "Error size" << std::endl;
                    std::cout << "Expected: " << result_cnt.size() << std::endl;
                    std::cout << "Obtained: " << result.size() << std::endl;
                    print_array(result_cnt, "result_cnt");
                    print_array(result, "result");
                    exit(0);
                }
                for(uint64_t j = 0; j < result_cnt.size(); ++j){
                    if(result_cnt[j] != result[j]){
                        std::cout << "Error solution" << std::endl;
                        std::cout << "Expected: " << result_cnt.size() << std::endl;
                        std::cout << "Obtained: " << result.size() << std::endl;
                        print_array(result_cnt, "result_cnt");
                        print_array(result, "result");
                        exit(0);
                    }
                }
            }
        }
    }


    /*std::vector<uint32_t> w1 = {2, 3, 5, 7, 1};
    std::vector<uint32_t> w2 = {2, 3, 7, 1};
    std::vector<uint32_t> w3 = {2, 3, 6, 8, 1};
    std::vector<uint32_t> w4 = {2, 3, 5, 9, 1};
    std::vector<uint32_t> w5 = {2, 4, 8, 1};
    std::vector<uint32_t> w6 = {2, 4, 5, 9, 1};
    std::vector<uint32_t> w7 = {2, 10, 5, 8, 1};
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

    xbw_trie::xbw<> xbwa(trie, 1);


    for(uint64_t i = 0; i < 25; ++i){
        size_t ini, end;
        xbwa.get_children(i, ini, end);
        std::cout << "n: " << i << " Children ini: " << ini << " end: " << end << std::endl;
        std::cout << "n: " << i << " parent: " << xbwa.get_parent(i) << std::endl;
    }

    //std::vector<uint32_t> pattern = {2, 3};
    std::vector<uint32_t> pattern = {9,1};
    auto solution = xbwa.sub_path_search(pattern);
    std::cout << "Solution" << std::endl;
    for(uint64_t i = 0; i < solution.size(); ++i){
        std::cout << solution[i] << std::endl;
    }*/

}