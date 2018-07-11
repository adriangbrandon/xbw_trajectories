//
// Created by adrian on 5/07/18.
//

//
// Created by adrian on 25/06/18.
//

#include <iostream>
#include <vector>
#include <trie.hpp>
#include <xbw.hpp>
#include <xbw_cnt.hpp>
#include <cinct_rml.hpp>

template<class t_array>
void print_array(const t_array &array, const std::string name){
    std::cout << name << ": ";
    for(size_t i = 0; i < array.size(); ++i){
        std::cout << array[i]  << ", ";
    }
    std::cout << std::endl;
}


void write_to_file(const std::vector<uint32_t> &v, const std::string file){
    std::ofstream out(file);
    out.write((char*) v.data(), v.size() * sizeof(uint32_t));
    out.close();
}


bool search_sol(const std::vector<uint32_t> &sol, const std::vector<std::vector<uint32_t>> & expected){
    return std::binary_search(expected.begin(), expected.end(), sol);
}

int main(int argc, char **argv) {

    size_t traj_size = (size_t) atoi(argv[1]);
    std::mt19937 rng;
    std::vector<uint32_t > words;
    {
        std::uniform_int_distribution<uint64_t> distribution(2,4);
        auto dice = bind(distribution, rng);
        for(size_t i = 0; i < traj_size; ++i){
            for (size_t j=0; j<40; ++j) {
                words.push_back(dice());
            }
            words.push_back(1);
        }
    }

    //TODO: remove this! only test of test
    //words = {7, 3, 5, 2, 1, 7, 3, 2, 1, 7, 3, 6, 8, 1, 7, 3, 5, 9, 1, 7, 4, 8, 1, 7, 4, 5, 9, 1, 7, 10, 5, 8, 1, 5, 10, 5, 8, 1, 3, 6, 1};

    write_to_file(words, "proba.txt");
    print_array(words, "words");

    std::vector<uint32_t > w1_rev;
    std::vector<uint32_t > aux;
    for(uint64_t i = 0; i < words.size(); ++i){
        if(words[i] == 1){
            std::reverse(aux.begin(), aux.end());
            aux.push_back(1);
            w1_rev.insert(w1_rev.end(), aux.begin(), aux.end());
            aux.clear();
        }else{
            aux.push_back(words[i]);
        }

    }

    write_to_file(w1_rev, "proba.txt.rev");
    xbw_trie::xbw_cnt<> xbwCnt("proba.txt");
    cinct::cinct_rml cinctRml("proba.txt");

    std::cout << "The structures were built" << std::endl;
    {
        std::uniform_int_distribution<uint64_t> distribution(2, 4);
        auto dice = bind(distribution, rng);
        for(size_t length = 1; length < 10; ++length){
            for(size_t i = 0; i < 1000; ++i){
                std::vector<uint32_t> pattern, pattern_rev;
                for(size_t j = 0; j < length; ++j){
                    pattern.push_back(dice());
                }
                //TODO: remove this! only test of test
                //pattern = {7,3};
                pattern_rev.resize(pattern.size());
                std::copy(pattern.begin(), pattern.end(), pattern_rev.begin());
                std::reverse(pattern_rev.begin(), pattern_rev.end());
                std::vector<std::vector<uint32_t >> expected;

                uint64_t start, end;
                if(cinctRml.labeled_search(pattern_rev.begin(), pattern_rev.end(), start, end)){
                    for(uint64_t pos = start; pos < end; ++pos){
                        std::vector<uint32_t > r;
                        r.resize(pattern.size());
                        std::copy(pattern.begin(), pattern.end(), r.begin());
                        cinctRml.extract_traj(pos, r);
                        expected.push_back(r);
                    }
                }
                std::cout << "expected.size: " << expected.size() << std::endl;
                //print_array(*expected.end(), "proando");
                std::sort(expected.begin(), expected.end());
                /*for(uint64_t j = 0; j < expected.size(); ++j){
                    std::cout << "expected[" << j << "]: ";
                    for(uint64_t jj = 0; jj < expected[j].size(); ++jj){
                        std::cout << expected[j][jj] << ", ";
                    }
                    std::cout << std::endl;
                }*/

                auto solution = xbwCnt.sub_path_search(pattern);
                auto sol = xbwCnt.get_trajectories(solution, pattern);
                for(uint64_t j = 0; j < sol.size(); ++j){
                    //print_array(sol[j], "sol["+std::to_string(j)+"]");
                    if(search_sol(sol[j], expected)){
                       // std::cout << "OK!" << std::endl;
                    }else{
                        std::cout << "Error" << std::endl;
                        exit(10);
                    }
                }

                //Double check
                std::sort(sol.begin(), sol.end());
                for(uint64_t j = 0; j < expected.size(); ++j){
                    //print_array(expected[j], "expected["+std::to_string(j)+"]");
                    if(search_sol(expected[j], sol)){
                        //std::cout << "OK!" << std::endl;
                    }else{
                        std::cout << "Error" << std::endl;
                        exit(10);
                    }
                }
            }
        }

        std::cout << "Size XBWT: " << sdsl::size_in_bytes(xbwCnt) << std::endl;
        std::cout << "Size CiNCT: " << sdsl::size_in_bytes(cinctRml) << std::endl;

        std::cout << "Everything is OK!" << std::endl;
    }


}