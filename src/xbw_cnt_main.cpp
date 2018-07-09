//
// Created by adrian on 3/07/18.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stack>
#include <algorithm>
#include <xbw_cnt.hpp>
#include <cinct_rml.hpp>


void write_to_file(const std::vector<uint32_t> &v, const std::string file){
    std::ofstream out(file);
    out.write((char*) v.data(), v.size() * sizeof(uint32_t));
    out.close();
}



int main(int argc, char **argv) {

    std::vector<uint32_t> w1 = {7, 3, 5, 2, 1, 7, 3, 2, 1, 7, 3, 6, 8, 1, 7, 3, 5, 9, 1, 7, 4, 8, 1,
                                7, 4, 5, 9, 1, 7, 10, 5, 8, 1, 5, 10, 5, 8, 1, 3, 6, 1};
    write_to_file(w1, "proba.txt");

    std::vector<uint32_t > w1_rev;
    std::vector<uint32_t > aux;
    for(uint64_t i = 0; i < w1.size(); ++i){
        if(w1[i] == 1){
            std::reverse(aux.begin(), aux.end());
            aux.push_back(1);
            w1_rev.insert(w1_rev.end(), aux.begin(), aux.end());
            aux.clear();
        }else{
            aux.push_back(w1[i]);
        }

    }
    write_to_file(w1_rev, "proba.txt.rev");

    xbw_trie::xbw_cnt<> xbwCnt("proba.txt");
    cinct::cinct_rml cinctRml("proba.txt");

    std::vector<uint32_t > pattern = {7,3};
    auto solution = xbwCnt.sub_path_search(pattern);
    for(uint64_t i = 0; i < xbwCnt.alpha.size(); ++i){
        std::cout << "alpha: " << xbwCnt.alpha[i] << std::endl;
    }
    for(uint64_t i = 0; i < solution.size(); ++i){
        std::cout << "xbw-values: " << solution[i] << std::endl;
        size_t ini, end;
        std::stack<std::pair<uint64_t, uint64_t >> stack_c;
        uint64_t w_i = pattern[pattern.size()-1];
        uint64_t v;
        xbwCnt.get_children(solution[i], ini, end);
        for(uint64_t j = end; j >= ini; --j){
            stack_c.push({j, w_i});
        }
        //std::cout << xbwCnt.alpha[solution[i]] << std::endl;
        //std::cout << xbwCnt.graph[3][xbwCnt.alpha[solution[i]]-1].value << std::endl;
        while(!stack_c.empty()){
            auto sol = stack_c.top();
            stack_c.pop();

            v = xbwCnt.alpha[sol.first];
            w_i = xbwCnt.graph[sol.second][v-1].value;
            //std::cout << xbwCnt.alpha[sol.first] << std::endl;
            //std::cout << sol.second << std::endl;
           // std::cout << "pos: " << sol.first << std::endl;
            std::cout << "sol: " << w_i << std::endl;
            //std::cout << "size: " << xbwCnt.graph[w_i].size() << std::endl;
            /*std::cout << "sol: " <<  w_i << std::endl;
            w_i = xbwCnt.graph[sol.second][xbwCnt.alpha[sol.first]-1].value;
            std::cout << "sol: " <<  w_i << std::endl;*/
            xbwCnt.get_children(sol.first, ini, end);
            for(uint64_t j = ini; j <= end; ++j){
                stack_c.push({j, w_i});
            }
            if(ini > end){
                std::cout << "solution" << std::endl;
            }
        }


    }

    uint64_t start, end;
    std::reverse(pattern.begin(), pattern.end());
    if(cinctRml.labeled_search(pattern.begin(), pattern.end(), start, end)){
        std::cout << "start: " << start << " end: " << end << std::endl;
        for(uint64_t pos = start; pos < end; ++pos){
            std::cout << "cinct-values: " << pos << std::endl;
            std::vector<uint32_t > r;
            cinctRml.extract_traj(pos, r);
            for(uint64_t i = 0; i < r.size(); ++i){
                std::cout << "i: " << r[i] << std::endl;
            }
            std::cout << std::endl;
        }

    }






}