//
// Created by adrian on 28/06/18.
//

#include <fstream>
#include <iostream>
#include <vector>
#include <trie.hpp>
#include <xbw.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>

int main(int argc, char **argv) {

    std::ifstream input(argv[1]);
    std::ifstream original(argv[2]);


    uint32_t x, orig;
    std::vector<uint32_t> w;

    xbw_trie::trie<> trie;
    while(input){
        input.read((char*) &x, sizeof(uint32_t));
        original.read((char*) &orig, sizeof(uint32_t));
        if(input.eof()) break;
        w.push_back(x);
        if(orig == 1){
            trie.insert(w);
            w.clear();
        }
    }
    input.close();
    original.close();

    //xbw_trie::xbw<uint32_t, std::vector<uint32_t>, sdsl::wt_huff_int<sdsl::rrr_vector<> >> xbwa(trie, 1);
    //xbw_trie::xbw<uint32_t, std::vector<uint32_t>, sdsl::wt_ap<sdsl::wt_huff<sdsl::rrr_vector<>> >> xbwa(trie, 1);
    //xbw_trie::xbw<uint32_t, std::vector<uint32_t>, sdsl::wt_huff_int<sdsl::rrr_vector<> >> xbwa(trie, 1);
    //xbw_trie::xbw<uint32_t, std::vector<uint32_t>,
    //        sdsl::wt_ap<sdsl::wt_huff<sdsl::rrr_vector<>>, sdsl::wt_huff_int<sdsl::rrr_vector<>>>> xbwa(trie, 1);
    xbw_trie::xbw<uint32_t, std::vector<uint32_t>, sdsl::wt_ap<sdsl::wt_huff< sdsl::hyb_vector<>> >> xbwa(trie, 1);
    sdsl::util::clear(trie);
    std::cout << "Size in bytes: " << sdsl::size_in_bytes(xbwa) << std::endl;

    std::ofstream out("size_xbw.html");
    sdsl::write_structure<sdsl::format_type::HTML_FORMAT>(xbwa, out);

}