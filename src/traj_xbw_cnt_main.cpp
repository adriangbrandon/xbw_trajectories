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
#include <xbw_cnt.hpp>

int main(int argc, char **argv) {

       //xbw_trie::xbw<uint32_t, std::vector<uint32_t>, sdsl::wt_huff_int<sdsl::rrr_vector<> >> xbwa(trie, 1);
    //xbw_trie::xbw<uint32_t, std::vector<uint32_t>, sdsl::wt_ap<sdsl::wt_huff<sdsl::rrr_vector<>> >> xbwa(trie, 1);
    //xbw_trie::xbw<uint32_t, std::vector<uint32_t>, sdsl::wt_huff_int<sdsl::rrr_vector<> >> xbwa(trie, 1);
    //xbw_trie::xbw<uint32_t, std::vector<uint32_t>,
    //        sdsl::wt_ap<sdsl::wt_huff<sdsl::rrr_vector<>>, sdsl::wt_huff_int<sdsl::rrr_vector<>>>> xbwa(trie, 1);
    //xbw_trie::xbw<uint32_t, std::vector<uint32_t>, sdsl::wt_ap<sdsl::wt_huff< sdsl::hyb_vector<>> >> xbwa(trie, 1);
    //xbw_trie::xbw_cnt<uint32_t, std::vector<uint32_t>, sdsl::wt_ap<sdsl::wt_huff< sdsl::rrr_vector<>> >> xbwa(argv[1]);
        xbw_trie::xbw_cnt<uint32_t, std::vector<uint32_t>, sdsl::wt_huff_int< sdsl::rrr_vector<> >> xbwa(argv[1]);
    std::cout << "Size in bytes: " << sdsl::size_in_bytes(xbwa) << std::endl;

    std::ofstream out("size_xbw_cnt.html");
    sdsl::write_structure<sdsl::format_type::HTML_FORMAT>(xbwa, out);

}