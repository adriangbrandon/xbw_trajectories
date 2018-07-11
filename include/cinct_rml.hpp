//
// Created by adrian on 3/06/18.
//

#ifndef CINCT_CINCT_RML_HPP
#define CINCT_CINCT_RML_HPP

#include <string>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <sdsl/config.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/csa_wt.hpp>
#include <rml_strategy.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <graph.hpp>

using namespace sdsl;
using namespace cnt;

namespace cinct {

    class cinct_rml {

    public:

        typedef uint64_t size_type;
        typedef uint32_t value_type;

    private:

        int_vector<> m_C;
        size_type m_sigma;
        graph m_graph;
        wt_huff_int<rrr_vector<> > m_wt_label_bwt;

        size_type _bytes_text(std::string input){
            std::ifstream i_stream(input, std::ifstream::ate | std::ifstream::binary);
            return static_cast<size_type >(i_stream.tellg());
        }

        void _construct_bwt(std::string input, std::string output, wt_int<> &aux_wt){
            uint8_t num_byte = 4; //32 bits integer (4 bytes)
            cache_config cc(false, "./", "construct_bwt");
            csa_wt_int<> fm;
            construct(fm, input, cc, num_byte);
            int_vector<> bwt;
            load_from_file(bwt, cache_file_name(conf::KEY_BWT_INT, cc));
            std::ofstream out(output);
            out.write((char*)bwt.data(), num_byte*bwt.size());
            out.close();
            m_C = fm.C;
            m_sigma = fm.sigma;
            aux_wt = fm.wavelet_tree;
            sdsl::util::delete_all_files(cc.file_map);
        }

        void _label_bwt(std::string input, std::string output){
            std::ifstream i_stream(input);
            std::ofstream o_stream(output);
            value_type value = 0, prev_value = -1,label;
            size_type pos = 0;
            while(i_stream){
                if(m_C[prev_value+1] == pos){
                    ++prev_value;
                }
                i_stream.read((char*) &value, sizeof(value_type));
                if(i_stream.eof()) break;
                label = m_graph.search_label(prev_value, value);
                o_stream.write((char*) &label, sizeof(value_type));
                pos++;
            }
            i_stream.close();
            o_stream.close();
        }

        void _construct_wt_label_bwt(std::string input){
            construct(m_wt_label_bwt, input, 4);
        }

        void _construct_graph(std::string input){
            rml_strategy rml_str(input);
            m_graph = graph(rml_str.mapping_edges);
            sdsl::util::clear(rml_str);
        }
        void _add_z_graph(const wt_int<> &aux_wt){
            for(value_type w_i = 0; w_i < m_graph.size(); ++w_i){
                for(size_type i = 0; i < m_graph[w_i].size(); ++i){
                    value_type w = m_graph[w_i][i].value;
                    int32_t z = m_wt_label_bwt.rank(m_C[w_i], i+1) - aux_wt.rank(m_C[w_i], w);
                    m_graph.set_z_at(i, w_i, z);
                    //std::cout << "z: " << z << std::endl;
                }
            }
        }






    public:

        const graph &graph_ptr = m_graph;
        cinct_rml(){};
        cinct_rml(std::string input ){
            std::string traj_rev = input + ".rev";
            std::string file_bwt = input + ".bwt";
            std::string file_label = input + ".lab";
            wt_int<> aux_wt;
            std::cout << "Building graph..."<< std::endl;
            _construct_graph(traj_rev);
            /* std::cout << "Graph" << std::endl;
             for(uint64_t i = 0; i < m_graph.size(); ++i){
                 std::cout << "w_i: " << i << " [";
                 for(uint64_t j = 0; j < m_graph[i].size(); ++j){
                     std::cout << "label: " << j << " value: " << m_graph[i][j].value << " z: " << m_graph[i][j].z << ",";
                 }
                 std::cout << "]" << std::endl;
             }
             std::cout << std::endl;*/
            std::cout << "Building BWT..."<< std::endl;
            _construct_bwt(traj_rev, file_bwt, aux_wt);
            std::cout << "Labeling BWT..."<< std::endl;
            _label_bwt(file_bwt, file_label);
            std::cout << "Bulding WT over BWT-labeling..."<< std::endl;
            _construct_wt_label_bwt(file_label);
            std::cout << "Computing Zs..."<< std::endl;
            _add_z_graph(aux_wt);
            std::cout << "CINCT has been built!!" << std::endl;
            sdsl::util::clear(aux_wt);

            /*std::cout << "C" << std::endl;
            for(uint64_t i = 0; i < m_C.size(); ++i){
                std::cout << m_C[i] << ", ";
            }
            std::cout << std::endl;



            std::cout << "BWT label" << std::endl;
            for(uint64_t i = 0; i < m_wt_label_bwt.size(); ++i){
                std::cout << m_wt_label_bwt[i] << ", ";
            }
            std::cout << std::endl;*/

            std::cout << "BWT label" << std::endl;
            uint64_t max_value = 0;
            for(uint64_t i = 0; i < m_wt_label_bwt.size(); ++i){
                uint64_t value = m_wt_label_bwt[i];
                if(max_value < value){
                    max_value = value;
                }
            }
            std::cout << "max_value: " << max_value << std::endl;



        }

        size_type pseudo_rank(const size_type pos, const value_type w, const value_type w_i){
            size_type label = m_graph.pos(w, w_i);
            if(label > 0 && m_C[w_i] <= pos && pos <= m_C[w_i+1]){
                return m_wt_label_bwt.rank(pos, label) - m_graph[w_i][label-1].z;
            }
            return 0;
        }

        size_type pseudo_rank_exist(const size_type pos, const size_type label, const value_type w_i){
            if(label > 0 && m_C[w_i] <= pos && pos <= m_C[w_i+1]){
                return m_wt_label_bwt.rank(pos, label) - m_graph[w_i][label-1].z;
            }
            return 0;
        }

        template <class iter_t>
        bool labeled_search(iter_t begin, iter_t end, size_type &l_res, size_type &r_res){
            size_type l, r;
            iter_t it = end;
            --it;
            value_type w = (*it), w_i;
            l = m_C[w]; r = m_C[w+1];
            while (begin < it and r+1-l > 0) {
                --it;
                w_i = w;
                w = (*it);
                size_type label = m_graph.pos(w, w_i);
                if(label == 0) return false;
                size_type c_begin = m_C[w];
                if (l == 0 and r+1 == m_wt_label_bwt.size()) {
                    l = c_begin;
                    r = m_C[w+1]-1;
                } else {
                    l = c_begin + pseudo_rank_exist(l, label, w_i);
                    r = c_begin + pseudo_rank_exist(r, label, w_i);
                }
            }
            if(l >= r) return false;
            l_res = l;
            r_res = r;
            return true;
        }

        template <class array_t>
        void extract(size_type start, size_type size, array_t &result){

            // binary search the character with C
            value_type upper_c = static_cast<value_type >(m_sigma), lower_c = 0; // lower_c inclusive, upper_c exclusive
            value_type w_i=0, w = 0;
            size_type z;
            do {
                w_i = (upper_c+lower_c)/2;
                if (start < m_C[w_i]) {
                    upper_c = w_i;
                } else if (start >= m_C[w_i+1]) {
                    lower_c = w_i+1;
                }
            } while (start< m_C[w_i] || start >= m_C[w_i+1]);
            size_type j = start;
            for(int64_t i = size-1; i >= 0; --i){
                size_type label = m_wt_label_bwt[j];
                w = m_graph[w_i][label-1].value;
                result[i] = w;
                j = m_C[w] + pseudo_rank_exist(j, label, w_i);
                w_i = w;
                if(i % 10000 == 0) std::cout << i << std::endl;
            }
        }


        template <class array_t>
        void extract_traj(const size_type start, array_t &result){

            // binary search the character with C
            value_type upper_c = static_cast<value_type >(m_sigma), lower_c = 0; // lower_c inclusive, upper_c exclusive
            value_type w_i=0, w = 2;
            size_type z;
            do {
                w_i = (upper_c+lower_c)/2;
                if (start < m_C[w_i]) {
                    upper_c = w_i;
                } else if (start >= m_C[w_i+1]) {
                    lower_c = w_i+1;
                }
            } while (start< m_C[w_i] || start >= m_C[w_i+1]);
            size_type j = start;
            while(true){
                size_type label = m_wt_label_bwt[j];
                w = m_graph[w_i][label-1].value;
                if(w < 2){
                    result.push_back(1);
                    break;
                };
                result.push_back(w);
                j = m_C[w] + pseudo_rank_exist(j, label, w_i);
                w_i = w;
            }
        }

        //! Serializes the cinct_rml index to a stream.
        /*!
         * \param out The outstream to which the data structure is written.
         * \return The number of bytes written to out.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            uint64_t written_bytes = 0;
            written_bytes += write_member(m_sigma, out, child, "sigma");
            written_bytes += m_C.serialize(out, child, "C");
            written_bytes += m_graph.serialize(out, child, "graph");
            written_bytes += m_wt_label_bwt.serialize(out, child, "hwt");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from istream
        /*! Serialize the cinct_rml from the given istream.
         *  \param istream Stream to load the cinct_rml from.
         */
        void load(std::istream& in)
        {
            read_member(m_sigma, in);
            m_C.load(in);
            m_graph.load(in);
            m_wt_label_bwt.load(in);
        }

        size_type size(){
            return m_wt_label_bwt.size();
        }


        void check_z_graph(){
            for(value_type w_i = 0; w_i < m_graph.size(); ++w_i){
                for(size_type i = 0; i < m_graph[w_i].size(); ++i){
                    if(m_graph[w_i][i].z < 0){
                        std::cout << "z: " <<  m_graph[w_i][i].z << std::endl;
                    }
                }
            }
        }

    };
}

#endif //CINCT_CINCT_RML_HPP
