//
// Created by adrian on 2/06/18.
//

#ifndef CINCT_RML_STRATEGY_HPP
#define CINCT_RML_STRATEGY_HPP

#include <util.hpp>
#include <string>
#include <fstream>
#include <vector>

namespace cinct {

    class rml_strategy {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;

    private:
        std::vector<std::vector<value_type >> m_mapping_edges;

        size_type _compute_nodes(const std::string input){
            std::ifstream i_stream(input);
            cnt::util::map_node_type map_nodes;
            cnt::util::map_node_type::iterator it;
            size_type number = 0;
            while(i_stream) {
                value_type x;
                i_stream.read((char *) &x, sizeof(int32_t));
                if (i_stream.eof()) break;
                if((it = map_nodes.find(x)) == map_nodes.end()){
                    map_nodes.insert(std::pair<value_type, char>(x, 'a'));
                    ++number;
                }
            }
            return number;
        }

        void _compute_frequencies(cnt::util::map_n_gram_type &n_gram_freq, const std::string input){
            cnt::util::map_n_gram_type::iterator it;
            std::ifstream i_stream(input);
            uint64_t traj_i = 0, max = 0;
            std::vector<uint32_t > v = {0,1};
            uint32_t first = 0;
            n_gram_freq.insert(std::pair<cnt::util::n_gram, uint32_t>(cnt::util::n_gram(v,2), 1));
            //First
            // i_stream.read((char*) &x[1], sizeof(value_type));

            int64_t pos = i_stream.tellg();
            i_stream.read((char*) &first, sizeof(value_type));
            std::cout << "first: " << first << std::endl;
            v[0] = first;
            v[1] = 0;
            n_gram_freq.insert(std::pair<cnt::util::n_gram, uint32_t>(cnt::util::n_gram(v,2), 1));

            uint32_t x[2] = {0, 0};
            i_stream.seekg(pos);
            while(i_stream){
                v.clear();
                i_stream.read((char*) &x, sizeof(value_type)*2);
                v.push_back(x[1]);
                v.push_back(x[0]);
                if(i_stream.eof()) break;
                pos += sizeof(value_type);
                i_stream.seekg(pos);
                cnt::util::n_gram n_g(v,2);
                if((it = n_gram_freq.find(n_g)) != n_gram_freq.end()){
                    ++it->second;
                }else{
                    n_gram_freq.insert(std::pair<cnt::util::n_gram, uint32_t>(n_g, 1));
                }
            }
            i_stream.close();
        }

        void _sort_frequencies(cnt::util::map_n_gram_type &n_gram_freq, cnt::util::keys_n_gram_type &keys){

            cnt::util::map_n_gram_type::iterator it, next_it;
            it = n_gram_freq.begin();
            uint64_t i = 0;
            while(it != n_gram_freq.end()){
                keys.push_back(*it);
                next_it = std::next(it);
                n_gram_freq.erase(it);
                it = next_it;
            }
            n_gram_freq.clear();

            std::sort(keys.begin(), keys.end(), cnt::util::compare_by_frequency);
        }

        void _build_mapping(cnt::util::keys_n_gram_type &keys, const size_type n_nodes){
            m_mapping_edges.resize(n_nodes+1, std::vector<value_type >()); //+1 because we start at 1
            for(size_type i = 0; i < keys.size(); ++i){
                uint32_t index = keys[i].first[0];
                uint32_t value = keys[i].first[1];
                m_mapping_edges[index].push_back(value);
            }
            keys.clear();
        }



    public:

        const std::vector<std::vector<value_type >> &mapping_edges = m_mapping_edges;

        rml_strategy(){};
        rml_strategy(std::string input){
            size_type n_nodes = _compute_nodes(input);
            std::cout << "Nodes: " << n_nodes << std::endl;
            cnt::util::map_n_gram_type n_gram_freq;
            _compute_frequencies(n_gram_freq, input);
            std::cout << "Frequencies: " << n_nodes << std::endl;
            cnt::util::keys_n_gram_type keys;
            _sort_frequencies(n_gram_freq, keys);
            std::cout << "Sort frequencies: " << n_nodes << std::endl;
            _build_mapping(keys, n_nodes);
            std::cout << "Mapping: " << n_nodes << std::endl;
        }

        void apply_labels(std::string input, std::string output){
            std::ifstream i_stream(input);
            std::ofstream o_stream(output);
            value_type value, prev_value = 0,label;
            while(i_stream){
                i_stream.read((char*) &value, sizeof(value_type));
                if(i_stream.eof()) break;
                label = 1;
                while(label-1 < m_mapping_edges[prev_value].size()
                      && m_mapping_edges[prev_value][label-1] != value){
                    ++label;
                }
                //if(label-1 == m_mapping_edges[prev_value].size()) exit(1);
                o_stream.write((char*) &label, sizeof(value_type));
                prev_value = value;
            }
            i_stream.close();
            o_stream.close();

        }


        rml_strategy& operator=(rml_strategy &&p){
            if (this != &p){
                m_mapping_edges = std::move(p.m_mapping_edges);
            }
            return *this;
        }

        //! Assignment operator
        rml_strategy& operator=(const rml_strategy& p)
        {
            if (this != &p) {
                m_mapping_edges = p.m_mapping_edges;
            }
            return *this;
        }

        //! Copy constructor
        rml_strategy(const rml_strategy& p)
        {
            m_mapping_edges = p.m_mapping_edges;
        }

        //! Move constructor
        rml_strategy(rml_strategy&& p)
        {
            *this = std::move(p);
        }

        //! Swap method
        void swap(rml_strategy& p)
        {
            if(this != &p){
                std::swap(m_mapping_edges, p.m_mapping_edges);
            }
        }



    };
}

#endif //CINCT_RML_STRATEGY_HPP
