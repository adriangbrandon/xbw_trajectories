//
// Created by adrian on 3/06/18.
//

#ifndef XBW_TRIE_GRAPH_HPP
#define XBW_TRIE_GRAPH_HPP

#include <vector>
#include <sdsl/int_vector.hpp>

namespace cnt {

    class node_graph {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;
        typedef int32_t z_type;

    private:
        value_type m_value;
        z_type m_z;

    public:
        const value_type &value = m_value;
        const z_type &z = m_z;

        node_graph(){};
        node_graph(value_type value, z_type z){
            m_value = value;
            m_z = z;
        }
        void set_z(const z_type z){
            m_z = z;
        }

        node_graph& operator=(node_graph &&p){
            if (this != &p){
                m_value = std::move(p.m_value);
                m_z = std::move(p.m_z);
            }
            return *this;
        }

        //! Assignment operator
        node_graph& operator=(const node_graph& p)
        {
            if (this != &p) {
                m_value = p.m_value;
                m_z = p.m_z;
            }
            return *this;
        }

        //! Copy constructor
        node_graph(const node_graph& p)
        {
            m_value = p.m_value;
            m_z = p.m_z;
        }

        //! Move constructor
        node_graph(node_graph&& p)
        {
            *this = std::move(p);
        }

        //! Swap method
        void swap(node_graph& p)
        {
            if(this != &p){
                std::swap(m_value, p.m_value);
                std::swap(m_z, p.m_z);
            }
        }

        //! Serializes the cinct_naive index to a stream.
        /*!
         * \param out The outstream to which the data structure is written.
         * \return The number of bytes written to out.
         */
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            uint64_t written_bytes = 0;
            written_bytes += sdsl::write_member(m_value, out, child, "value");
            written_bytes += sdsl::write_member(m_z, out, child, "z");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in)
        {
            sdsl::read_member(m_value, in);
            sdsl::read_member(m_z, in);
        }

       /* uint32_t operator[](uint64_t i) const {
            return m_values[i];
        }*/
    };

    class graph {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;
        typedef node_graph::z_type z_type;
        typedef std::vector<std::vector<node_graph >> values_type;

    private:
        values_type m_values;


    public:

        const values_type &values = m_values;

        graph(){};
        graph(const std::vector<std::vector<value_type >> &mapping_edges){
            m_values.resize(mapping_edges.size());
            for(size_type i = 0; i < mapping_edges.size(); ++i){
                for(size_type j = 0; j < mapping_edges[i].size(); ++j){
                    m_values[i].push_back(node_graph(mapping_edges[i][j], 0));
                }
                m_values[i].shrink_to_fit();
            }
            m_values.shrink_to_fit();
        }

        size_type size(){
            return m_values.size();
        }

        graph& operator=(graph &&p){
            if (this != &p){
                m_values = std::move(p.m_values);
            }
            return *this;
        }

        //! Assignment operator
        graph& operator=(const graph& p)
        {
            if (this != &p) {
                m_values = p.m_values;
            }
            return *this;
        }

        //! Copy constructor
        graph(const graph& p)
        {
            m_values = p.m_values;
        }

        //! Move constructor
        graph(graph&& p)
        {
            *this = std::move(p);
        }

        //! Swap method
        void swap(graph& p)
        {
            if(this != &p){
                std::swap(m_values, p.m_values);
            }
        }

        const values_type::value_type& operator[](size_type i) const {
            return m_values[i];
        }

        bool set_z(const value_type w, const value_type w_i, const z_type z){

            for(size_type i = 0; i < m_values[w_i].size(); ++i){
                if(m_values[w_i][i].value == w){
                    m_values[w_i][i].set_z(z);
                    return true;
                }
            }
            return false;
        }

        void set_z_at(const size_type  i, const value_type w_i, const z_type z){
            m_values[w_i][i].set_z(z);
        }

        size_type pos(const value_type w, const value_type w_i) const{
            for(size_type i = 0; i < m_values[w_i].size(); ++i){
                if(m_values[w_i][i].value == w){
                    return i+1;
                }
            }
            return 0;
        }

        value_type search_label(value_type prev_value, value_type value){
            value_type label = 1;
            while(label-1 < m_values[prev_value].size()
                  && m_values[prev_value][label-1].value != value){
                ++label;
            }
            if(label-1 == m_values[prev_value].size()){
                std::cerr << "No label found" << std::endl;
                exit(1);
            };
            return label;
        }

        //! Serializes the cinct_naive index to a stream.
        /*!
         * \param out The outstream to which the data structure is written.
         * \return The number of bytes written to out.
         */
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            uint64_t written_bytes = 0;
            sdsl::write_member(m_values.size(), out);
            for(uint64_t i = 0; i < m_values.size(); ++i){
                sdsl::write_member(m_values[i].size(), out);
                written_bytes += sdsl::serialize_vector(m_values[i], out, child, "values");
            }
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from istream
        /*! Serialize the graph from the given istream.
         *  \param istream Stream to load the cinct_rml from.
         */
        void load(std::istream& in)
        {
            size_type size;
            sdsl::read_member(size, in);
            m_values.resize(size);
            for(uint64_t i = 0; i < m_values.size(); ++i){
                sdsl::read_member(size, in);
                m_values[i].resize(size);
                sdsl::load_vector(m_values[i], in);
            }
        }

    };
}

#endif //CINCT_GRAPH_HPP
