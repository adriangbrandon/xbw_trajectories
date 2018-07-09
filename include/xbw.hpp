//
// Created by adrian on 23/06/18.
//

#ifndef XBW_TRIE_XBW_HPP
#define XBW_TRIE_XBW_HPP

#define VERBOSE 1


#include <trie.hpp>
#include <stack>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/sd_vector.hpp>


namespace xbw_trie {


    template<class t_array>
    void print_array(const t_array &array, const std::string name){
        for(size_t i = 0; i < array.size(); ++i){
            std::cout << name << "[" << i << "]: " << array[i] << std::endl;
        }
    }

    template <class t_value = uint32_t,
              class t_word = std::vector<t_value> >
    class xbw_node {

    public:
        typedef t_value value_type;
        typedef t_word  word_type;

    private:
        bool m_is_right_most = false;
        value_type m_alpha;
        word_type m_pi;

    public:
        const bool& is_right_most = m_is_right_most;
        const value_type& alpha = m_alpha;
        const word_type& pi = m_pi;

        xbw_node(){};
        xbw_node(const bool rm, const value_type a, const word_type &p){
            m_is_right_most = rm;
            m_alpha = a;
            m_pi = p;
        }


        bool operator<(const xbw_node &p) const {
            size_t t_i = m_pi.size();
            size_t n_i = p.pi.size();
            //Lenth substring
            size_t n = t_i;
            bool r_aux = true;
            if(t_i > n_i){
                n = n_i;
                r_aux = false;
            }
            //Compare substrings
            for(size_t i = 0; i < n; ++i){
                if(m_pi[t_i-i-1] < p.pi[n_i-i-1]) return true;
                if(m_pi[t_i-i-1] > p.pi[n_i-i-1]) return false;
            }
            //Compare last element when the substrings are identical
            if(t_i == n_i) return m_is_right_most < p.is_right_most;
            //If substrings are different
            return r_aux;
        }

        value_type pi_first(){
            return m_pi[m_pi.size()-1];
        }

        //4,2,6
        bool contains(const word_type &pattern){
            size_t length = pattern.size();
            if(length > m_pi.size()+1) return false;
            if(m_alpha != pattern[length-1]){
                return false;
            }
            size_t start = m_pi.size() - (length-1);
            for(size_t i = 0; i < length-1; ++i){
                if(m_pi[start+i] != pattern[i]){
                    return false;
                }
            }
            return true;
        }

        xbw_node& operator=(xbw_node &&p){
            if (this != &p){
                m_is_right_most = std::move(p.m_is_right_most);
                m_alpha = std::move(p.m_alpha);
                m_pi = std::move(p.m_pi);
            }
            return *this;
        }

        //! Assignment operator
        xbw_node& operator=(const xbw_node& p)
        {
            if (this != &p) {
                m_is_right_most = p.m_is_right_most;
                m_alpha = p.m_alpha;
                m_pi = p.m_pi;
            }
            return *this;
        }

        //! Copy constructor
        xbw_node(const xbw_node& p)
        {
            m_is_right_most = p.m_is_right_most;
            m_alpha = p.m_alpha;
            m_pi = p.m_pi;
        }

        //! Move constructor
        xbw_node(xbw_node&& p)
        {
            *this = std::move(p);
        }

        //! Swap method
        void swap(xbw_node& p)
        {
            if(this != &p){
                std::swap(m_is_right_most, p.m_is_right_most);
                std::swap(m_alpha, p.m_alpha);
                std::swap(m_pi, p.m_pi);
            }
        }



        void print(){
            std::cout << "right: " << m_is_right_most << " alpha: " << m_alpha << " pi: [";
            for(uint64_t j = 0; j < m_pi.size(); ++j){
                std::cout << m_pi[j] << ",";
            }
            std::cout << "]" << std::endl;
        }
    };

    template <class t_value = uint32_t,
              class t_word = std::vector<t_value>,
              class t_wt = sdsl::wt_blcd_int< sdsl::rrr_vector<>>,
              class t_bv = sdsl::bit_vector>
    class xbw {

    public:
        typedef t_value value_type;
        typedef t_word word_type;
        typedef size_t size_type;
        typedef std::vector<xbw_node<t_value, t_word>> nodes_type;

    private:

        t_wt m_alpha;
        t_bv m_last;
        sdsl::rank_support_v5<> m_rank1_last;
        sdsl::select_support_mcl<> m_select1_last;
        sdsl::sd_vector<> m_a;
        sdsl::sd_vector<>::select_1_type m_select1_a;
        sdsl::sd_vector<>::rank_1_type m_rank1_a;
        value_type m_leaf;

    public:
        const t_wt& alpha = m_alpha;
        const t_bv& last = m_last;
        const sdsl::rank_support_v5<>& rank1_last = m_rank1_last;
        const sdsl::select_support_mcl<>& select1_last = m_select1_last;
        const sdsl::sd_vector<>& a = m_a;
        const sdsl::sd_vector<>::select_1_type& select1_a = m_select1_a;
        const sdsl::sd_vector<>::rank_1_type& rank1_a = m_rank1_a;
        const value_type& leaf = m_leaf;

        struct stack_element{
            bool is_right_most_child;
            value_type alpha;
            const trie_node<value_type, word_type>* node;
            value_type level;
        };


        void _build_nodes(const trie<t_value, t_word>& trie, nodes_type &nodes) {
            std::stack<stack_element> trie_nodes_stack;
            trie_nodes_stack.push({0, 0, &trie.root, 0});
            word_type pi;
            size_type prev_level = 0;
            //Pre-order
            while(!trie_nodes_stack.empty()){
                stack_element se = trie_nodes_stack.top();
                trie_nodes_stack.pop();
                if(se.level){
                    if(prev_level >= se.level){
                        pi.resize(se.level-1);
                    }
                    xbw_node<value_type, word_type> node(se.is_right_most_child, se.alpha, pi);
                    nodes.push_back(node);
                    pi.push_back(se.alpha);
                }
                for(auto it = se.node->children.begin(); it != se.node->children.end(); ++it){
                    bool is_right_most = se.node->right_most_child == it->first;
                    trie_nodes_stack.push({is_right_most, it->first, &it->second, se.level+1});
                }
                prev_level = se.level;
            }


        }

        void _xbw_construction(const trie<value_type, word_type>& t, value_type leaf, nodes_type &nodes) {
            m_leaf = leaf;
            _build_nodes(t, nodes);
            std::cout << "Nodes built" << std::endl;
            std::cout << "nodes_size: " << nodes.size() << std::endl;
            std::sort(nodes.begin(), nodes.end());
            std::cout << "Nodes sorted" << std::endl;

#if VERBOSE
            std::cout << "nodes_size: " << nodes.size() << std::endl;
            for(uint64_t i =0; i < nodes.size(); ++i){
                //nodes[i].print();
                std::cout << "right: " << nodes[i].is_right_most << " alpha: " << nodes[i].alpha << " pi: [";
                for(uint64_t j = 0; j < nodes[i].pi.size(); ++j){
                    std::cout << nodes[i].pi[j] << ",";
                }
                std::cout << "]" << std::endl;
            }
            std::cout << std::endl;
#endif
            sdsl::int_vector<> alpha_int(nodes.size());
            t_bv last_aux(nodes.size(), 0);
            t_bv a_aux(nodes.size()+1, 0);
            sdsl::util::assign(m_last, last_aux);
            //sdsl::util::assign(m_a, a_aux);

            a_aux[0] = 1;
            a_aux[nodes.size()]=1;
            bool empty = true;
            value_type prev_value;
            for(size_type i = 0; i < nodes.size(); ++i){
                alpha_int[i] = nodes[i].alpha;
                if(nodes[i].is_right_most) m_last[i] = 1;
                if(!nodes[i].pi.empty()){
                    if(empty){
                        a_aux[i]=1;
                    }else if(prev_value != nodes[i].pi_first()){
                        a_aux[i]=1;
                    }
                    prev_value = nodes[i].pi_first();
                    empty = false;
                }
            }

            std::cout << "El" << std::endl;

            m_a = sdsl::sd_vector<>(a_aux);
            sdsl::construct_im(m_alpha, alpha_int);
            sdsl::util::init_support(m_rank1_last, &m_last);
            sdsl::util::init_support(m_select1_last, &m_last);
            sdsl::util::init_support(m_select1_a, &m_a);
            sdsl::util::init_support(m_rank1_a, &m_a);

#if VERBOSE
            print_array(m_alpha, "alpha");
            print_array(m_last, "last");
            print_array(m_a, "a");
#endif
        }

    public:
        xbw(){};

        xbw(const trie<value_type, word_type>& t, const value_type leaf){
            nodes_type nodes;
            _xbw_construction(t, leaf, nodes);
            nodes.clear();

        }

        xbw(const trie<value_type, word_type>& t, const value_type leaf, nodes_type &nodes){
            _xbw_construction(t, leaf, nodes);
        }

        bool is_leaf(const size_type n) const{
            return m_leaf == m_alpha[n];
        }

        void get_children(const size_type n, size_type &ini, size_type &end) const{
            if(is_leaf(n)){
                end = 0;
                ini = 1;
                return;
            }
            value_type c = m_alpha[n];
            //std::cout << "c: " << c << std::endl;
            size_type r = m_alpha.rank(n+1, c);
            //std::cout << "r: " << r << std::endl;
            size_type y = m_select1_a(c);
            //std::cout << "y: " << y<< std::endl;
            size_type z = m_rank1_last(y);
            //std::cout << "z: " << z << std::endl;
            ini = m_select1_last(z+r-1)+1;
            end = m_select1_last(z+r);
        }

        size_type get_ranked_child(const size_type n, const size_type k) const {
            size_type ini, end;
            get_children(n, ini, end);
            if(k > end-ini+1) return (size_type)-1;
            return ini+k-1;
        }

        size_type get_char_ranked_child(const size_type n, const value_type c, const size_type k) const {
            size_type ini, end;
            get_children(n, ini, end);
            if(ini > end) return (size_type)-1;
            size_type y1 = m_alpha.rank(ini, c);
            size_type y2 = m_alpha.rank(end+1, c);
            if(k > y2-y1) return (size_type)-1;
            return m_alpha.select(y1+k, c);
        }

        size_type get_degree(const size_type n) const{
            size_type ini, end;
            get_children(n, ini, end);
            return end-ini+1;
        }

        size_type get_char_degree(const size_type n, const value_type c) const{
            size_type ini, end;
            get_children(n, ini, end);
            if(ini > end) return 0;
            return m_alpha.rank(end+1, c) - m_alpha.rank(ini, c);
        }

        size_type get_parent(const size_type n){
            size_type c = m_rank1_a(n+1);
            if(c <= 1) return (size_type)-1;
            size_type y = m_select1_a(c);
            size_type k = m_rank1_last(n)-m_rank1_last(y);
            size_type p = m_alpha.select(k+1, c);
            return p;
        }

        void sub_path_search(const word_type &pattern, size_type &ini, size_type &end) const {
            if(pattern.size()<=1){
                ini = 0;
                end = m_alpha.size()-1;
                return;
            }
            //Range of pattern[0]
            if(pattern[0] > m_alpha.sigma){
                ini = 1;
                end = 0;
                return;
            }
            size_type first = m_select1_a(pattern[0]);
            size_type last = m_select1_a(pattern[0]+1)-1;
            for(size_type i = 1; i < pattern.size()-1; ++i){
                //Number of pattern[i] before first and until last
                size_type k1 = m_alpha.rank(first, pattern[i]);
                size_type k2 = m_alpha.rank(last+1, pattern[i]);
                std::cout << "k1: " << k1 << std::endl;
                std::cout << "k2: " << k2 << std::endl;
                if(k1 == k2){ //Not found
                    ini = 1;
                    end = 0;
                    return;
                }
                //Pos in alpha of the k-th pattern[i]
                size_type z1 = m_alpha.select(k1+1, pattern[i]);
                size_type z2 = m_alpha.select(k2, pattern[i]);
                first = get_ranked_child(z1, 1); //first child of z1
                last = get_ranked_child(z2, get_degree(z2));//last child of z2
            }
            ini = first;
            end = last;
        }


        std::vector<size_type> sub_path_search(const word_type &pattern) const {
            std::vector<size_type> ret;
            size_type ini, end;
            sub_path_search(pattern, ini, end);
            if(ini > end){
                return ret;
            }
            size_type c_last = pattern[pattern.size()-1];
            size_type k1 = m_alpha.rank(ini, c_last);
            size_type k2 = m_alpha.rank(end+1, c_last);
            std::cout << "k1: " << k1 << std::endl;
            std::cout << "k2: " << k2 << std::endl;
            size_type len = k2-k1;
            ret.resize(len);
            for(size_type i = k1+1; i <= k2; ++i){
                ret[i-k1-1] = m_alpha.select(i,c_last);
            }
            return ret;
        }


        xbw& operator=(xbw &&p){
            if (this != &p){
                m_alpha = std::move(p.m_alpha);
                m_last = std::move(p.m_last);
                m_a = std::move(p.m_a);
                m_rank1_last = std::move(p.m_rank1_last);
                m_rank1_last.set_vector(&m_last);
                m_select1_last = std::move(p.m_select1_last);
                m_select1_last.set_vector(&m_last);
                m_select1_a = std::move(p.m_select1_a);
                m_select1_a.set_vector(&m_a);
                m_rank1_a = std::move(p.m_rank1_a);
                m_rank1_a.set_vector(&m_a);
                m_leaf = std::move(p.m_leaf);
            }
            return *this;
        }

        //! Assignment operator
        xbw& operator=(const xbw& p)
        {
            if (this != &p) {
                m_alpha = p.m_alpha;
                m_last = p.m_last;
                m_a = p.m_a;
                m_leaf = p.m_leaf;
                m_rank1_last = p.m_rank1_last;
                m_rank1_last.set_vector(&m_last);
                m_select1_last = p.m_select1_last;
                m_select1_last.set_vector(&m_last);
                m_select1_a = p.m_select1_a;
                m_select1_a.set_vector(&m_a);
                m_rank1_a = p.m_rank1_a;
                m_rank1_a.set_vector(&m_a);
            }
            return *this;
        }

        //! Copy constructor
        xbw(const xbw& p)
        {
            m_alpha = p.m_alpha;
            m_last = p.m_last;
            m_a = p.m_a;
            m_leaf = p.m_leaf;
            m_rank1_last = p.m_rank1_last;
            m_rank1_last.set_vector(&m_last);
            m_select1_last = p.m_select1_last;
            m_select1_last.set_vector(&m_last);
            m_select1_a = p.m_select1_a;
            m_select1_a.set_vector(&m_a);
            m_rank1_a = p.m_rank1_a;
            m_rank1_a.set_vector(&m_a);
        }

        //! Move constructor
        xbw(xbw&& p)
        {
            *this = std::move(p);
        }

        //! Swap method
        void swap(xbw& p)
        {
            if(this != &p){
                m_last.swap(p.m_last);
                m_alpha.swap(p.m_alpha);
                m_a.swap(p.m_a);
                m_rank1_last.swap(p.m_rank1_last);
                m_select1_last.swap(p.m_select1_last);
                m_rank1_a.swap(p.m_rank1_a);
                m_select1_a.swap(p.m_select1_a);
                std::swap(m_leaf, p.m_leaf);
            }
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr,
                            std::string name="") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;

            written_bytes += m_last.serialize(out, child, "last");
            written_bytes += m_alpha.serialize(out, child, "alpha");
            written_bytes += m_a.serialize(out, child, "a");
            written_bytes += m_rank1_last.serialize(out, child, "rank_last");
            written_bytes += m_select1_last.serialize(out, child, "select_last");
            written_bytes += m_rank1_a.serialize(out, child, "rank_a");
            written_bytes += m_select1_a.serialize(out, child, "select_a");
            written_bytes += write_member(m_leaf, out, child, "leaf");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }


    };
}

#endif //XBW_TRIE_XBW_HPP
