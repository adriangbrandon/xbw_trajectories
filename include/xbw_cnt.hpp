//
// Created by adrian on 3/07/18.
//

#ifndef XBW_TRIE_XBW_CNT_HPP
#define XBW_TRIE_XBW_CNT_HPP

#include <graph.hpp>
#include <trie.hpp>
#include <xbw.hpp>
#include <util.hpp>
#include <sdsl/util.hpp>
#include <sdsl/int_vector.hpp>


namespace xbw_trie {


    template<class t_value = uint32_t,
             class t_word = std::vector<t_value>,
             class t_wt = sdsl::wt_blcd_int< sdsl::rrr_vector<>>,
             class t_bv = sdsl::bit_vector>
    class xbw_cnt {

    public:

        typedef t_value value_type;
        typedef t_word word_type;
        typedef uint64_t size_type;
        typedef xbw<t_value, t_word, t_wt, t_bv> xbw_type;

    private:

        //int_vector<> m_C;
        value_type m_leaf;
        size_type m_sigma;
        cnt::graph m_graph;
        t_wt m_alpha;
        t_bv m_last;
        sdsl::rank_support_v5<> m_rank1_last;
        sdsl::select_support_mcl<> m_select1_last;
        sdsl::sd_vector<> m_a;
        sdsl::sd_vector<>::select_1_type m_select1_a;
        sdsl::sd_vector<>::rank_1_type m_rank1_a;
        //wt_huff_int<rrr_vector<> > m_wt_label_bwt;

    public:
        const cnt::graph &graph = m_graph;
        const t_wt &alpha = m_alpha;

    private:

        void _construct_xbw(const std::string &in, xbw_type &xbwa){
            uint32_t x, orig;
            std::ifstream input(in);
            trie<> trie;
            word_type w;
            while(input){
                input.read((char*) &x, sizeof(t_value));
                if(input.eof()) break;
                w.push_back(x);
                if(x == 1){
                    trie.insert(w);
                    w.clear();
                }
            }
            input.close();
            xbwa =  xbw_type(trie, 1);
            sdsl::util::clear(trie);
        }

        void _construct_graph(xbw_type &xbwa){

            cnt::util::map_n_gram_type::iterator it, next_it;
            cnt::util::map_n_gram_type n_gram_freq;
            uint64_t traj_i = 0, max = 0;
            std::vector<uint32_t > v(2);
            uint32_t first = 0;
            //1. Collect the frequencies of the n-grams
            value_type prev_value = 0;
            for(size_type i = 0; i < xbwa.alpha.size(); ++i){
                if(xbwa.a[i]) ++prev_value;
                v[0]=prev_value;
                v[1]=xbwa.alpha[i];
                cnt::util::n_gram n_g(v,2);
                if((it = n_gram_freq.find(n_g)) != n_gram_freq.end()){
                    ++it->second;
                }else{
                    n_gram_freq.insert(std::pair<cnt::util::n_gram, uint32_t>(n_g, 1));
                }
            }

            //2. Store the frequencies in an array
            it = n_gram_freq.begin();
            cnt::util::keys_n_gram_type keys;
            while(it != n_gram_freq.end()){
                keys.push_back(*it);
                next_it = std::next(it);
                n_gram_freq.erase(it);
                it = next_it;
            }
            n_gram_freq.clear();

            //3. Sort the frequencies
            std::sort(keys.begin(), keys.end(), cnt::util::compare_by_frequency);

            //4. Build the graph
            std::vector<std::vector<value_type>> adj_list(prev_value+1, std::vector<value_type>());
            for(size_type i = 0; i < keys.size(); ++i){
                uint32_t index = keys[i].first[0];
                uint32_t value = keys[i].first[1];
                adj_list[index].push_back(value);
            }
            keys.clear();
            m_graph = cnt::graph(adj_list);

        }


        void _label_xbw(xbw_type &xbwa){
            value_type prev_value = 0, label;
            sdsl::int_vector<> alpha_int_vector(xbwa.alpha.size());
            for(size_type i = 0; i < xbwa.alpha.size(); ++i){
                if(xbwa.a[i]) ++prev_value;
                label = m_graph.search_label(prev_value, xbwa.alpha[i]);
                alpha_int_vector[i] = label;
            }
            sdsl::construct_im(m_alpha, alpha_int_vector);
        }

        void _add_z_graph(const t_wt &aux_wt){
            for(value_type w_i = 0; w_i < m_graph.size(); ++w_i){
                for(size_type i = 0; i < m_graph[w_i].size(); ++i){
                    value_type w = m_graph[w_i][i].value;
                    size_type p_c = m_select1_a(w_i);
                    int32_t z = m_alpha.rank(p_c, i+1) - aux_wt.rank(p_c, w);
                    m_graph.set_z_at(i, w_i, z);
                    //std::cout << "z: " << z << std::endl;
                }
            }
        }

        size_type _pseudo_rank(const size_type pos, const value_type w, const value_type w_i) const{
            size_type label = m_graph.pos(w, w_i);
            if(label > 0 && m_select1_a(w_i) <= pos && pos <= m_select1_a(w_i+1)){
                return m_alpha.rank(pos, label) - m_graph[w_i][label-1].z;
            }
            return 0;
        }

        size_type _pseudo_rank_exist(const size_type pos, const size_type label, const value_type w_i){
            if(label > 0 && m_select1_a(w_i) <= pos && pos <= m_select1_a(w_i+1)){
                return m_alpha.rank(pos, label) - m_graph[w_i][label-1].z;
            }
            return 0;
        }

    public:

        xbw_cnt(){};

        xbw_cnt(const std::string &traj){
            xbw_type xbwa;
            _construct_xbw(traj, xbwa);
            _construct_graph(xbwa);
            m_last = xbwa.last;
            m_rank1_last = xbwa.rank1_last;
            m_select1_last = xbwa.select1_last;
            m_rank1_last.set_vector(&m_last);
            m_select1_last.set_vector(&m_last);
            m_a = xbwa.a;
            m_rank1_a = xbwa.rank1_a;
            m_select1_a = xbwa.select1_a;
            m_rank1_a.set_vector(&m_a);
            m_select1_a.set_vector(&m_a);
            m_leaf = xbwa.leaf;
            m_sigma = xbwa.alpha.sigma;
            _label_xbw(xbwa);
            _add_z_graph(xbwa.alpha);
            sdsl::util::clear(xbwa);
        }


        bool is_leaf(const size_type n) const{
            return m_leaf == m_alpha[n];
        }

        bool is_leaf_value(const value_type v) const {
            return m_leaf == v;
        }

        void get_children(const size_type n, const value_type w_i, size_type &ini, size_type &end) const{

            value_type label = m_alpha[n];
            value_type c = m_graph[w_i][label-1].value;
            if(is_leaf_value(c)){
                end = 0;
                ini = 1;
                return;
            }
            //std::cout << "c: " << c << std::endl;
            size_type r = _pseudo_rank(n+1, c, w_i);
            //size_type r = m_alpha.rank(n+1, c);
            //std::cout << "r: " << r << std::endl;
            size_type y = m_select1_a(c);
            //std::cout << "y: " << y<< std::endl;
            size_type z = m_rank1_last(y);
            //std::cout << "z: " << z << std::endl;
            ini = m_select1_last(z+r-1)+1;
            end = m_select1_last(z+r);
        }

        void get_children(const size_type n, size_type &ini, size_type &end) const{
            value_type w_i = m_rank1_a(n+1);
            get_children(n, w_i, ini, end);
        }

        size_type get_ranked_child(const size_type n, const value_type w_i, const size_type k) const {
            size_type ini, end;
            get_children(n, w_i, ini, end);
            if(k > end-ini+1) return (size_type)-1;
            return ini+k-1;
        }

        size_type get_char_ranked_child(const size_type n, const value_type w_i, const value_type c, const size_type k) const {
            size_type ini, end;
            get_children(n, w_i, ini, end);
            if(ini > end) return (size_type)-1;
            size_type y1 = m_alpha.rank(ini, c);
            size_type y2 = m_alpha.rank(end+1, c);
            if(k > y2-y1) return (size_type)-1;
            return m_alpha.select(y1+k, c);
        }

        size_type get_degree(const size_type n, const value_type w_i) const{
            size_type ini, end;
            get_children(n, w_i, ini, end);
            return end-ini+1;
        }

        size_type get_char_degree(const size_type n, const value_type w_i, const value_type c) const{
            size_type ini, end;
            get_children(n, w_i, ini, end);
            if(ini > end) return 0;
            return m_alpha.rank(end+1, c) - m_alpha.rank(ini, c);
        }


        //TODO: it doesn't work
        size_type get_parent(const size_type n){
            size_type c = m_rank1_a(n+1);
            if(c <= 1) return (size_type)-1;
            size_type y = m_select1_a(c);
            size_type k = m_rank1_last(n)-m_rank1_last(y);
            size_type p = m_alpha.select(k+1, c);
            return p;
        }

        void sub_path_search(const word_type &pattern, size_type &ini, size_type &end) const {
            //Range of pattern[0]
            if(pattern[0] > m_sigma){
                ini = 1;
                end = 0;
                return;
            }
            size_type first = m_select1_a(pattern[0]);
            size_type last = m_select1_a(pattern[0]+1)-1;
            size_type  len;
            for(size_type i = 1; i < pattern.size()-1; ++i){
                //Number of pattern[i] before first and until last
                size_type k1 = _pseudo_rank(first, pattern[i], pattern[i-1]);
                size_type k2 = _pseudo_rank(last+1, pattern[i], pattern[i-1]);
                //TODO:
               /* size_type label = m_graph.pos(pattern[i], pattern[i-1]);
                size_type k1 = m_alpha.rank(first, label);
                size_type k2 = m_alpha.rank(last+1, label);*/
                if(k1 == k2){ //Not found
                    ini = 1;
                    end = 0;
                    return;
                }
                len = k2-k1;
                //Pos in alpha of the k-th pattern[i]
                size_type label = m_graph.pos(pattern[i], pattern[i-1]);
                size_type rank_label = m_alpha.rank(first, label);
                size_type z1 = m_alpha.select(rank_label +1, label);
                size_type z2 = m_alpha.select(rank_label + len, label);
                first = get_ranked_child(z1, pattern[i-1], 1); //first child of z1
                last = get_ranked_child(z2, pattern[i-1], get_degree(z2, pattern[i-1]));//last child of z2
            }
            ini = first;
            end = last;
        }


        std::vector<size_type> sub_path_search(const word_type &pattern) const {
            std::vector<size_type> ret;
            size_type ini= 0, end = m_alpha.size()-1;
            if(pattern.size() == 1){
                value_type w_i = 1, w_i_last = m_sigma;
                value_type c = pattern[0];;
                size_type l_rank_c = 0, u_rank_c = 0;
                size_type pos_c = 0, pos_next_c = 0;
                //size_type pos_ret = 0;
                while(w_i <= w_i_last){
                    size_type label = m_graph.pos(c, w_i);
                    pos_next_c = m_select1_a(w_i+ 1);
                    l_rank_c = pos_c ? m_alpha.rank(pos_c, label) : 0;
                    u_rank_c = m_alpha.rank(pos_next_c, label);
                    for(size_type i = l_rank_c+1; i <= u_rank_c; ++i){
                       ret.push_back(m_alpha.select(i, label));
                    }
                    ++w_i;
                    pos_c = pos_next_c;
                }
            }else{
                sub_path_search(pattern, ini, end);
                if(ini > end) return ret;
                value_type w_i = pattern[pattern.size()-2];
                value_type c = pattern[pattern.size()-1];
                size_type k1 = _pseudo_rank(ini, c, w_i);
                size_type k2 = _pseudo_rank(end+1, c, w_i);
                size_type len = k2-k1;
                ret.resize(len);
                size_type label = m_graph.pos(c, w_i);
                size_type rank_label = m_alpha.rank(ini, label);
                for(size_type i = 1; i <= len; ++i){
                    ret[i-1] = m_alpha.select(rank_label + i,label);
                }
            }
            return ret;

        }

        std::vector<size_type> sub_path_search_v2(const word_type &pattern) const {
            std::vector<size_type> ret;
            size_type ini, end;
            sub_path_search(pattern, ini, end);
            if(ini > end){
                return ret;
            }
            value_type w_i = 1, w_i_last = m_sigma;
            if(pattern.size() > 1){
                w_i = pattern[pattern.size()-2];
            }else{
                std::cout << "ini: " << ini << std::endl;
                std::cout << "end: " << end << std::endl;
            }
            value_type c_last = pattern[pattern.size()-1];
            //size_type k1 = m_alpha.rank(ini, c_last);
            size_type k1 = _pseudo_rank(ini, c_last, w_i);
            size_type k2 = _pseudo_rank(end+1, c_last, w_i_last);
            size_type label = m_graph.pos(c_last, w_i);
            /*size_type k1 = m_alpha.rank(ini, label);
            size_type k2 = m_alpha.rank(end+1, label);*/
            size_type len = k2-k1;
            ret.resize(len);
            std::cout
            << "k1: " << k1 << std::endl;
            std::cout << "k2: " << k2 << std::endl;
            size_type rank_label = m_alpha.rank(ini, label);
            std::cout << "rank_label: " << rank_label << std::endl;
            std::cout << "label: " << label << std::endl;
            /*size_type z1 = m_alpha.select(rank_label + k1+1, label);
            size_type z2 = m_alpha.select(rank_label + k2, label);*/
            for(size_type i = 1; i <= len; ++i){
            //for(size_type i = k1+1; i <= k2; ++i){
                //std::cout << m_alpha.select(rank_label + i,label) << std::endl;
                ret[i-1] = m_alpha.select(rank_label + i,label);
            }
            return ret;
        }


        xbw_cnt& operator=(xbw_cnt &&p){
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
                m_sigma = std::move(p.m_sigma);
                m_graph = std::move(p.m_graph);
            }
            return *this;
        }

        //! Assignment operator
        xbw_cnt& operator=(const xbw_cnt& p)
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
                m_sigma = p.m_sigma;
                m_graph = p.m_graph;
            }
            return *this;
        }

        //! Copy constructor
        xbw_cnt(const xbw_cnt& p)
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
            m_sigma = p.m_sigma;
            m_graph = p.m_graph;
        }

        //! Move constructor
        xbw_cnt(xbw_cnt&& p)
        {
            *this = std::move(p);
        }

        //! Swap method
        void swap(xbw_cnt& p)
        {
            if(this != &p){
                m_last.swap(p.m_last);
                m_alpha.swap(p.m_alpha);
                m_a.swap(p.m_a);
                m_rank1_last.swap(p.m_rank1_last);
                m_select1_last.swap(p.m_select1_last);
                m_rank1_a.swap(p.m_rank1_a);
                m_select1_a.swap(p.m_select1_a);
                m_graph.swap(p.m_graph);
                std::swap(m_sigma, p.m_sigma);
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
            written_bytes += m_graph.serialize(out, child, "graph");
            written_bytes += write_member(m_leaf, out, child, "leaf");
            written_bytes += write_member(m_sigma, out, child, "sigma");

            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
    };
}

#endif //XBW_TRIE_XBW_CNT_HPP
