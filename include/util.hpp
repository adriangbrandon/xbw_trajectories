//
// Created by adrian on 2/06/18.
//

#ifndef XBW_TRIE_UTIL_HPP
#define XBW_TRIE_UTIL_HPP

#include <unordered_map>
#include <vector>

namespace cnt {

    namespace util{

        class n_gram {

        public:
            typedef std::vector<uint32_t> data_type;
        private:
            uint64_t m_n = 0;
            data_type  m_values;
        public:
            const uint64_t &n = m_n;
            n_gram(){};

            n_gram(data_type array, uint64_t n){
                m_values = array;
                m_n = n;
            }

            n_gram& operator=(n_gram &&p){
                if (this != &p){
                    m_values = std::move(p.m_values);
                    m_n = std::move(p.m_n);
                }
                return *this;
            }

            //! Assignment operator
            n_gram& operator=(const n_gram& p)
            {
                if (this != &p) {
                    m_values = p.m_values;
                    m_n = p.m_n;
                }
                return *this;
            }

            //! Copy constructor
            n_gram(const n_gram& p)
            {
                m_values = p.m_values;
                m_n = p.m_n;
            }

            //! Move constructor
            n_gram(n_gram&& p)
            {
                *this = std::move(p);
            }

            //! Swap method
            void swap(n_gram& p)
            {
                if(this != &p){
                    std::swap(m_values, p.m_values);
                    std::swap(m_n, p.m_n);
                }
            }

            uint32_t operator[](uint64_t i) const {
                return m_values[i];
            }
        };

        class n_gram_hash {
        public:
            std::size_t operator()(const n_gram& k) const
            {

                uint64_t result = 0;
                for(uint64_t i = 0; i < k.n; i++){
                    result = result ^ (std::hash<int32_t>()(k[i]) << i);
                }
                return result;
            }
        };

        class n_gram_equal {
        public:
            bool operator()(const n_gram& lhs, const n_gram& rhs) const
            {
                if(lhs.n != rhs.n) return false;
                for(uint64_t i = 0; i < lhs.n; ++i){
                    if(lhs[i] != rhs[i]) return false;
                }
                return true;
            }
        };

        static struct {
            bool operator()(std::pair<util::n_gram , uint32_t> a, std::pair<util::n_gram , uint32_t> b) const
            {
                return a.second > b.second;
            }
        } compare_by_frequency;

        typedef std::unordered_map<uint32_t , char> map_node_type;

        typedef std::unordered_map<util::n_gram, uint32_t, util::n_gram_hash, util::n_gram_equal> map_n_gram_type;
        typedef std::vector<std::pair<util::n_gram , uint32_t >> keys_n_gram_type;
    }


}

#endif //CINCT_UTIL_HPP
