//
// Created by adrian on 23/06/18.
//

#ifndef XBW_TRIE_TRIE_HPP
#define XBW_TRIE_TRIE_HPP

#include <vector>
#include <map>
#include <queue>
#include <iostream>

namespace xbw_trie {


    /**
     * Node of a trie
     */
    template <class value_t = uint32_t,
              class word_t = std::vector<value_t>>
    class trie_node {


    public:

        typedef value_t value_type;
        typedef uint64_t size_type;
        typedef std::map<value_type, trie_node> children_type;
        typedef typename children_type::iterator iterator;

    private:
        value_type m_right_most_child = 0;
        children_type m_children;

    public:
        const children_type &children = m_children;

        const value_t &right_most_child = m_right_most_child;

        trie_node(){};

        iterator insert_child(const value_type c){
            iterator it;
            if((it = m_children.find(c)) == m_children.end()){
                auto ret = m_children.insert({c, trie_node()});
                it = ret.first;
                if(c > m_right_most_child) m_right_most_child = c;

            }
            return it;
        }

        trie_node& operator=(trie_node &&p){
            if (this != &p){
                m_children = std::move(p.m_children);
                m_right_most_child = std::move(p.m_right_most_child);
            }
            return *this;
        }

        //! Assignment operator
        trie_node& operator=(const trie_node& p)
        {
            if (this != &p) {
                m_children = p.m_children;
                m_right_most_child = p.m_right_most_child;
            }
            return *this;
        }

        //! Copy constructor
        trie_node(const trie_node& p)
        {
            m_children = p.m_children;
            m_right_most_child = p.m_right_most_child;
        }

        //! Move constructor
        trie_node(trie_node&& p)
        {
            *this = std::move(p);
        }

        //! Swap method
        void swap(trie_node& p)
        {
            if(this != &p){
                std::swap(m_children, p.m_children);
                std::swap(m_right_most_child, p.m_right_most_child);
            }
        }
    };


    /*
     * Trie
     */
    template <class value_t = uint32_t,
              class word_t = std::vector<value_t>>
    class trie {

    public:
        typedef value_t value_type;
        typedef word_t word_type;
        typedef uint64_t size_type;

    private:

        trie_node<value_type, word_type> m_root;

    public:
        const trie_node<value_type, word_type> &root = m_root;

        trie(){}
        void insert(const word_type &word){

            //1. Start on the root
            trie_node<value_type, word_type>* actual_node = &m_root;
            for(size_type i = 0; i < word.size(); ++i){
                auto it = actual_node->insert_child(word[i]);
                actual_node = &(it->second);
            }
        }



        trie& operator=(trie &&p){
            if (this != &p){
                m_root = std::move(p.m_root);
            }
            return *this;
        }

        //! Assignment operator
        trie& operator=(const trie& p)
        {
            if (this != &p) {
                m_root = p.m_root;
            }
            return *this;
        }

        //! Copy constructor
        trie(const trie& p)
        {
            m_root = p.m_root;
        }

        //! Move constructor
        trie(trie&& p)
        {
            *this = std::move(p);
        }

        //! Swap method
        void swap(trie& p)
        {
            if(this != &p){
                std::swap(m_root, p.m_root);
            }
        }

        void print(){
            std::queue<const trie_node<value_type, word_type >*> queue;
            queue.push(&m_root);
            while(!queue.empty()){
                auto node = queue.front();
                queue.pop();
                for(auto it = node->children.begin(); it != node->children.end(); ++it){
                    std::cout << it->first << " ";
                    queue.push(&(it->second));
                }
            }
            std::cout << std::endl;
        }



    };

}

#endif //XBW_TRIE_TRIE_HPP
