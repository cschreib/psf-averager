#ifndef TREE_INCLUDED
#define TREE_INCLUDED

#include <vif.hpp>

using namespace vif;

template<typename KeyType, typename DataType>
class tree {
    vec<1,vec<1,KeyType>> bins_low;
    vec<1,vec<1,KeyType>> bins_up;

public :
    using leaf = DataType;
    struct branch;

    struct node {
        uint_t id = npos;

        explicit node(uint_t i) : id(i) {}

        bool operator < (const node& n) const {
            return id < n.id;
        }

        bool operator < (uint_t i) const {
            return id < i;
        }

        void print(std::string ind = "") {
            if (b) {
                vif::print(ind, "- ", id);
                for (auto& n : b->nodes) {
                    n.print(ind+"  ");
                }
            } else if (l) {
                vif::print(ind, "- ", id, ": ", l->data, " (", l->weight, ")");
            }
        }

        std::unique_ptr<branch> b;
        std::unique_ptr<leaf> l;
    };

    using ElemType = meta::data_type_t<DataType>;
    static const uint_t DataDim = meta::vec_dim<DataType>::value;

    struct branch {
        std::vector<node> nodes;
    };

    node root;
    uint_t nleaf = 0;
    uint_t nnode = 0;
    uint_t depth = 0;

    vec1u bin_key(const vec<1,KeyType>& key) const {
        vec1u k;
        k.resize(depth);

        for (uint_t i : range(depth)) {
            auto& b = bins_up[i];
            uint_t p = upper_bound(b, key.safe[i]);
            if (p == npos) p = b.size()-1;
            k.safe[i] = p;
        }

        return k;
    }

    explicit tree() : root(0) {}

    void reset_() {
        if (depth == 0) {
            root.l = std::unique_ptr<leaf>(new leaf);
            ++nleaf;
        } else {
            root.b = std::unique_ptr<branch>(new branch);
        }
    }

    void setup(const vec<1,vec<2,KeyType>>& tbins) {
        clear();

        depth = tbins.size();

        bins_low.resize(depth);
        bins_up.resize(depth);
        for (uint_t i : range(depth)) {
            bins_low[i] = tbins[i](0,_);
            bins_up[i] = tbins[i](1,_);
        }

        reset_();
    }

    leaf& insert_binned(const vec1u& k) {
        // Locate (and if needed, create) the leaf corresponding to this object
        node* n = &root;
        for (uint_t i : range(depth)) {
            branch& b = *n->b;

            auto p = std::lower_bound(b.nodes.begin(), b.nodes.end(), k.safe[i]);
            if (p == b.nodes.end() || p->id != k.safe[i]) {
                p = b.nodes.insert(p, node(k.safe[i]));
                ++nnode;

                if (i != depth-1) {
                    p->b = std::unique_ptr<branch>(new branch);
                } else {
                    // Leaf didn't exist, create it and return
                    p->l = std::unique_ptr<leaf>(new leaf);
                    ++nleaf;
                    return *p->l;
                }
            }

            n = &(*p);
        }

        // Found the leaf, return it
        return *n->l;
    }

    leaf& insert(const vec<1,KeyType>& key) {
        // Locate (and if needed, create) the leaf corresponding to this object
        return insert_binned(bin_key(key));
    }

    template<typename F>
    void for_each_leaf_(const F& f, branch& b, vec1u& key, uint_t d) {
        if (d == depth-1) {
            for (auto& n : b.nodes) {
                key[d] = n.id;
                f(key, *n.l);
            }
        } else {
            for (auto& n : b.nodes) {
                key[d] = n.id;
                for_each_leaf_(f, *n.b, key, d+1);
            }
        }
    }

    template<typename F>
    void for_each_leaf(const F& f) {
        if (depth == 0) {
            vec1u key;
            f(key, *root.l);
        } else {
            vec1u key(depth);
            for_each_leaf_(f, *root.b, key, 0);
        }
    }

    const leaf& first_leaf() const {
        // Locate the first leaf in the tree
        const node* n = &root;
        for (uint_t i : range(depth)) {
            auto& b = *n->b;
            vif_check(!b.nodes.empty(), "cannot query the first leaf in an empty tree");
            n = &b.nodes[0];
        }

        const leaf& l = *n->l;
        return l;
    }

    void clear() {
        if (depth == 0) {
            root.l = nullptr;
        } else {
            root.b = nullptr;
        }

        nleaf = 0;
        nnode = 1;
    }

    uint_t size() const {
        return nleaf;
    }
};

#endif
