#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>
#include <queue>
#include <iostream>
#include <memory>
#include <utility>
#include <cassert>
#include <unordered_set>

class Graph
{
    typedef long long Distance;
    typedef int Vertex;

    
    int N;  // Number of nodes
    int s, t;   // Source and target
    static constexpr int INFINITY = std::numeric_limits<int>::max() / 2;
    Distance estimate = INFINITY;    // Estimate of the distance from s to t
    std::vector<std::vector<std::pair<Vertex, Distance>>> outgoing_edges;   // Lists of edges outgoing from each node 
    std::vector<std::vector<std::pair<Vertex, Distance>>> incoming_edges;   // Lists of edges incoming to each node
    std::vector<Vertex> level;     // Levels of nodes for node ordering
    std::vector<Vertex> rank;      // Ranks of nodes - positions in the node ordering
    std::vector<std::vector<Distance>> bidistance;      // Distance to node v, bidistance[0][v] - forward search, bidistance[1][v] - backward search.
    std::vector<bool> contracted;       // record of nodes that have been contracted
    std::vector<bool> witness_found;        // witness path found or not
    std::vector<bool> shortcut_cover;       // shortcut covers added
    std::vector<Vertex> hops;    // number of hops
    std::unordered_set<Vertex> workset;      // unordered set to keep track of nodes affected

    // Wrapper around STL priority_queue
    class StlHeap
    {
    public:
        using T = std::pair<Distance, Vertex>;
        using Queue = std::priority_queue<T, std::vector<T>, std::greater<T>>;

        StlHeap() {
            queue.reset(new Queue());
        }

        bool empty() const {
            return queue->empty();
        }

        void update(Vertex v, Distance d) {
            queue->push(std::make_pair(d,v));
        }

        void clear() {
            queue.reset(new Queue());
        }

        std::pair<Distance, Vertex> pop() {
            std::pair<Distance, Vertex> top = queue->top();
            queue->pop();
            return top;
        }

    private:
        std::unique_ptr<Queue> queue;
    };

    // Priority queues for forward and backward searches
    StlHeap diqueue[2];
    // std::vector<std::priority_queue<std::pair<long long, long long>, std::vector<std::pair<long long, long long>>, std::greater<std::pair<long long, long long>>>;
public:
    Graph() {
        read_stdin();
        bidistance.resize(2, std::vector<long long>(N, INFINITY));
        contracted.resize(N, false);
        witness_found.resize(N, true);
        shortcut_cover.resize(N, false);
        hops.resize(N, 0);
        level.resize(N, 0);
        rank.resize(N, 0);
    }

    void workset_clear(bool clear_rev = true) {

        for (Vertex node : workset) {
            bidistance[0][node] = INFINITY ;
            if (!clear_rev) { hops[node] = 0; shortcut_cover[node] = false; witness_found[node] = true; }
            if (clear_rev) { bidistance[1][node] = INFINITY ; }
        }
        workset.clear();
        diqueue[0].clear();
        if (clear_rev) diqueue[1].clear();
    }

    void dijkstra(Vertex start_node, Vertex avoid_node, Distance max_cost) {
        workset_clear(false);
        bidistance[0][start_node] = 0;
        diqueue[0].update(start_node, 0);
        workset.insert(start_node);
        std::pair<Distance, Vertex> pair;
        while(!diqueue[0].empty()) {
            pair = diqueue[0].pop();
            if (bidistance[0][pair.second] < pair.first || hops[pair.second] >= 5) continue;
            if (bidistance[0][pair.second] > max_cost) break;
            for (std::pair<Vertex, Distance>& node : outgoing_edges[pair.second]) {
                if (contracted[node.first] || node.first == avoid_node) continue;
                if (bidistance[0][node.first] > bidistance[0][pair.second] + node.second) {
                    workset.insert(node.first);
                    bidistance[0][node.first] = bidistance[0][pair.second] + node.second;
                    hops[node.first] = hops[pair.second] + 1;
                    if (hops[node.first] < 5 && bidistance[0][node.first] < max_cost) { diqueue[0].update(node.first, bidistance[0][node.first]); }
                    if (!witness_found[node.first]) { witness_found[node.first] = true; }
                }
            }
        }
    }

    void contract(Vertex node) {
        contracted[node] = true;
        workset_clear(false);
        long long max_out = 0;
        for (std::pair<Vertex, Distance>& pair_out : outgoing_edges[node]) { 
            if (contracted[pair_out.first]) { continue; }
            max_out = std::max(max_out, pair_out.second); 
            witness_found[pair_out.first] = false;
            workset.insert(pair_out.first);
        }
        
        for (std::pair<Vertex, Distance>& pair_in : incoming_edges[node]) {
            if (contracted[pair_in.first]) { continue; }
            workset.insert(pair_in.first);
            dijkstra(pair_in.first, node, pair_in.second + max_out); 
            for (std::pair<Vertex, Distance>& pair_out : outgoing_edges[node]) {
                if (contracted[pair_out.first]) { continue; }
                Distance edge_cost = pair_in.second + pair_out.second;
                if (bidistance[0][pair_out.first] > edge_cost) { 
                    outgoing_edges[pair_in.first].emplace_back(std::make_pair(pair_out.first, edge_cost));
                    incoming_edges[pair_out.first].emplace_back(std::make_pair(pair_in.first, edge_cost));
                }
            }
        }
    }

    int compute_importance(Vertex node) {  
        Distance max_out = 0;
        int imp = 0, ed = -1*incoming_edges[node].size() - outgoing_edges[node].size(), cn = 0, sc = 0;
        workset_clear(false);
        for (std::pair<Vertex, Distance>& pair_out : outgoing_edges[node]) {   
            if (contracted[pair_out.first]) { cn += 1; level[node] = std::max(level[node], level[pair_out.first] + 1); continue; }
            max_out = std::max(max_out, pair_out.second); 
            witness_found[pair_out.first] = false;
            shortcut_cover[pair_out.first] = true;
        }
        
        for (std::pair<Vertex, Distance>& pair_in : incoming_edges[node]) {
            if (contracted[pair_in.first]) { cn += 1; level[node] = std::max(level[node], level[pair_in.first] + 1); continue; }
            shortcut_cover[pair_in.first] = true;
            dijkstra(pair_in.first, node, pair_in.second + max_out); 
            for (std::pair<Vertex, Distance>& pair_out : outgoing_edges[node]) {
                if (contracted[pair_out.first]) { continue; }
                if (!witness_found[pair_out.first]) { 
                    ed += 1;
                    if (shortcut_cover[pair_in.first]) { sc += 1; shortcut_cover[pair_in.first] = false; }
                    if (shortcut_cover[pair_out.first]) { sc += 1; shortcut_cover[pair_out.first] = false; }
                }
                else { 
                    witness_found[pair_out.first] = false; 
                    if (bidistance[0][pair_out.first] > pair_in.second + pair_out.second) { 
                        ed += 1; 
                        if (shortcut_cover[pair_in.first]) { sc += 1; shortcut_cover[pair_in.first] = false; }
                        if (shortcut_cover[pair_out.first]) { sc += 1; shortcut_cover[pair_out.first] = false; }
                    }
                }
            }
        }
        return (ed + 2*cn + 3*sc + 1*level[node]); // heuristics
    }

    void preprocess() {
        std::priority_queue<std::pair<int, Vertex>, std::vector<std::pair<int, Vertex>>, std::greater<std::pair<int, Vertex>>> queue;

        // Implement the rest of the algorithm yourself
        for (int i = 0; i < N; i++) {
            queue.emplace(std::make_pair(compute_importance(i), i));
        }
        int rnk = 1;
        int new_imp;
        std::pair<int, Vertex> selected_node;
        while(!queue.empty()) {
            selected_node = queue.top(); // min imp is preferred
            queue.pop();
            if (queue.size() != N-1) { // if size is N-1, it is the first iteration
                new_imp = compute_importance(selected_node.second);
                if (selected_node.first < new_imp) {
                    queue.emplace(std::make_pair(new_imp, selected_node.second)); 
                    continue;
                } 
            }
            contract(selected_node.second);
            rank[selected_node.second] = rnk++;
        }
    }
    
    // Returns distance from s to t in the graph
    Distance query(Vertex u, Vertex w) {

        // Implement the rest of the algorithm yourself
        workset_clear(true);
        if (u == w) return 0;
        diqueue[0].update(u, 0);
        diqueue[1].update(w, 0);
        bidistance[0][u] = 0;
        bidistance[1][w] = 0;
        workset.insert(u);
        workset.insert(w);
        std::pair<Distance, Vertex> pair;
        estimate = INFINITY ;

        while(!diqueue[0].empty() || !diqueue[1].empty()) {

            // forward graph search
            if (!diqueue[0].empty()) {
                pair = diqueue[0].pop();
                
                if (bidistance[0][pair.second] == pair.first) {
                    if (pair.first <= estimate) {
                        update(pair.second, true);
                        estimate = std::min(estimate, bidistance[0][pair.second] + bidistance[1][pair.second]);
                    }
                }
            }
            // reverse graph search
            if (!diqueue[1].empty()) {
                pair = diqueue[1].pop();
                
                if (bidistance[1][pair.second] == pair.first) {
                    if (pair.first <= estimate) {
                        update(pair.second, false);
                        estimate = std::min(estimate, bidistance[0][pair.second] + bidistance[1][pair.second]);
                    }
                }
            }
            if (diqueue[0].empty() && diqueue[1].empty()) break;
        }
        if (estimate >= INFINITY ) return -1;
        return estimate;
    }

private:
    // Try to relax the node v using distance d either in the forward or in the backward search
    void update(Vertex v, bool forward) {
        // Implement this method yourself
        if (forward){
            for (std::pair<Vertex, Distance> node : outgoing_edges[v]){
                if (bidistance[0][node.first] > bidistance[0][v] + node.second && rank[node.first] > rank[v]) {
                    if (bidistance[0][node.first] == INFINITY ) { workset.insert(node.first); }
                    bidistance[0][node.first] = bidistance[0][v] + node.second;
                    diqueue[0].update(node.first, bidistance[0][node.first]);
                }
            } 
        }
        else {
            for (std::pair<Vertex, Distance> node : incoming_edges[v]){
                if (bidistance[1][node.first] > bidistance[1][v] + node.second && rank[node.first] > rank[v]) { 
                    if (bidistance[1][node.first] == INFINITY ) { workset.insert(node.first); }
                    bidistance[1][node.first] = bidistance[1][v] + node.second;
                    diqueue[1].update(node.first, bidistance[1][node.first]);
                }
            } 
        }
    }

    void set_n(Vertex n) {
        N = n;
        outgoing_edges.resize(n);
        incoming_edges.resize(n);
    }


    void add_edge_to_list(std::vector<std::pair<Vertex,Distance>>& list, Vertex w, Distance c) {
        for (long long i = 0; i < list.size(); ++i) {
            std::pair<Vertex, Distance>& p = list[i];
            if (p.first == w) {
                if (p.second > c) {
                    p.second = c;
                }
                return;
            }
        }
        list.push_back(std::make_pair(w, c));
    }

    void add_directed_edge(Vertex u, Vertex v, Distance c) {
        add_edge_to_list(outgoing_edges[u], v, c);
        add_edge_to_list(incoming_edges[v], u, c);
    }

    void add_edge(Vertex u, Vertex v, Distance c) {
        add_directed_edge(u, v, c);
    }

    bool read_stdin() {
        int u,v,n,m;
        long long c;
        assert(scanf("%d %d", &n, &m) == 2);
        set_n(n);
        for (long long i = 0; i < m; ++i) {
            assert(scanf("%d %d %lld", &u, &v, &c) == 3);
            add_edge(u-1, v-1, c);
        }
        return true;
    }
};

int main() {
    Graph g;
    g.preprocess();
    std::cout << "Ready" << std::endl;

    int t;
    assert(scanf("%d", &t) == 1);
    int u, v;
    for (long long i = 0; i < t; ++i) {
        assert(scanf("%d %d", &u, &v) == 2);
        printf("%lld\n", g.query(u-1, v-1));
    }
}
