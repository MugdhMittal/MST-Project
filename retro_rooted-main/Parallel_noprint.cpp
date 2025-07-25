#include <bits/stdc++.h>
using namespace std;
const double INF = 10000000;
const int MAX_EDGE_ID = 1000000;
#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif // Number of threads for parallel processing

struct Edge
{
    int x, y, id; // Nodes x and y connected by this edge, id is the index of the edge
    double w;     // Weight of the edge
    // replacedId indicates which changed edge (by index) might replace this edge
    // -1 meaning no replacement
    int replacedId = -1;
    Edge() : x(-1), y(-1), w(0.0), id(-1), replacedId(-1) {} // constructor
    Edge(int u, int v, double w, int id_) : x(u), y(v), w(w), id(id_), replacedId(-1) {}
    bool operator==(const Edge &other) const
    {
        return (x == other.x && y == other.y && w == other.w && replacedId == other.replacedId);
    }

    bool operator<(const Edge &other) const
    {
        return w < other.w;
    }
};

// Structure to store Vertex (Information for Rooted Tree)
struct Vertex
{
    int parent;
    int root;
    // maxE: The maximum weighted edge on the path from this vertex to the root
    Edge maxE;
    Vertex()
    {
        parent = -1;
        root = -1;
        maxE = Edge();
    }
};

enum class EdgeStatus
{
    NONE = 0, // added to remainder
    INS = 1,  // inserted to key edges
    DEL = 2,  // deleted from key edges
    RPL = 3   // replaced an existing key edge
};

enum class OperationType
{
    INSERT = 0,
    DELETE = 1
};

// Disjoint Set Union-Find with path compression and union by rank
class DSU
{
public:
    vector<int> parent, rank;

    DSU(int n)
    {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; ++i)
            parent[i] = i;
    }

    int find(int x)
    {
        if (parent[x] != x)
            parent[x] = find(parent[x]);
        return parent[x];
    }

    bool union_sets(int x, int y)
    {
        int rootX = find(x), rootY = find(y);
        if (rootX == rootY)
            return false;

        if (rank[rootX] > rank[rootY])
            parent[rootY] = rootX;
        else if (rank[rootX] < rank[rootY])
            parent[rootX] = rootY;
        else
        {
            parent[rootY] = rootX;
            rank[rootX]++;
        }
        return true;
    }
};

// Function to find the MST and divide edges into Ex and Er
pair<vector<Edge>, vector<Edge>> findMST(vector<Edge> &edges, int n)
{
    vector<Edge> filtered;
    for (Edge &e : edges)
    {
        if (e.w != INF)
            filtered.push_back(e);
    }

    sort(filtered.begin(), filtered.end());

    DSU dsu(n);
    vector<Edge> Ex, Er;
    int edgeCounter = 0;

    for (Edge &edge : filtered)
    {
        if (edge.id == -1)
            edge.id = edgeCounter++;

        if (dsu.union_sets(edge.x, edge.y))
            Ex.push_back(edge);
        else
            Er.push_back(edge);
    }

    return {Ex, Er};
}

int findMinHeightVertex(vector<vector<pair<int, int>>> &adjEx, vector<bool> &visited, int start)
{
    int V = adjEx.size();
    vector<int> degree(V, 0);
    queue<int> q;

    vector<bool> localVisited(V, false);
    queue<int> bfs;
    bfs.push(start);
    localVisited[start] = true;

    // First, gather component starting from 'start'
    while (!bfs.empty())
    {
        int u = bfs.front();
        bfs.pop();
        for (auto &[v, _] : adjEx[u])
        {
            if (!localVisited[v])
            {
                localVisited[v] = true;
                bfs.push(v);
            }
        }
    }

    // Mark global visited[]
    for (int i = 0; i < V; ++i)
        if (localVisited[i])
            visited[i] = true;

    // Build degree[] and identify leaves
    for (int i = 0; i < V; ++i)
    {
        if (localVisited[i])
        {
            degree[i] = adjEx[i].size();
            if (degree[i] <= 1)
                q.push(i);
        }
    }

    int remaining = std::count(localVisited.begin(), localVisited.end(), true);
    while (remaining > 2 && !q.empty())
    {
        int sz = q.size();
        remaining -= sz;
        for (int i = 0; i < sz; ++i)
        {
            int leaf = q.front();
            q.pop();
            for (auto &[nbr, _] : adjEx[leaf])
            {
                if (--degree[nbr] == 1)
                    q.push(nbr);
            }
        }
    }

    // Return one of the remaining centroids
    while (!q.empty())
    {
        int r = q.front();
        q.pop();
        if (localVisited[r])
            return r;
    }

    cerr << "[ERROR] Could not find root for component starting at " << start << "\n";
    return -1;
}

vector<vector<pair<int, int>>> buildExAdj(vector<Edge> &Ex, int V)
{

    vector<vector<pair<int, int>>> adjEx(V); // (neighbor, index_of_edge_in_Ex)
    for (int i = 0; i < Ex.size(); i++)
    {
        Edge &e = Ex[i];
        if (e.x < 0 || e.y < 0)
            continue;
        adjEx[e.x].push_back({e.y, i});
        adjEx[e.y].push_back({e.x, i});
    }

    return adjEx;
}

class RootedTree
{
public:
    int V;
    vector<Vertex> RootedT;
    vector<Edge> Ex;
    vector<Edge> Er;
    double MSTWeight = 0.0;             // Total weight of the MST
    std::unordered_set<int> deletedIds; // Track deleted edge ids
    RootedTree() : V(0) {}
    RootedTree(int V)
    {
        this->V = V;
        RootedT.resize(V);
    }
    RootedTree(const RootedTree &other)
    {
        V = other.V;
        RootedT = other.RootedT;
        Ex = other.Ex;
        Er = other.Er;
        deletedIds = other.deletedIds;
        MSTWeight = other.MSTWeight;
    }

    // Copy assignment operator
    RootedTree &operator=(const RootedTree &other)
    {
        if (this != &other)
        {
            V = other.V;
            RootedT = other.RootedT;
            Ex = other.Ex;
            Er = other.Er;
            deletedIds = other.deletedIds;
            MSTWeight = other.MSTWeight;
        }
        return *this;
    }

    void reverseParentPath(int from, int to)
    {
        vector<int> path;
        int x = from;
        while (x != to)
        {
            path.push_back(x);
            x = RootedT[x].parent;
        }
        path.push_back(to);

        for (int i = 0; i + 1 < path.size(); i++)
        {
            RootedT[path[i]].parent = path[i + 1];
        }

        RootedT[path.back()].parent = to; // keep root pointing to itself
    }

    // We Create the Rooted Tree (Algorithm 3)
    //  We find a vertex with Minimum Height vertex as root and run BFS from it.
    //  If multiple components exist, we repeat.

    // BFS begins. Use parallelisaztion for BFS
    void bfsFromRoot(int root, vector<vector<pair<int, int>>> &adjEx)
    {
        queue<int> q;
        // Push all neighbors of root
        for (auto &nbr : adjEx[root])
        {
            int v = nbr.first;
            int eidx = nbr.second;
            if (RootedT[v].parent == -1)
            {
                RootedT[v].parent = root;
                RootedT[v].root = root;
                RootedT[v].maxE = Ex[eidx]; // edge (root->v)
                q.push(v);
            }
        }

        while (!q.empty())
        {
            int cur = q.front();
            q.pop();
            for (auto &nbr : adjEx[cur])
            {
                int nxt = nbr.first;
                int eidx = nbr.second;
                if (RootedT[nxt].parent == -1)
                {
                    RootedT[nxt].parent = cur;
                    RootedT[nxt].root = RootedT[cur].root;
                    // Check if this edge is max on path
                    if (Ex[eidx].w > RootedT[cur].maxE.w)
                    {
                        RootedT[nxt].maxE = Ex[eidx];
                    }
                    else
                    {
                        RootedT[nxt].maxE = RootedT[cur].maxE;
                    }
                    q.push(nxt);
                }
            }
        }
    }
    /// BFS ends.

    void Create_Tree(vector<vector<pair<int, int>>> &adjEx)
    {
        V = adjEx.size();
        RootedT.assign(V, Vertex());

        vector<bool> visited(V, false);

        for (int i = 0; i < V; i++)
        {
            if (!visited[i] && RootedT[i].parent == -1)
            {
                int r = findMinHeightVertex(adjEx, visited, i);
                if (r == -1)
                {
                    cerr << "[ERROR] Could not find a valid root. Component starting at " << i << " might be isolated.\n";
                    continue;
                }

                RootedT[r].parent = r;
                RootedT[r].root = r;
               //cout << "[DEBUG] Rooted tree built with root = " << r << "\n";
                bfsFromRoot(r, adjEx);
            }
        }
    }

    int findRoot(int v)
    {
        return RootedT[v].root;
    }

    // Utility function to get path to root
    std::vector<int> getPathToRoot(int x)
    {
        std::vector<int> path;
        while (x != this->RootedT[x].parent)
        {
            path.push_back(x);
            x = RootedT[x].parent;
        }
        path.push_back(x);
        return path;
    }

    // Function to find the maximum edge on the path from u to v
    Edge findMaxEdgeOnPath(int u, int v)
    {
        // If in different components
        if (RootedT[u].root != RootedT[v].root)
        {
            // Path passes through roots or is disconnected
            // If they are different roots => no direct path in the tree
            //  a dummy edge with w = -1
            return Edge();
        }

        int root = RootedT[u].root;
        // If the maximum edges on u->root and v->root differ, O(1)
        // else O(h) by traversing up.

        Edge umx = RootedT[u].maxE;
        Edge vmx = RootedT[v].maxE;

        if (umx.id == vmx.id)
        {
            std::vector<int> upath = getPathToRoot(u);
            std::vector<int> vpath = getPathToRoot(v);
            std::reverse(upath.begin(), upath.end());
            std::reverse(vpath.begin(), vpath.end());
            vector<vector<pair<int, int>>> adjEx = buildExAdj(Ex, V);

            int len = (int)std::min(upath.size(), vpath.size());
            int lca = -1;
            for (int i = 0; i < len; i++)
            {
                if (upath[i] == vpath[i])
                    lca = upath[i];
                else
                    break;
            }

            auto maxOnPath = [&](int start, int end)
            {
                double mxw = -1;
                Edge mxE(-1, -1, -1, -1);
                int cur = start;
                while (cur != end)
                {
                    int p = RootedT[cur].parent;
                    if (p == cur)
                        break;
                    // Find edge cur-p
                    bool found = false;
                    for (auto &nb : adjEx[cur])
                    {
                        if (nb.first == p)
                        {
                            Edge &ed = Ex[nb.second];
                            if (ed.w > mxw)
                            {
                                mxw = ed.w;
                                mxE = ed;
                            }
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        break;
                    cur = p;
                }
                return mxE.w >= 0 ? mxE : Edge();
            };
            Edge mxE1 = maxOnPath(u, lca);
            Edge mxE2 = maxOnPath(v, lca);

            if (mxE1.w > mxE2.w)
                return mxE1;
            return mxE2;
        }
        return umx.w > vmx.w ? umx : vmx;
    }

    void Classify_Edges(vector<Edge> &CE,
                        vector<EdgeStatus> &Status,
                        vector<Edge> &Marked,
                        vector<bool> &operation)
    {

        struct Replacement
        {
            std::atomic<int> ce_index;
            std::atomic<double> weight;
            Replacement() : ce_index(-1), weight(INF) {}
        };
        int n = CE.size();
        std::vector<Replacement> bestReplacements(MAX_EDGE_ID);

        Status.resize(n);
        Marked.resize(n);

        auto classifyWorker = [&](int start, int end)
        {
            for (int i = start; i < end; ++i)
            {
                Status[i] = EdgeStatus::NONE;
                Marked[i] = Edge();
                Edge &E = CE[i];
                bool op = operation[i];

                if (!op)
                {
                    Status[i] = EdgeStatus::DEL;
                    continue;
                }

                int ru = RootedT[E.x].root;
                int rv = RootedT[E.y].root;
                if (ru != rv)
                {
                    Status[i] = EdgeStatus::INS;
                    continue;
                }

                Edge umx = RootedT[E.x].maxE;
                Edge vmx = RootedT[E.y].maxE;

                Edge MaxW;
                if (umx.id != vmx.id)
                {
                    MaxW = (umx.w > vmx.w) ? umx : vmx;
                }
                else
                {
                    MaxW = findMaxEdgeOnPath(E.x, E.y); // fallback for Case 2b
                }
                if (MaxW.id != -1 && MaxW.w > E.w)
                {
                    int maxWid = MaxW.id;
                    if (maxWid >= MAX_EDGE_ID)
                        continue; // Avoid overflow

                    double oldW = bestReplacements[maxWid].weight.load();
                    int oldI = bestReplacements[maxWid].ce_index.load();

                    // Atomic write with tie-break
                    if (E.w < oldW || (E.w == oldW && i < oldI))
                    {
                        bestReplacements[maxWid].weight.store(E.w);
                        bestReplacements[maxWid].ce_index.store(i);
                        Marked[i] = MaxW;
                    }
                }
            }
        };

        std::vector<std::thread> threads;
        int chunkSize = (n + NUM_THREADS - 1) / NUM_THREADS;

        for (int t = 0; t < NUM_THREADS; ++t)
        {
            int start = t * chunkSize;
            int end = std::min(start + chunkSize, n);
            threads.emplace_back(classifyWorker, start, end);
        }

        for (auto &th : threads)
            th.join();

        // Finalize best replacements
        for (int maxWid = 0; maxWid < MAX_EDGE_ID; ++maxWid)
        {
            int i = bestReplacements[maxWid].ce_index.load();
            if (i == -1)
                continue;

            Status[i] = EdgeStatus::RPL;

            auto it = std::find_if(Ex.begin(), Ex.end(),
                                   [&](const Edge &e)
                                   { return e.id == Marked[i].id; });

            if (it != Ex.end())
            {
                Marked[i].replacedId = std::distance(Ex.begin(), it);
            }
            else
            {
                Marked[i].replacedId = -1;
                cerr << "[ERROR] Could not locate edge id=" << Marked[i].id << " in Ex!" << endl;
            }

            //cout << "[DEBUG] Marking edge for replacement: "<< "MaxW.id = " << Marked[i].id<< ", rpl_idx = " << Marked[i].replacedId<< ", E.id = " << CE[i].id<< ", MaxW.w = " << Marked[i].w<< ", E.w = " << CE[i].w << endl;
        }
    }

    void Process_Status(vector<Edge> &CE, vector<EdgeStatus> &Status, vector<Edge> &Marked)
    {
        //cout << "[DEBUG] → Process_Status: parallel, paper-aligned version\n";

        if (Status.size() != CE.size() || Marked.size() != CE.size())
        {
            cerr << "[FATAL] Size mismatch in Process_Status" << endl;
            exit(1);
        }

        // Shared post-merge containers
        std::unordered_set<int> deletedNow;
        std::unordered_set<int> replacedInEx;
        std::vector<Edge> insertToEr;
        std::vector<Edge> insertToEx;

        // Thread-local containers
        std::vector<std::unordered_set<int>> thread_deletedNow(NUM_THREADS);
        std::vector<std::unordered_set<int>> thread_replacedInEx(NUM_THREADS);
        std::vector<std::vector<Edge>> thread_insertToEr(NUM_THREADS);
        std::vector<std::vector<Edge>> thread_insertToEx(NUM_THREADS);
        std::vector<std::vector<std::pair<int, Edge>>> thread_exReplacements(NUM_THREADS); // <index, edge>

        int chunkSize = (CE.size() + NUM_THREADS - 1) / NUM_THREADS;
        std::vector<std::thread> threads;

        for (int t = 0; t < NUM_THREADS; ++t)
        {
            threads.emplace_back([&, t]()
                                 {
            int start = t * chunkSize;
            int end = std::min((int)CE.size(), start + chunkSize);

            for (int i = start; i < end; ++i)
            {
                Edge &E = CE[i];
                EdgeStatus S = Status[i];

                if (deletedIds.count(E.id))
                    continue;

                switch (S)
                {
                case EdgeStatus::DEL:
                    thread_deletedNow[t].insert(E.id);
                    break;

                case EdgeStatus::INS:
                    thread_insertToEx[t].push_back(E);
                    break;

                case EdgeStatus::NONE:
                    thread_insertToEr[t].push_back(E);
                    break;

                case EdgeStatus::RPL:
                {
                    Edge &Erpl = Marked[i];

                    if (deletedIds.count(Erpl.id))
                        break;

                    if (Erpl.replacedId >= 0 &&
                        Erpl.replacedId < (int)Ex.size() &&
                        Ex[Erpl.replacedId].id == Erpl.id)
                    {
                        thread_exReplacements[t].emplace_back(Erpl.replacedId, E);
                        thread_replacedInEx[t].insert(Erpl.id);
                    }
                    else
                    {
                        thread_insertToEr[t].push_back(E);
                        thread_insertToEr[t].push_back(Erpl);
                        cerr << "[WARN] RPL fallback to Er: E.id=" << E.id
                             << ", Erpl.id=" << Erpl.id
                             << ", replacedId=" << Erpl.replacedId << "\n";
                    }

                    thread_insertToEr[t].push_back(Erpl); // keep for safety
                    break;
                }

                default:
                    break;
                }
            }

            //cout << "[DEBUG] Thread " << t << " done.\n"; 
            });
        }

        for (auto &th : threads)
            th.join();

        // Safe merging into shared state
        for (int t = 0; t < NUM_THREADS; ++t)
        {
            deletedNow.insert(thread_deletedNow[t].begin(), thread_deletedNow[t].end());
            for (int id : thread_deletedNow[t])
                deletedIds.insert(id);

            replacedInEx.insert(thread_replacedInEx[t].begin(), thread_replacedInEx[t].end());
            insertToEr.insert(insertToEr.end(), thread_insertToEr[t].begin(), thread_insertToEr[t].end());
            insertToEx.insert(insertToEx.end(), thread_insertToEx[t].begin(), thread_insertToEx[t].end());

            // Apply replacements safely to Ex
            for (auto &[idx, newEdge] : thread_exReplacements[t])
            {
                if (idx >= 0 && idx < (int)Ex.size())
                    Ex[idx] = newEdge;
            }
        }

        bool rebuildRequired = !deletedNow.empty() || !replacedInEx.empty();

        if (rebuildRequired)
        {
            std::vector<Edge> validEdges;

            for (const Edge &e : Ex)
            {
                if (!deletedIds.count(e.id) &&
                    !replacedInEx.count(e.id) &&
                    e.w != INF && e.w != -1)
                {
                    validEdges.push_back(e);
                }
            }

            for (const Edge &e : Er)
            {
                if (!deletedIds.count(e.id) && e.w != INF && e.w != -1)
                    validEdges.push_back(e);
            }

            for (const Edge &e : insertToEx)
                if (e.w != INF && e.w != -1)
                    validEdges.push_back(e);

            for (const Edge &e : insertToEr)
                if (e.w != INF && e.w != -1)
                    validEdges.push_back(e);

            tie(Ex, Er) = findMST(validEdges, V);

            //cout << "[DEBUG] Rebuilding MST: Ex.size = " << Ex.size() << ", Er.size = " << Er.size() << "\n";
        }
        else
        {
            for (const Edge &e : insertToEx)
                if (!deletedIds.count(e.id) && e.w != INF && e.w != -1)
                    Ex.push_back(e);

            for (const Edge &e : insertToEr)
                if (!deletedIds.count(e.id) && e.w != INF && e.w != -1)
                    Er.push_back(e);
        }

        Er.erase(remove_if(Er.begin(), Er.end(),
                           [&](const Edge &e)
                           {
                               return deletedIds.count(e.id) || e.w == INF || e.w == -1;
                           }),
                 Er.end());

        auto adj = buildExAdj(Ex, V);
        Create_Tree(adj);

        MSTWeight = 0.0;
        for (const Edge &e : Ex)
            if (e.w != INF && e.w != -1)
                MSTWeight += e.w;

        //cout << "[DEBUG] Process_Status done: Ex.size = " << Ex.size()<< ", Er.size = " << Er.size()<< ", MSTWeight = " << MSTWeight << "\n";
    }

    void Repair_Tree(vector<Edge> &Er, vector<EdgeStatus> &Status, vector<Edge> &Marked)
    {
        //cout << "[DEBUG] Repair_Tree called\n";

        bool changed;
        int loopCount = 0;

        struct Replacement
        {
            std::atomic<int> er_index;
            std::atomic<double> weight;
            Replacement() : er_index(-1), weight(INF) {}
        };

        do
        {
            changed = false;
            loopCount++;

            if (loopCount > 20)
            {
                cerr << "[FATAL] Repair_Tree exceeded safe loop limit. Exiting.\n";
                break;
            }

            Status.assign(Er.size(), EdgeStatus::NONE);
            Marked.assign(Er.size(), Edge());
            std::vector<Replacement> bestReplacements(MAX_EDGE_ID);

            // [1] Threaded worker for classification
            auto repairWorker = [&](int start, int end)
            {
                for (int i = start; i < end; ++i)
                {
                    Edge &E = Er[i];

                    if (deletedIds.count(E.id) || E.w == INF)
                        continue;

                    int ru = RootedT[E.x].root;
                    int rv = RootedT[E.y].root;

                    if (ru != rv)
                    {
                        Status[i] = EdgeStatus::INS;
                        changed = true;
                        continue;
                    }

                    Edge umx = RootedT[E.x].maxE;
                    Edge vmx = RootedT[E.y].maxE;
                    Edge MaxW;

                    if (umx.id != vmx.id)
                    {
                        MaxW = (umx.w > vmx.w) ? umx : vmx;
                    }
                    else
                    {
                        MaxW = findMaxEdgeOnPath(E.x, E.y); // fallback
                    }

                    if (MaxW.id != -1 && MaxW.w > E.w)
                    {
                        int maxWid = MaxW.id;
                        if (maxWid >= MAX_EDGE_ID)
                            continue;

                        double oldW = bestReplacements[maxWid].weight.load();
                        int oldI = bestReplacements[maxWid].er_index.load();

                        if (E.w < oldW || (E.w == oldW && i < oldI))
                        {
                            bestReplacements[maxWid].weight.store(E.w);
                            bestReplacements[maxWid].er_index.store(i);
                            Marked[i] = MaxW;
                            changed = true;
                        }
                    }
                }
            };

            // [2] Launch threads
            std::vector<std::thread> threads;
            int chunkSize = (Er.size() + NUM_THREADS - 1) / NUM_THREADS;

            for (int t = 0; t < NUM_THREADS; ++t)
            {
                int start = t * chunkSize;
                int end = std::min(start + chunkSize, (int)Er.size());
                threads.emplace_back(repairWorker, start, end);
            }
            for (auto &th : threads)
                th.join();

            // [3] Finalize replacements
            for (int maxWid = 0; maxWid < MAX_EDGE_ID; ++maxWid)
            {
                int i = bestReplacements[maxWid].er_index.load();
                if (i == -1)
                    continue;

                Status[i] = EdgeStatus::RPL;

                auto it = std::find_if(Ex.begin(), Ex.end(),
                                       [&](const Edge &e)
                                       { return e.id == Marked[i].id; });

                if (it != Ex.end())
                {
                    Marked[i].replacedId = std::distance(Ex.begin(), it);
                }
                else
                {
                    Marked[i].replacedId = -1;
                    cerr << "[ERROR] Could not locate edge id=" << Marked[i].id << " in Ex!" << endl;
                }
            }

            // [4] Apply Status to update Ex and Er (no tree repair here)
            for (int i = 0; i < (int)Er.size(); i++)
            {
                Edge &E = Er[i];
                if (Status[i] == EdgeStatus::INS)
                {
                    Ex.push_back(E);
                    //cout << "[DEBUG] Repair inserted edge (" << E.x << "," << E.y << ") id=" << E.id << endl;
                }
                else if (Status[i] == EdgeStatus::RPL)
                {
                    Edge &Erpl = Marked[i];
                    if (Erpl.replacedId >= 0 && Erpl.replacedId < (int)Ex.size())
                    {
                        Ex.erase(Ex.begin() + Erpl.replacedId);
                        Ex.push_back(E);

                        //cout << "[DEBUG] Replaced edge " << Erpl.id<< " with " << E.id << " → rebuilding rooted tree\n";

                        auto adjEx = buildExAdj(Ex, V);
                        Create_Tree(adjEx); // immediately fix the tree
                    }
                    else
                    {
                        cerr << "[ERROR] Invalid replacedId in repair: " << Erpl.replacedId << endl;
                    }

                    Er.push_back(Erpl); // Retain replaced edge in Er
                }
            }

            // Clean Er
            Er.erase(remove_if(Er.begin(), Er.end(),
                               [&](const Edge &e)
                               {
                                   return deletedIds.count(e.id) || e.w == INF ||
                                          any_of(Ex.begin(), Ex.end(),
                                                 [&](const Edge &ex)
                                                 { return ex.id == e.id; });
                               }),
                     Er.end());

            // Deferred tree rebuild
            if (changed)
            {
                auto adjEx = buildExAdj(Ex, V);
                Create_Tree(adjEx);
                //cout << "[DEBUG] Tree rebuilt after batched replacements\n";
            }

        } while (changed);

        //cout << "[DEBUG] Repair_Tree completed\n";
    }

    void ProcessAllEdges(vector<Edge> &CE, vector<bool> &operation)
    {
        //std::cout << "[DEBUG] ProcessAllEdges called\n";

        if (CE.empty())
        {
            //std::cout << "[DEBUG] No edges to process\n";
            return;
        }

        vector<EdgeStatus> Status;
        vector<Edge> Marked;

        int loopCounter = 0;
        while (!CE.empty())
        {
            if (++loopCounter > 20)
            {
                std::cerr << "[FATAL] ProcessAllEdges exceeded safe loop limit. Exiting loop.\n";
                break;
            }

            //std::cout << "[DEBUG] Processing CE loop. CE.size = " << CE.size() << "\n";

            // Step 1: Assign default statuses
            Status.assign(CE.size(), EdgeStatus::NONE);
            Marked.assign(CE.size(), Edge());

            // Step 2: Classify edges
            Classify_Edges(CE, Status, Marked, operation);

            // Step 3: Filter out deleted edges before processing
            std::vector<Edge> newCE;
            std::vector<bool> newOp;
            std::vector<EdgeStatus> newStatus;
            std::vector<Edge> newMarked;

            for (int i = 0; i < (int)CE.size(); ++i)
            {
                if (!deletedIds.count(CE[i].id))
                {
                    newCE.push_back(CE[i]);
                    newOp.push_back(operation[i]);
                    newStatus.push_back(Status[i]);
                    newMarked.push_back(Marked[i]);
                }
            }

            CE = std::move(newCE);
            operation = std::move(newOp);
            Status = std::move(newStatus);
            Marked = std::move(newMarked);

            operation.resize(CE.size());
            Status.resize(CE.size());
            Marked.resize(CE.size());

            // Step 4: Process status changes
            Process_Status(CE, Status, Marked);

            // Step 5: Do NOT re-add already processed CE entries
            // Clear CE and operation so that loop exits once all edges are handled
            CE.clear();
            operation.clear();
        }

        // Final sanity check
        if (Ex.size() > V - 1)
        {
            std::cerr << "[ERROR] MST has too many edges: Ex.size = " << Ex.size()
                      << ", expected = " << (V - 1) << std::endl;
        }

        // Final repair to ensure full MST
        if (Ex.size() < V - 1)
        {
            //std::cout << "[DEBUG] Final MST repair triggered from ProcessAllEdges()\n";
            vector<EdgeStatus> dummyStatus;
            vector<Edge> dummyMarked;
            Repair_Tree(Er, dummyStatus, dummyMarked);
        }

        //std::cout << "[DEBUG] ProcessAllEdges completed. Ex.size = " << Ex.size() << "\n";
    }
};

const int T = 50011;         // Timeline size
const int Tc = sqrt(T) + 10; // Timeline size for block decomposition (+10 for optimisation and incase N is small +10 reduces complications)

// Global variables
int tp; // Key Blocks of Timeline
int n;  // Number of nodes
int m;  // Number of edges
int q;  // Number of queries
vector<vector<pair<Edge, OperationType>>> que;

vector<Edge> getMST(int t, vector<RootedTree> &rt)
{
    //cout << "[DEBUG] getMST called with t = " << t << endl;

    if (t < 0)
    {
        //cout << "[DEBUG] Time t is negative. Returning empty MST." << endl;
        return vector<Edge>();
    }

    if (t >= tp * Tc)
    {
        //cout << "[DEBUG] Time t exceeds maximum range (tp * M). Returning empty MST." << endl;
        return vector<Edge>();
    }

    int bkt = t / Tc;
    //cout << "[DEBUG] Processing block #" << bkt << " (bkt = " << bkt << ")" << endl;

    if (bkt >= rt.size())
    {
        //cout << "[ERROR] Block index out of range for rt. Returning empty MST." << endl;
        return vector<Edge>();
    }

    vector<Edge> CE;
    vector<bool> operation;

    for (int i = bkt * Tc; i <= t && i < que.size(); ++i)
    {
        if (!que[i].empty())
        {
            //cout << "[DEBUG] Processing que[" << i << "] with " << que[i].size() << " updates." << endl;
            for (auto &entry : que[i])
            {
                if (i > t)
                    continue;
                Edge edge = entry.first;
                OperationType op = entry.second;

                if (op == OperationType::DELETE && rt[bkt].deletedIds.count(edge.id))
                {
                    //cout << "[DEBUG] Skipping already deleted edge id = " << edge.id << " from que[" << i << "]" << endl;
                    continue;
                }

                if (op == OperationType::INSERT)
                {
                    auto &Ex = rt[bkt].Ex;
                    if (std::any_of(Ex.begin(), Ex.end(), [&](const Edge &e)
                                    { return e.id == edge.id; }))
                    {
                       // cout << "[DEBUG] Skipping already inserted edge id = " << edge.id << " from que[" << i << "]" << endl;
                        continue;
                    }

                    if (rt[bkt].deletedIds.count(edge.id))
                    {
                        //cout << "[DEBUG] Skipping reinsertion of deleted edge id = " << edge.id << " from que[" << i << "]" << endl;
                        continue;
                    }
                }

                CE.push_back(edge);
                operation.push_back(op == OperationType::INSERT);
                //cout << "  [DEBUG] Scheduled " << (op == OperationType::INSERT ? "INSERT" : "DELETE")<< " of edge (" << edge.x << ", " << edge.y<< ", w=" << edge.w << ", id=" << edge.id << ")" << endl;
            }
        }
    }

    //cout << "[DEBUG] Total edges to process: " << CE.size() << endl;

    RootedTree working = rt[bkt]; // ← Important: use a copy to avoid state contamination
    working.ProcessAllEdges(CE, operation);
    vector<Edge> mst = working.Ex;

    //cout << "[DEBUG] getMST returning MST with " << mst.size() << " edges." << endl;
    //cout << "[DEBUG] Total MST weight: " << working.MSTWeight << endl;
    double recomputed = 0;
    for (const auto &e : mst)
        if (e.w != INF && e.w != -1)
            recomputed += e.w;

    if (std::abs(recomputed - working.MSTWeight) > 1e-6)
    {
        cerr << "[WARNING] MST weight mismatch: expected " << recomputed
             << ", actual = " << working.MSTWeight << "\n";
    }
    return mst;
}

void addAndDeleteEdge(int x, int y, int w, int t, int id,
                      vector<RootedTree> &rt, OperationType op)
{
    //cout << "[DEBUG] addAndDeleteEdge called with x=" << x<< ", y=" << y << ", w=" << w << ", t=" << t<< ", id=" << id << ", op=" << (op == OperationType::INSERT ? "INSERT" : "DELETE") << endl;

    if (t < 0 || t >= T)
    {
        //cout << "[ERROR] Time t is out of bounds. Skipping." << endl;
        return;
    }

    // Ensure que is large enough
    if ((int)que.size() <= t)
        que.resize(T);

    // Propagate to all remaining times in this block
    int bkt = t / Tc;
    int blockEnd = std::min(T - 1, (bkt + 1) * Tc - 1);
    for (int i = t; i <= blockEnd; ++i)
    {
        que[i].push_back({Edge(x, y, w, id), op});
    }

    // Then propagate to future blocks
    for (int b = bkt + 1; b < tp; ++b)
    {
        int blockTime = b * Tc;
        if ((int)que.size() <= blockTime)
            que.resize(blockTime + 1);
        que[blockTime].push_back({Edge(x, y, w, id), op});
    }
}

int main()
{
    cout << "[DEBUG] Starting Parallel MST program\n";
    // Start measuring time
    clock_t tStart = clock();
    // Redirect cout to file (except final time output)

   // std::ofstream debugOut("debug_output.txt");
    //if (!debugOut.is_open())
    //{
     //   std::cerr << "Error opening debug output file.\n";
      //  return 1;
    //}
    //std::streambuf *coutbuf = std::cout.rdbuf(); // Save old buffer
    //std::cout.rdbuf(debugOut.rdbuf());           // Redirect cout to file

    const std::string filename = "in000.txt";
    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        cerr << "Error opening input file." << endl;
        return 1;
    }
    infile >> n >> m >> q;
    //cout << "[DEBUG] Read n = " << n << ", m = " << m << ", q = " << q << endl;

    vector<Edge> edges(m);
    tp = (T / Tc) + 2;
    vector<RootedTree> RT(tp, RootedTree(n));

    for (int i = 0; i < m; ++i)
    {
        int x, y, w;
        infile >> x >> y >> w;
        edges[i] = Edge(--x, --y, w, i);
    }

    auto [baseEx, baseEr] = findMST(edges, n);
    RT[0].Ex = baseEx;
    RT[0].MSTWeight = 0.0;
    for (auto &e : baseEx)
        RT[0].MSTWeight += e.w;
    RT[0].Er = baseEr;
    vector<vector<pair<int, int>>> adjEx_0 = buildExAdj(baseEx, n);
    RT[0].Create_Tree(adjEx_0);

    for (int b = 1; b < tp; ++b)
    {
        RT[b] = RootedTree(RT[b - 1]);
    }

    int queryId = n + m;

    for (int qi = 0; qi < q; ++qi)
    {
        int op;
        infile >> op;

        if (op == 0)
        {
            int t;
            infile >> t;
            //cout << "[DEBUG] Performing MST query at time " << t << endl;
            vector<Edge> mst = getMST(t, RT);
            if (mst.empty())
            {
              //  cout << "[WARNING] MST is empty at time t = " << t << ".\n";
                continue;
            }
           // cout << "MST at time " << t << ":\n";
        }
        else if (op == 1)
        {
            int x, y, w, t;
            infile >> x >> y >> w >> t;
            x--;
            y--;
            addAndDeleteEdge(x, y, w, t, queryId++, RT, OperationType::INSERT);
        }
        else if (op == 2)
        {
            int x, y, w, t;
            infile >> x >> y >> w >> t;
            x--;
            y--;

            int realId = queryId;
            bool found = false;

            for (int b = tp - 1; b >= 0 && !found; --b)
            {
                for (auto &e : RT[b].Ex)
                {
                    if (((e.x == x && e.y == y) || (e.x == y && e.y == x)) && std::abs(e.w - w) < 1e-6)
                    {
                        realId = e.id;
                        found = true;
                        break;
                    }
                }
                for (auto &e : RT[b].Er)
                {
                    if (((e.x == x && e.y == y) || (e.x == y && e.y == x)) && std::abs(e.w - w) < 1e-6)
                    {
                        realId = e.id;
                        found = true;
                        break;
                    }
                }
            }

            if (!found)
            {
                for (int b = t / Tc; b < tp && !found; ++b)
                {
                    int blockTime = b * Tc;
                    if (blockTime < (int)que.size())
                    {
                        for (auto &[e, opType] : que[blockTime])
                        {
                            if (opType == OperationType::INSERT &&
                                ((e.x == x && e.y == y) || (e.x == y && e.y == x)) &&
                                std::abs(e.w - w) < 1e-6)
                            {
                                realId = e.id;
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }

            if (!found)
            {
              // cout << "[WARNING] Could not find edge (" << x << "-" << y << ", w=" << w << ") for deletion.\n";
            }

            addAndDeleteEdge(x, y, w, t, realId, RT, OperationType::DELETE);
            queryId++;
        }
        else
        {
            //cout << "[ERROR] Invalid operation type in input: " << op << endl;
        }
    }

    infile.close();

   // std::cout.rdbuf(coutbuf); // Restore cout to console
    //debugOut.close();         // Close the debug output file
    //cout << "[DEBUG] Execution completed. Total queries processed: " << queryId - n - m << endl;
        printf("Execution Time: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    return 0;
}
