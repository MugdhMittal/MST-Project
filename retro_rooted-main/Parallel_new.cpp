#include <bits/stdc++.h>
#include <execution>
using namespace std;
const int INF = 1e8; // A large value to represent infinity
const int MAX_EDGE_ID = 10000000;

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif

struct Edge
{
    int x, y;                                                                          // End Points of Edge
    int id;                                                                            // Edge ID
    double w;                                                                          // Weight of Edge
    std::atomic<int> replacedID = -1;                                                  // ID of the edge that replaces this edge in the MST
    Edge() : x(0), y(0), w(0.0), id(-1), replacedID(-1) {}                             // Constructor
    Edge(int x, int y, double w, int id) : x(x), y(y), w(w), id(id), replacedID(-1) {} // Constructor with parameters

    // Custom copy constructor
    Edge(const Edge &other)
        : x(other.x), y(other.y), id(other.id), w(other.w), replacedID(other.replacedID.load()) {}

    // Custom copy assignment operator
    Edge &operator=(const Edge &other)
    {
        if (this != &other)
        {
            x = other.x;
            y = other.y;
            id = other.id;
            w = other.w;
            replacedID.store(other.replacedID.load());
        }
        return *this;
    }

    bool operator<(const Edge &e) const
    {
        return w < e.w; // Comparison based on weight
    }

    bool operator==(const Edge &e) const
    {
        return (x == e.x && y == e.y && w == e.w && replacedID == e.replacedID); // Comparison based on all attributes
    }
};

struct vertex
{
    int parent;                                      // Parent of the vertex in the MST
    int root;                                        // Root of the vertex in the MST
    Edge MaxE;                                       // Maximum edge in the path to the root
    vertex() : parent(-1), root(-1), MaxE(Edge()) {} // Constructor
};

enum class EdgeStatus
{
    NONE = 0 /*Added to Er*/,
    INS = 1 /*Added to Ex*/,
    DEL = 2 /*Removed from graph*/,
    RPL = 3 /*Replaced in MST*/
};

enum class OperationType
{
    INSERT = 0 /*Insert Edge*/,
    DELETE = 1 // Delete Edge
};

// Disjoint Set Union-Find structure with path compression and union by rank
class DSU
{
public:
    vector<int> parent, rank;

    DSU(int n)
    {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; ++i)
        {
            parent[i] = i; // Initialize each vertex as its own parent
        }
    };

    int find(int x)
    {
        if (parent[x] != x)
        {
            parent[x] = find(parent[x]); // Path compression
        }
        return parent[x];
    }

    bool unionSets(int x, int y)
    {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX != rootY)
        {
            if (rank[rootX] < rank[rootY])
            {
                parent[rootX] = rootY; // Union by rank
            }
            else if (rank[rootX] > rank[rootY])
            {
                parent[rootY] = rootX; // Union by rank
            }
            else
            {
                parent[rootY] = rootX;
                rank[rootX]++;
            }
            return true; // Return true if union was successful
        }
        else
        {
            return false; // Return false if they are already in the same set
        }
    }
};

// Function to find MST and divide it into Ex and Er
pair<vector<Edge>, vector<Edge>> findMST(vector<Edge> &edges, int n)
{
    sort(edges.begin(), edges.end()); // Sort edges by weight
    DSU dsu(n);                       // Initialize DSU
    vector<Edge> Ex;                  // Edges in the MST
    vector<Edge> Er;                  // Edges not in the MST

    for (const auto &edge : edges)
    {
        int rootX = dsu.find(edge.x);
        int rootY = dsu.find(edge.y);
        if (dsu.unionSets(edge.x, edge.y))
        {
            Ex.push_back(edge); // Add edge to MST
        }
        else
        {
            Er.push_back(edge); // Edge not in MST
        }
    }
    return {Ex, Er}; // Return the two sets of edges
}

// Build an Adjacency list for the MST edges
vector<vector<pair<int, int>>> buildExAdj(const vector<Edge> &Ex, int V)
{
    vector<vector<pair<int, int>>> adjEx(V); // Adjacency list for edges in the MST (neighbour, edge_id)
    for (const auto &edge : Ex)
    {
        adjEx[edge.x].emplace_back(edge.y, edge.id); // Add edge to adjacency list
        adjEx[edge.y].emplace_back(edge.x, edge.id); // Add reverse edge to adjacency list
    }
    return adjEx; // Return the adjacency list
}

// Function to find a vertex with minimum height in the tree- Root of the tree
int findMinHeightVertex(vector<vector<pair<int, int>>> &adjEx, vector<bool> &visited, int start)
{
    int V = adjEx.size();

    // Step 1: Extract connected component
    vector<int> component;
    queue<int> q;
    q.push(start);
    visited[start] = true;

    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        component.push_back(u);

        for (auto &[v, _] : adjEx[u])
        {
            if (!visited[v])
            {
                visited[v] = true;
                q.push(v);
            }
        }
    }

    // Step 2: Compute degree only for component nodes
    unordered_map<int, int> degree;
    for (int u : component)
        degree[u] = adjEx[u].size();
    // Step 3: Leaf pruning to find center(s)
    queue<int> leafQ;
    for (int u : component)
    {
        if (degree[u] <= 1)
            leafQ.push(u);
    }
    int lastPopped = -1;
    int remaining = component.size();
    while (remaining > 2)
    {
        int sz = leafQ.size();
        remaining -= sz;

        while (sz--)
        {
            int u = leafQ.front();
            leafQ.pop();
            degree[u] = 0; // Mark as removed

            for (auto &[v, _] : adjEx[u])
            {
                if (--degree[v] == 1)
                    leafQ.push(v);
            }
        }
    }
    // Step 4: Return one of the remaining centers
    while (!leafQ.empty())
    {
        lastPopped = leafQ.front();
        leafQ.pop();
    }

    return lastPopped; // Likely to be the best root
}

class RootedTree
{
public:
    int V;                                             // Number of vertices in the tree
    vector<vertex> RootedT;                            // Vector of vertices in the rooted tree
    vector<Edge> Ex;                                   // Edges in the MST
    vector<Edge> Er;                                   // Edges not in the MST
    double MSTWeight;                                  // Total weight of the MST
    RootedTree() : V(0), Ex(), Er(), MSTWeight(0.0) {} // Constructor

    RootedTree(int V)
    {
        this->V = V;       // Initialize number of vertices
        RootedT.resize(V); // Resize the vector of vertices
    }

    RootedTree(const RootedTree &other)
    {
        V = other.V;                 // Copy number of vertices
        RootedT = other.RootedT;     // Copy vector of vertices
        Ex = other.Ex;               // Copy edges in the MST
        Er = other.Er;               // Copy edges not in the MST
        MSTWeight = other.MSTWeight; // Copy total weight of the MST
    }
    RootedTree &operator=(const RootedTree &other)
    {
        if (this != &other)
        {
            V = other.V;                 // Copy number of vertices
            RootedT = other.RootedT;     // Copy vector of vertices
            Ex = other.Ex;               // Copy edges in the MST
            Er = other.Er;               // Copy edges not in the MST
            MSTWeight = other.MSTWeight; // Copy total weight of the MST
        }
        return *this; // Return the current object
    }

    // We Create the Rooted Tree (Algorithm 3)
    //  We find a vertex with Minimum Height vertex as root and run BFS from it.
    //  If multiple components exist, we repeat.
    void bfsbranch(int root, vector<vector<pair<int, int>>> &adjEx)
    {
        std::queue<int> frontier; // Queue for BFS
        frontier.push(root);      // Start BFS from the root
        while (!frontier.empty())
        {
            int u = frontier.front();        // Get the current vertex
            frontier.pop();                  // Remove it from the queue
            for (auto &[v, eidx] : adjEx[u]) // Iterate over adjacent vertices
            {
                if (RootedT[v].parent == -1) // If not visited
                {
                    RootedT[v].parent = u;             // Set parent
                    RootedT[v].root = RootedT[u].root; // Set root
                    if (Ex[eidx].w > RootedT[u].MaxE.w)
                        RootedT[v].MaxE = Ex[eidx]; // Update max edge
                    else
                        RootedT[v].MaxE = RootedT[u].MaxE; // Keep max edge
                    frontier.push(v);                      // Add to queue for further exploration
                }
            }
        }
    }

    void bfsfromRoot(int root, vector<vector<pair<int, int>>> &adjEx)
    {
        std::vector<std::thread> threads;   // Vector to hold threads
        std::vector<int> rootchildren;      // Vector to hold children of the root
        for (auto &[v, eidx] : adjEx[root]) // Iterate over adjacent vertices of the root
        {
            if (RootedT[v].parent == -1) // If not visited
            {
                RootedT[v].parent = root;   // Set parent
                RootedT[v].root = root;     // Set root
                RootedT[v].MaxE = Ex[eidx]; // Set max edge
                rootchildren.push_back(v);  // Add child to the list
            }
        }
        size_t total = rootchildren.size();                          // Total number of children
        size_t chunk_size = (total + NUM_THREADS - 1) / NUM_THREADS; // Size of each chunk
        for (size_t i = 0; i < total; i += chunk_size)               // Iterate over chunks
        {
            size_t end = std::min(i + chunk_size, total); // End index for the chunk
            threads.emplace_back([&, i, end]()
                                 {
                for (size_t j = i; j < end; ++j) // Iterate over children in the chunk
                {
                    int u = rootchildren[j]; // Get the child vertex
                    for (auto &[v, eidx] : adjEx[u]) // Iterate over adjacent vertices
                    {
                        if (RootedT[v].parent == -1) // If not visited
                        {
                            RootedT[v].parent = u; // Set parent
                            RootedT[v].root = RootedT[u].root; // Set root
                            if (Ex[eidx].w > RootedT[u].MaxE.w)
                                RootedT[v].MaxE = Ex[eidx]; // Update max edge
                            else
                                RootedT[v].MaxE = RootedT[u].MaxE; // Keep max edge
                        }
                    }
                } });
        }
        for (auto &t : threads) // Join all threads
            t.join();
    }

    // Create the Rooted Tree from the edges in Ex
    //  We find a vertex with Minimum Height vertex as root and  run BFS from it.
    void Create_Tree(vector<vector<pair<int, int>>> &adjEx)
    {
        V = adjEx.size();

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
                cout << "[DEBUG] Rooted tree built with root = " << r << "\n";
                bfsfromRoot(r, adjEx);
            }
        }
    }
    std::vector<int> getPathToRoot(int x)
    {
        std::vector<int> path;
        while (x != RootedT[x].parent)
        {
            path.push_back(x);
            x = RootedT[x].parent;
        }
        path.push_back(x);
        return path; // Return the path from x to its root
    }

    Edge FindMaxEdgeOnPath(int u, int v, const vector<vector<pair<int, int>>> &adjEx)
    {
        // If in different components
        if (RootedT[u].root != RootedT[v].root)
        {
            // Path passes through roots or is disconnected
            // If they are different roots => no direct path in the tree
            //  a dummy edge with w = -1
            return Edge(-1, -1, -1, -1);
        }

        int root = RootedT[u].root;
        // If the maximum edges on u->root and v->root differ, O(1)
        // else O(h) by traversing up.

        Edge umx = RootedT[u].MaxE;
        Edge vmx = RootedT[v].MaxE;

        if (umx.id == vmx.id)
        {
            std::vector<int> upath = getPathToRoot(u);
            std::vector<int> vpath = getPathToRoot(v);
            std::reverse(upath.begin(), upath.end());
            std::reverse(vpath.begin(), vpath.end());

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

    void classifyEdges(vector<Edge> &CE, vector<EdgeStatus> &Status, vector<Edge> &Marked, vector<bool> &operations, vector<vector<pair<int, int>>> &adjEx)
    {
        if (Status.size() != CE.size() || Marked.size() != CE.size() || operations.size() != CE.size())
        {
            cerr << "[FATAL] Size mismatch in classifyEdges\n";
            exit(1);
        }

        size_t total = CE.size();
        size_t chunk_size = (total + NUM_THREADS - 1) / NUM_THREADS;
        std::vector<std::thread> threads;

        for (int t = 0; t < NUM_THREADS; ++t)
        {
            size_t start = t * chunk_size;
            size_t end = std::min(start + chunk_size, total);

            threads.emplace_back([&, start, end]()
                                 {
                for (size_t i = start; i < end; ++i)
                {
                    const Edge &E = CE[i];
                    Status[i] = EdgeStatus::NONE;
                    Marked[i] = Edge(); // dummy edge

                    if (operations[i] == static_cast<int>(OperationType::DELETE))
                    {
                        Status[i] = EdgeStatus::DEL;
                        continue;
                    }

                    int u = E.x;
                    int v = E.y;

                    if (RootedT[u].root != RootedT[v].root)
                    {
                        Status[i] = EdgeStatus::INS;
                        continue;
                    }

                    Edge MaxW = FindMaxEdgeOnPath(u, v, adjEx);

                    if (MaxW.w > E.w && MaxW.id >= 0)
                    {
                        Edge &sharedMax = Ex[MaxW.id];

                        int prev = sharedMax.replacedID.load();
                        while (true)
                        {
                            if (prev == -1 || (prev >= 0 && CE[prev].w > E.w))
                            {
                                if (sharedMax.replacedID.compare_exchange_weak(prev, i))
                                {
                                    Marked[i] = MaxW;
                                    Status[i] = EdgeStatus::RPL;
                                    break;
                                }
                                // CAS failed, retry with new prev
                            }
                            else break;
                        }
                    }
                } });
        }

        for (auto &th : threads)
            th.join();
    }

    std::unordered_set<int> Process_Status(
        vector<Edge> &CE,
        vector<EdgeStatus> &Status,
        vector<Edge> &Marked,
        const vector<vector<pair<int, int>>> &adjEx)
    {
        if (Status.size() != CE.size() || Marked.size() != CE.size())
        {
            cerr << "[FATAL] Size mismatch in Process_Status\n";
            exit(1);
        }

        std::vector<Edge> localEr[NUM_THREADS];
        std::vector<std::pair<int, Edge>> localExInsertions[NUM_THREADS];
        std::vector<Edge> localNewChangedEdges[NUM_THREADS];
        std::unordered_set<int> updatedVertices[NUM_THREADS];

        size_t total = CE.size();
        size_t chunk_size = (total + NUM_THREADS - 1) / NUM_THREADS;
        std::vector<std::thread> threads;

        for (int t = 0; t < NUM_THREADS; ++t)
        {
            size_t start = t * chunk_size;
            size_t end = std::min(start + chunk_size, total);

            threads.emplace_back([&, start, end, t]()
                                 {
                for (size_t i = start; i < end; ++i)
                {
                    Edge &E = CE[i];
                    EdgeStatus S = Status[i];

                    if (S == EdgeStatus::DEL)
                    {
                        E.w = INF;
                    }

                    if (S == EdgeStatus::NONE)
                    {
                        localEr[t].push_back(E);
                    }
                    else if (S == EdgeStatus::INS)
                    {
                        localExInsertions[t].emplace_back(E.id, E);
                        updatedVertices[t].insert(E.x);
                        updatedVertices[t].insert(E.y);
                    }
                    else if (S == EdgeStatus::RPL)
                    {
                        Edge &Erpl = Marked[i];

                        if (Erpl.id < 0)
                        {
                            localNewChangedEdges[t].push_back(E);
                            continue;
                        }

                        auto it = std::find_if(Ex.begin(), Ex.end(), [&](const Edge &edge) {
                            return edge.id == Erpl.id;
                        });

                        if (it == Ex.end())
                        {
                            localNewChangedEdges[t].push_back(E);
                            continue;
                        }

                        Edge &edgeInMST = *it;

                        if (edgeInMST.replacedID.load() == i)
                        {
                            localExInsertions[t].emplace_back(E.id, E);
                            edgeInMST.w = -1;
                            updatedVertices[t].insert(E.x);
                            updatedVertices[t].insert(E.y);
                        }
                        else
                        {
                            localNewChangedEdges[t].push_back(E);
                        }
                    }
                } });
        }

        for (auto &th : threads)
            th.join();

        // Merge Er, Ex, and updatedVertices
        for (int t = 0; t < NUM_THREADS; ++t)
        {
            for (Edge &e : localEr[t])
                Er.push_back(e);
            for (auto &[id, e] : localExInsertions[t])
                Ex.push_back(e);
        }

        // Save CE ← new changed edges
        CE.clear();
        for (int t = 0; t < NUM_THREADS; ++t)
            CE.insert(CE.end(), localNewChangedEdges[t].begin(), localNewChangedEdges[t].end());

        std::unordered_set<int> allUpdated;
        for (int t = 0; t < NUM_THREADS; ++t)
            allUpdated.insert(updatedVertices[t].begin(), updatedVertices[t].end());

        return allUpdated;
    }

    void Repair_Subtrees(const std::unordered_set<int> &roots, const vector<vector<pair<int, int>>> &adjEx)
    {
        std::vector<int> updated(roots.begin(), roots.end());
        size_t chunk_size = (updated.size() + NUM_THREADS - 1) / NUM_THREADS;
        std::vector<std::thread> threads;

        for (int t = 0; t < NUM_THREADS; ++t)
        {
            size_t start = t * chunk_size;
            size_t end = std::min(start + chunk_size, updated.size());

            threads.emplace_back([&, start, end]()
                                 {
                for (size_t i = start; i < end; ++i)
                {
                    int root = updated[i];

                    std::queue<int> q;
                    q.push(root);

                    RootedT[root].parent = root;
                    RootedT[root].root = root;
                    RootedT[root].MaxE = Edge();

                    std::unordered_set<int> visited;
                    visited.insert(root);

                    while (!q.empty())
                    {
                        int u = q.front();
                        q.pop();

                        for (const auto &[v, eid] : adjEx[u])
                        {
                            if (visited.count(v))
                                continue;

                            auto it = std::find_if(Ex.begin(), Ex.end(), [&](const Edge &e) {
                                return e.id == eid && e.w != INF && e.w >= 0;
                            });

                            if (it == Ex.end())
                                continue;

                            visited.insert(v);
                            q.push(v);

                            RootedT[v].parent = u;
                            RootedT[v].root = root;
                            RootedT[v].MaxE = (it->w > RootedT[u].MaxE.w) ? *it : RootedT[u].MaxE;
                        }
                    }
                } });
        }

        for (auto &th : threads)
            th.join();
    }

    void Process_Edges(vector<Edge> &CE, vector<OperationType> &operations, vector<vector<pair<int, int>>> &adjEx)
    {
        if (CE.empty())
            return;

        vector<EdgeStatus> Status(CE.size(), EdgeStatus::NONE);
        vector<Edge> Marked(CE.size(), Edge());

        // Convert OperationType vector to vector<bool> for classifyEdges
        vector<bool> opBools;
        opBools.reserve(operations.size());
        for (auto &op : operations)
            opBools.push_back(op == OperationType::INSERT);

        // Repeat until no new changed edges are produced
        while (!CE.empty())
        {
            classifyEdges(CE, Status, Marked, opBools, adjEx);

            std::unordered_set<int> updatedVertices = Process_Status(CE, Status, Marked, adjEx);

            if (!updatedVertices.empty())
            {
                Repair_Subtrees(updatedVertices, adjEx);
            }
        }
    }
};

const int T = 50011;         // Timeline size
const int Tc = sqrt(T) + 10; // Timeline size for block decomposition (+10 for optimisation and incase N is small +10 reduces complications)
// Global variables
int tp;                                               // Key Blocks of Timeline
int n;                                                // Number of nodes
int m;                                                // Number of edges
int q;                                                // Number of queries
vector<vector<pair<Edge, OperationType>>> que(T + 1); // Query queue for edges and operations
vector<int> lastCheckpoint;

vector<Edge> getMST(int t, int Tc, int n, const vector<vector<pair<Edge, OperationType>>> &que, vector<RootedTree> &rt)
{
    int bkt = t / Tc;
    vector<Edge> CE;
    vector<OperationType> operations;

    for (int i = bkt * Tc; i <= t; ++i)
    {
        if (!que[i].empty())
        {
            for (const auto &entry : que[i])
            {
                Edge edge = entry.first;
                OperationType op = entry.second;
                CE.push_back(edge);
                operations.push_back(op);
            }
        }
    }
    RootedTree &working = rt[bkt];
    auto adjEx = buildExAdj(working.Ex, n);
    working.Process_Edges(CE, operations, adjEx);
    // Return the current MST edges after processing
    return working.Ex;
}

void addAndDeleteEdge(int x, int y, int w, int t, int id, vector<RootedTree> &rt, OperationType op)
{
    que[t].emplace_back(Edge(x, y, w, id), op);
    int bkt = t / Tc;
    for (int b = bkt + 1; b < tp; ++b)
    {
        int blockTime = b * Tc;
        if ((int)que.size() <= blockTime)
            que.resize(blockTime + 1);
        que[blockTime].emplace_back(Edge(x, y, w, id), op);
    }
}

int main()
{
    cout << "Starting Parallel MST computation..." << endl;
    clock_t start = clock();
    const string filename = "in001.txt"; // Input file name
    ifstream infile(filename);
    if (!infile.is_open())
    {
        cerr << "Error opening input file." << endl;
        return 1;
    }
    infile >> n >> m >> q;
    cout << "Read n = " << n << ", m = " << m << ", q = " << q << endl;
    vector<Edge> edges(m);
    tp = (T / Tc) + 2;                        // Total number of time blocks
    vector<RootedTree> RT(tp, RootedTree(n)); // Vector of Rooted Trees for each time block
    for (int i = 0; i < m; ++i)
    {
        int x, y, w;
        infile >> x >> y >> w;
        edges[i] = Edge(--x, --y, w, i); // Store edges with 0-based indexing
    }
    auto [baseEx, baseEr] = findMST(edges, n); // Find the initial MST
    RT[0].Ex = baseEx;                         // Set the edges in the first Rooted Tree
    RT[0].Er = baseEr;                         // Set the edges not in the MST
    RT[0].MSTWeight = 0.0;                     // Initialize MST weight
    for (const auto &e : baseEx)
        RT[0].MSTWeight += e.w;                                   // Calculate total weight of the MST
    vector<vector<pair<int, int>>> adjEx = buildExAdj(baseEx, n); // Build adjacency list for the MST edges
    RT[0].Create_Tree(adjEx);                                     // Create the Rooted Tree from the edges in Ex
    for (int b = 1; b < tp; ++b)
    {
        RT[b] = RT[0]; // Copy the initial Rooted Tree to all time blocks
        RT[b].V = n;   // Ensure the number of vertices is set correctly
    }
    int queryId = n + m; // Start query ID after the edges
    for (int qi = 0; qi < q; ++qi)
    {
        int op;
        infile >> op;

        if (op == 0)
        {
            int t;
            infile >> t;
            cout << "[DEBUG] Performing MST query at time " << t << endl;
            vector<Edge> mst = getMST(t, Tc, n, que, RT);
            if (mst.empty())
            {
                cout << "[WARNING] MST is empty at time t = " << t << ".\n";
                continue;
            }
            cout << "MST at time " << t << ":\n";
            double mstWeight = 0.0;
            for (const Edge &e : mst)
            {
                if (e.w < 0 || e.w >= INF)
                    continue; // skip invalid/deleted edges
                mstWeight += e.w;
                cout << " Edge (" << e.x << ", " << e.y << "), w = " << e.w << ", id = " << e.id << "\n";
            }

            cout << "Total MST Weight: " << mstWeight << "\n";
        }
        else if (op == 1)
        {
            int x, y, w, t;
            infile >> x >> y >> w >> t;
            x--;
            y--;
            addAndDeleteEdge(x, y, w, t, queryId++, RT, OperationType::INSERT);
            cout << "[DEBUG] Inserted edge (" << x << ", " << y << ") with weight " << w << " at time " << t << ". Query ID: " << queryId - 1 << "\n";
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
                cout << "[WARNING] Could not find edge (" << x << "-" << y << ", w=" << w << ") for deletion.\n";
            }
            // Perform the deletion and print its single, unique query ID
            addAndDeleteEdge(x, y, w, t, realId, RT, OperationType::DELETE);
            cout << "[DEBUG] Deleted edge (" << x << ", " << y << ") with weight " << w << " at time " << t << ". Query ID: " << queryId - 1 << endl;
            queryId++;
        }
        else
        {
            cout << "[ERROR] Invalid operation type in input: " << op << endl;
        }
    }
    infile.close();
    cout << "[DEBUG] Execution completed. Total queries processed: " << queryId - n - m << endl;
    printf("Execution Time: %.2fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    return 0;
}