#include <bits/stdc++.h> 
using namespace std;
const double INF = 1e15;

struct Edge
{
    int x, y, id; // Nodes x and y connected by this edge, id is the index of the edge
    double w; // Weight of the edge
    // replacedId indicates which changed edge (by index) might replace this edge
    // -1 meaning no replacement
    int replacedId = -1; 
    Edge() : x(-1), y(-1), w(0.0), id(-1), replacedId(-1) {}//constructor
    Edge(int u, int v, double w, int id_) : x(u), y(v), w(w), id(id_), replacedId(-1) {}
    bool operator==(const Edge &other) const
        {
            return (x == other.x && y== other.y && w == other.w && replacedId == other.replacedId);
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
    NONE=0, // added to remainder
    INS=1,  // inserted to key edges
    DEL=2,  // deleted from key edges
    RPL=3   // replaced an existing key edge
};

enum class OperationType {
    INSERT = 0,
    DELETE = 1
};


int findMinHeightVertex(vector<vector<pair<int, int>>> &adjEx)
{
    {
    int V = adjEx.size();
    vector<int> degree(V, 0); // Degree of each vertex
    queue<int> q;             // Queue to store leaves

    // Initialize the degree array and find initial leaves
    for (int i = 0; i < V; i++)
    {
        degree[i] = adjEx[i].size();
        if (degree[i] == 1)
        {
            q.push(i); // Add leaves to the queue
        }
    }

    // Trim leaves level by level
    int remainingVertices = V;
    while (remainingVertices > 2)
    {
        int sz = q.size();
        remainingVertices -= sz; // Remove current leaves

        for (int i = 0; i < sz; i++)
        {
            int leaf = q.front();
            q.pop();

            // Reduce the degree of neighbors
            for (const auto &[neighbor, w] : adjEx[leaf])
            {
                degree[neighbor]--;
                if (degree[neighbor] == 1)
                {
                    q.push(neighbor); // Add new leaf
                }
            }
        }
    }

    // The remaining vertices are the centers
    if (q.size() == 1)
    {
        return q.front();
    }
    else if (q.size() == 2)
    {
        int first = q.front();
        q.pop();
        int second = q.front();

        // Return either vertex (degree comparison is unnecessary as both are balanced)
        return first;
    }

    return -1; // Should not reach here for a valid tree
}
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
    std::unordered_set<int> deletedIds; // Track deleted edge ids
    RootedTree() : V(0) {}
    RootedTree(int V)
    {
        this->V = V;
        RootedT.resize(V);
    }
    RootedTree(const RootedTree& other)
{
    V = other.V;
    RootedT = other.RootedT;
    Ex = other.Ex;
    Er = other.Er;
    deletedIds = other.deletedIds;
}

    

    //We Create the Rooted Tree (Algorithm 3)
    // We find a vertex with Minimum Height vertex as root and run BFS from it.
    // If multiple components exist, we repeat.

        
    //BFS begins. Use parallelisaztion for BFS
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
    ///BFS ends.
    
    
    void Create_Tree(vector<vector<pair<int, int>>> &adjEx)
    {
        V = adjEx.size();

        for (int i = 0; i < V; i++)
        {
            if (RootedT[i].parent == -1)
            {
                int r = findMinHeightVertex(adjEx);
                RootedT[r].parent = r;
                RootedT[r].root = r;
                // BFS from r
                bfsFromRoot(r, adjEx);
                break;
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
    for (int i = 0; i < (int)CE.size(); i++)
    {
        Status[i] = EdgeStatus::NONE;
        Marked[i] = Edge();
        Edge &E = CE[i];
        bool op = operation[i];

        if (op)
        {
            // op == true now means “INSERT a new edge”
            int ru = findRoot(E.x), rv = findRoot(E.y);
            if (ru != rv)
            {
                // Endpoints in different components → insert
                Status[i] = EdgeStatus::INS;
            }
            else
            {
                // Same component → maybe replace a heavier edge
                Edge MaxW = findMaxEdgeOnPath(E.x, E.y);
                if (MaxW.id == -1)
                {
                    cout << "[DEBUG] No replaceable edge found: MaxW.id == -1" << endl;
                }
                else if (MaxW.w > E.w && MaxW.id!=-1)
                {
                    // New edge is lighter → replace
                    Status[i] = EdgeStatus::RPL;
                    Marked[i] = MaxW;

                    // Find MaxW’s index in Ex
                    auto it = std::find_if(Ex.begin(), Ex.end(),
                                           [&](const Edge &e) { return e.id == MaxW.id; });
                    if (it != Ex.end()) {
                        Marked[i].replacedId = std::distance(Ex.begin(), it);
                    } else {
                        Marked[i].replacedId = -1;
                        cerr << "[ERROR] Could not locate edge id="
                             << MaxW.id << " in Ex!" << endl;
                    }

                    cout << "[DEBUG] Marking edge for replacement: "
                         << "MaxW.id = " << MaxW.id
                         << ", rpl_idx = " << Marked[i].replacedId
                         << ", E.id = " << E.id
                         << ", MaxW.w = " << MaxW.w
                         << ", E.w = " << E.w << endl;
                }
                else
                {
                    cout << "[DEBUG] Found heavier or equal edge on path, won't replace: "
                         << "MaxW.w = " << MaxW.w
                         << ", NewEdge.w = " << E.w << endl;
                }
            }
        }
        else
        {
            // op == false now means “DELETE the edge”
            Status[i] = EdgeStatus::DEL;
        }
    }
}




void Process_Status(vector<Edge> &CE, vector<EdgeStatus> &Status, vector<Edge> &Marked)
{
    std::cout << "[DEBUG] → Process_Status: entering, CE.size="
              << CE.size() << ", Status.size=" << Status.size()
              << ", Marked.size=" << Marked.size() << "\n";

    if (Status.size() != CE.size() || Marked.size() != CE.size()) {
        std::cerr << "[FATAL] Size mismatch in Process_Status" << std::endl;
        exit(1);
    }

    bool rebuildtree = false;

    std::unordered_map<int, int> edgeVisitCount;

    for (int i = 0; i < (int)CE.size(); ++i)
    {
        Edge &E = CE[i];
        EdgeStatus S = Status[i];

        edgeVisitCount[E.id]++;
        if (edgeVisitCount[E.id] > 5) {
            std::cerr << "[LOOP WARNING] Edge id=" << E.id << " processed > 5 times.\n";
        }

        std::cout << "[DEBUG] Edge #" << i << " status: " << static_cast<int>(S)
                  << " (" << E.x << "-" << E.y << ", w=" << E.w << ", id=" << E.id << ")" << std::endl;

        if (deletedIds.count(E.id)) {
            std::cout << "[DEBUG] Edge id=" << E.id << " already processed. Skipping." << std::endl;
            continue;
        }

        switch (S)
        {
        case EdgeStatus::DEL:
        {
            std::cout << "[DEBUG] Edge marked for deletion" << std::endl;
            bool deleted = false;

            auto itEx = std::remove_if(Ex.begin(), Ex.end(), [&](const Edge &e) {
                return e.id == E.id;
            });

            if (itEx != Ex.end()) {
                Ex.erase(itEx, Ex.end());
                rebuildtree = true;
                deleted = true;
                std::cout << "[DEBUG] Deleted edge from Ex: id=" << E.id << std::endl;
            }

            auto itEr = std::remove_if(Er.begin(), Er.end(), [&](const Edge &e) {
                return e.id == E.id;
            });

            if (itEr != Er.end()) {
                Er.erase(itEr, Er.end());
                deleted = true;
                std::cout << "[DEBUG] Deleted edge from Er: id=" << E.id << std::endl;
            }

            if (!deleted) {
                std::cerr << "[WARNING] No matching edge found for deletion: ("
                          << E.x << "-" << E.y << ", w=" << E.w << ", id=" << E.id << ")" << std::endl;
            }

            deletedIds.insert(E.id);
            break;
        }

        case EdgeStatus::NONE:
            Er.push_back(E);
            std::cout << "[DEBUG] Edge moved to remainder edges (Er)" << std::endl;
            break;

        case EdgeStatus::INS:
        {
            bool exists = std::any_of(Ex.begin(), Ex.end(), [&](const Edge &e) {
                return e.id == E.id;
            });

            if (exists) {
                std::cout << "[DEBUG] Edge with id=" << E.id << " already exists in Ex. Skipping insertion.\n";
            } else {
                Ex.push_back(E);
                rebuildtree = true;
                std::cout << "[DEBUG] Edge inserted into Ex\n";
            }

            // Remove this edge from Er to prevent reprocessing
            Er.erase(std::remove_if(Er.begin(), Er.end(), [&](const Edge &e) {
                return deletedIds.count(e.id) || e.w == INF || e.w == -1;
            }), Er.end());

            break;
        }

        case EdgeStatus::RPL:
{
    Edge &Erpl = Marked[i];

    std::cout << "[DEBUG] Edge marked for replacement, original id = " << Erpl.id << std::endl;

    if (deletedIds.count(Erpl.id)) {
        std::cout << "[DEBUG] Replacement target already removed: id=" << Erpl.id << std::endl;

        // Avoid inserting the new edge if old one is already gone
        Er.erase(std::remove_if(Er.begin(), Er.end(), [&](const Edge &e) {
            return e.id == E.id;
        }), Er.end());
        continue;
    }

    auto it = std::find_if(Ex.begin(), Ex.end(), [&](const Edge &e) {
        return e.id == Erpl.id;
    });

    if (it != Ex.end()) {
        Ex.erase(it);
        std::cout << "[DEBUG] Replaced edge in Ex: id=" << Erpl.id << std::endl;
        deletedIds.insert(Erpl.id);
        rebuildtree = true;
    } else {
        std::cerr << "[WARNING] RPL: Edge to be replaced not found in Ex: id=" << Erpl.id << std::endl;
    }

    // Now insert the new edge
    auto alreadyPresent = std::any_of(Ex.begin(), Ex.end(), [&](const Edge &e) {
        return e.id == E.id;
    });
    if (!alreadyPresent) {
        Ex.push_back(E);
        rebuildtree = true;
        std::cout << "[DEBUG] Inserted replacement edge into Ex: id=" << E.id << std::endl;
    }

    // Remove from Er
    Er.erase(std::remove_if(Er.begin(), Er.end(), [&](const Edge &e) {
        return e.id == E.id;
    }), Er.end());

    break;
}

        }
    }

    std::cout << "[DEBUG] ← Process_Status: exiting, CE.size=" << CE.size()
              << ", Ex.size=" << Ex.size() << "\n";

    CE.clear();

    Ex.erase(std::remove_if(Ex.begin(), Ex.end(), [](const Edge &e) {
        return e.w == INF || e.w == -1;
    }), Ex.end());

    Er.erase(std::remove_if(Er.begin(), Er.end(), [](const Edge &e) {
        return e.w == INF || e.w == -1;
    }), Er.end());

    // Deduplicate Er
    std::unordered_set<int> seenIds;
    Er.erase(std::remove_if(Er.begin(), Er.end(), [&](const Edge &e) {
        if (seenIds.count(e.id)) return true;
        seenIds.insert(e.id);
        return false;
    }), Er.end());

    // Remove edges marked as deleted
    Er.erase(std::remove_if(Er.begin(), Er.end(), [&](const Edge &e) {
        return deletedIds.count(e.id);
    }), Er.end());

    if (rebuildtree)
    {
        std::cout << "[DEBUG] Rebuilding tree due to changes in Ex" << std::endl;
        vector<vector<pair<int, int>>> adjEx = buildExAdj(Ex, V);
        Create_Tree(adjEx);
    }
    else
    {
        std::cout << "[DEBUG] No changes in Ex, no need to rebuild tree" << std::endl;
    }

    if (Ex.size() < V - 1 && !Er.empty()) {
    std::cout << "[DEBUG] Incomplete MST detected. Attempting repair.\n";
    Repair_Tree(Er, Status, Marked);
}

    std::cout << "[DEBUG] Process_Status completed, CE.size = "
              << CE.size() << ", Ex.size = " << Ex.size() << std::endl;
}




void Repair_Tree(vector<Edge> &Er, vector<EdgeStatus> &Status, vector<Edge> &Marked)
{
    std::cout << "[DEBUG] Repair_Tree called" << std::endl;

    bool changed;
    int loopCount = 0;

    do {
        changed = false;

        Status.resize(Er.size(), EdgeStatus::NONE);
        Marked.resize(Er.size());

        for (int i = 0; i < (int)Er.size(); i++) {
            Edge &E = Er[i];
            int ru = findRoot(E.x), rv = findRoot(E.y);
            if (ru != rv) {
                Status[i] = EdgeStatus::INS;
                changed = true;
                std::cout << "[DEBUG] Repair candidate: " << E.x << "-" << E.y
                          << " (w=" << E.w << ", id=" << E.id << ")\n";
            }
        }

        if (changed) {
            Er.erase(std::remove_if(Er.begin(), Er.end(), [&](const Edge &e) {
                return deletedIds.count(e.id);
            }), Er.end());

            Process_Status(Er, Status, Marked);

            // Remove processed edges from Er
            Er.erase(std::remove_if(Er.begin(), Er.end(), [&](const Edge &e) {
                return deletedIds.count(e.id) ||
                       std::any_of(Ex.begin(), Ex.end(), [&](const Edge &ex) { return ex.id == e.id; });
            }), Er.end());
        }

        // Fallback repair
        if (!changed && Ex.size() < V - 1 && !Er.empty()) {
            std::cout << "[DEBUG] Fallback repair: attempting to reconnect disconnected components.\n";
            std::sort(Er.begin(), Er.end());  // Try lightest edge first
            for (auto it = Er.begin(); it != Er.end(); ++it) {
                Edge &e = *it;
                int ru = findRoot(e.x), rv = findRoot(e.y);
                if (ru != rv) {
                    std::cout << "[DEBUG] Attempting fallback insert: "
                              << e.x << "-" << e.y << ", w=" << e.w << ", id=" << e.id << "\n";
                    Ex.push_back(e);
                    deletedIds.insert(e.id);
                    Er.erase(it);  // Remove used edge
                    vector<vector<pair<int, int>>> adjEx = buildExAdj(Ex, V);
                    Create_Tree(adjEx);
                    changed = true;
                    std::cout << "[DEBUG] Successfully repaired by fallback insert: id=" << e.id << "\n";
                    break;
                }
            }
        }

        loopCount++;
        if (loopCount > 20) {
            std::cerr << "[FATAL] Repair_Tree exceeded loop limit. Likely infinite loop. Breaking.\n";
            break;
        }

    } while (changed);

    std::cout << "[DEBUG] Repair_Tree completed\n";
}



void ProcessAllEdges(vector<Edge> &CE, vector<bool> &operation)
{
    std::cout << "[DEBUG] ProcessAllEdges called\n";

    if (CE.empty()) {
        std::cout << "[DEBUG] No edges to process\n";
        return;
    }

    vector<EdgeStatus> Status;
    vector<Edge> Marked;

    while (!CE.empty())
    {
        std::cout << "[DEBUG] Processing CE loop. CE.size = " << CE.size() << "\n";

        // Step 1: Assign default statuses
        Status.assign(CE.size(), EdgeStatus::NONE);
        Marked.assign(CE.size(), Edge());

        // Step 2: Classify edges (no operation[] here!)
        Classify_Edges(CE, Status, Marked, operation);   // correct for your code

        // Step 3: Filter out deleted edges before processing
        CE.erase(std::remove_if(CE.begin(), CE.end(), [&](const Edge &e) {
            return deletedIds.count(e.id);
        }), CE.end());

        // Keep auxiliary vectors in sync
        operation.resize(CE.size());
        Status.resize(CE.size());
        Marked.resize(CE.size());

        // Step 4: Process status changes
        Process_Status(CE, Status, Marked);

        // Step 5: Remove from CE any edges that were:
        // - Deleted
        // - Already inserted into Ex
        CE.erase(std::remove_if(CE.begin(), CE.end(), [&](const Edge &e) {
            return deletedIds.count(e.id) ||
                   std::any_of(Ex.begin(), Ex.end(), [&](const Edge &ex) {
                       return ex.id == e.id;
                   });
        }), CE.end());

        // Step 6: Resize 'operation' to match CE size again
        operation.resize(CE.size());
    }

    // Final sanity check
    if (Ex.size() > V - 1) {
        std::cerr << "[ERROR] MST has too many edges: Ex.size = " << Ex.size()
                  << ", expected = " << (V - 1) << std::endl;
    }

    std::cout << "[DEBUG] ProcessAllEdges completed. Ex.size = " << Ex.size() << "\n";
}





    void updateNode(int id, int x, int y, int w, OperationType op)
    {
        cout << "[DEBUG] updateNode called with id = " << id << ", x = " << x << ", y = " << y << ", w = " << w << ", op = " << (op == OperationType::INSERT ? "INSERT" : "DELETE")<< endl;
        if (deletedIds.count(id)) {
    std::cout << "[DEBUG] Edge id=" << id << " already deleted. Skipping update.\n";
    return;
}
        
        for (const auto &e : Ex) {
            if (e.id == id && (e.w == INF || e.w == -1)) {
                cout << "[WARNING] Attempt to update already deleted edge id = " << id << endl;
                return; // Skip processing
            }
        }

        if (op==OperationType::INSERT) { // insertion
        for (const auto &e : Ex) {
            if (e.id == id) {
                cout << "[INFO] Edge already exists in Ex, skipping insertion: id = " << id << endl;
                return;
            }
        }
    }
        // Create a new edge with the given parameters
        vector<Edge> CE(1, Edge(x, y, w, id));
        vector<bool> operation(1, op==OperationType::INSERT); // true for insert, false for delete
        ProcessAllEdges(CE, operation);
    }

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
        if (rootX == rootY) return false;

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
    sort(edges.begin(), edges.end()); // Requires Edge::operator<

    DSU dsu(n);
    vector<Edge> Ex, Er;
    int edgeCounter = 0;

    for (Edge &edge : edges)
    {
        edge.id = edgeCounter++;
        if (dsu.union_sets(edge.x, edge.y))
            Ex.push_back(edge);
        else
            Er.push_back(edge);
    }

    return {Ex, Er};
}


const int T = 50011; //Timeline size
const int Tc = sqrt(T) + 10; //Timeline size for block decomposition (+10 for optimisation and incase N is small +10 reduces complications)

// Global variables
int tp;             //Key Blocks of Timeline
int n;              // Number of nodes
int m;              // Number of edges
int q;              // Number of queries
vector<vector<pair<Edge, OperationType>>> que;


vector<Edge> getMST(int t, vector<RootedTree> &rt)
{
    cout << "[DEBUG] getMST called with t = " << t << endl;

    // Time boundary checks
    if (t < 0)
    {
        cout << "[DEBUG] Time t is negative. Returning empty MST." << endl;
        return vector<Edge>();
    }

    if (t >= tp * Tc)
    {
        cout << "[DEBUG] Time t exceeds maximum range (tp * M). Returning empty MST." << endl;
        return vector<Edge>();
    }

    if (t < Tc)
    {
        cout << "[DEBUG] Time t is before first modification block. Returning initial MST (rt[0])" << endl;
        return rt[0].Ex;
    }

    int bkt = t / Tc;
    cout << "[DEBUG] Processing block #" << bkt << " (bkt = " << bkt << ")" << endl;

    if (bkt >= rt.size())
    {
        cout << "[ERROR] Block index out of range for rt. Returning empty MST." << endl;
        return vector<Edge>();
    }

    vector<Edge> CE; // Edge changes to process
    vector<bool> operation;

    for (int i = bkt * Tc + 1; i <= t && i < que.size(); ++i)
    {
        if (!que[i].empty())
        {
            cout << "[DEBUG] Processing que[" << i << "] with " << que[i].size() << " updates." << endl;
            for (int j = 0; j < que[i].size(); ++j)
            {
                const auto &entry = que[i][j];
                CE.push_back(entry.first);
                operation.push_back(entry.second == OperationType::INSERT);
                cout << "  [DEBUG] Scheduled " << (entry.second == OperationType::INSERT ? "INSERT" : "DELETE")
                     << " of edge (" << entry.first.x << ", " << entry.first.y
                     << ", w=" << entry.first.w << ")" << endl;
            }
        }
    }

    cout << "[DEBUG] Total edges to process: " << CE.size() << endl;

    rt[bkt].ProcessAllEdges(CE, operation);
    vector<Edge> mst = rt[bkt].Ex;  // ← fetch the MST after processing

    cout << "[DEBUG] getMST returning MST with " << mst.size() << " edges." << endl;
    return mst;
}


 void addAndDeleteEdge(int x, int y, int w, int t, int id,
                      vector<RootedTree> &rt, OperationType op)
{
    cout << "[DEBUG] addAndDeleteEdge called with x=" << x
         << ", y=" << y << ", w=" << w << ", t=" << t
         << ", id=" << id << ", op=" << (op == OperationType::INSERT ? "INSERT" : "DELETE") << endl;

    if (t < 0 || t >= T)
    {
        cout << "[ERROR] Time t is out of bounds. Skipping." << endl;
        return;
    }

    if ((int)que.size() <= t)
        que.resize(T);  // Ensure que is large enough

    // Schedule for the query-time algorithm
    que[t].push_back({ Edge(x, y, w, id), op });

    int block = t / Tc;
    for (int b = block; b < (int)rt.size(); ++b)
    {
        rt[b].updateNode(id, x, y, w, op);
    }

    cout << "[DEBUG] addAndDeleteEdge completed\n";
}




int main()
{
    // Start measuring time
    clock_t tStart = clock();
    ///freopen("debug_output.txt", "w", stdout);
    cout << "Timer started: " << tStart << endl;

    // Read input from file
    std::ifstream infile("C:\\Personal\\Studies\\Projects\\Parallel and Concurrent DS\\retro_rooted-main\\in000");
    if (!infile.is_open()) {
        cerr << "Error opening input file." << endl;
        return 1;
    }

    infile >> n >> m >> q;
    cout << "[DEBUG] Read n = " << n << ", m = " << m << ", q = " << q << endl;

    vector<Edge> edges(m);
    tp = (T / Tc) + 2;
    vector<RootedTree> RT(tp, RootedTree(n));

    // Load original edges and build base MST in RT[0]
    for (int i = 0; i < m; ++i) {
        int x, y, w;
        infile >> x >> y >> w;
        edges[i] = Edge(--x, --y, w, i);
    }
    auto [baseEx, baseEr] = findMST(edges, n);
    RT[0].Ex = baseEx;
    RT[0].Er = baseEr;

    // Build adjacency and rooted‐tree for block 0
    vector<vector<pair<int, int>>> adjEx_0 = buildExAdj(baseEx, n);
    RT[0].Create_Tree(adjEx_0);

    // Copy block 0’s state into all subsequent blocks
    for (int b = 1; b < tp; ++b) {
        RT[b] = RootedTree(RT[b - 1]);
    }

    int queryId = n + m;
    const int totalQueries = 5;

    for (int qi = 0; qi < totalQueries; ++qi) {
        int op;
        cout << "\nEnter operation (0 = MST query, 1 = Add Edge, 2 = Delete Edge): ";
        cin >> op;

        if (op == 0) {
            // — MST query at time t
            int t;
            cout << "Enter time t to query MST: ";
            cin >> t;
            if (t < 0 || t >= T) {
                cout << "[ERROR] Invalid time t = " << t << ". Must be in [0, " << T - 1 << "]" << endl;
                --qi;
                continue;
            }
            cout << "[DEBUG] Performing MST query at time " << t << endl;
            vector<Edge> mst = getMST(t, RT);

            double sum = 0.0;
            for (auto &e : mst) {
                sum += e.w;
                cout << "[MST] Edge: " << e.x << "-" << e.y << ", w=" << e.w << ", id=" << e.id << endl;
            }
            cout << "MST weight: " << sum << endl;
        }
        else if (op == 1 || op == 2) {
            // — Add or Delete an edge
            int x, y, w, t;
            cout << "Enter x y weight w and time t: ";
            cin >> x >> y >> w >> t;
            if (t < 0 || t >= T) {
                cout << "[ERROR] Invalid time t = " << t << ". Must be in [0, " << T - 1 << "]" << endl;
                --qi;
                continue;
            }

            x--; 
            y--;
            bool isInsert = (op == 1);
            int realId = queryId;

            if (!isInsert) {
                // — Deletion: scan every block RT[i] to find matching (x,y,w)
                bool found = false;
                for (int b = 0; b < tp && !found; ++b) {
                    // scan Ex in block b
                    for (auto &e : RT[b].Ex) {
                        if (((e.x == x && e.y == y) || (e.x == y && e.y == x)) && e.w == w) {
                            realId = e.id;
                            found = true;
                            break;
                        }
                    }
                    if (found) break;
                    // scan Er in block b
                    for (auto &e : RT[b].Er) {
                        if (((e.x == x && e.y == y) || (e.x == y && e.y == x)) && e.w == w) {
                            realId = e.id;
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    cout << "[WARNING] Could not find edge (" 
                         << x << "-" << y << ", w=" << w << ") in any RT[].Ex or RT[].Er." << endl;
                }
            }

            // Schedule and apply the update
            OperationType opType = isInsert ? OperationType::INSERT : OperationType::DELETE;
            addAndDeleteEdge(x, y, w, t, realId, RT, opType);
            queryId++;
        }
        else {
            cout << "[ERROR] Invalid operation type. Must be 0, 1, or 2." << endl;
            --qi;
        }
    }

    printf("Execution Time: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    return 0;
}






/*Sample Input:
20 40 10
12 16 4033
6 16 1870
6 4 4967
14 15 7277
9 8 9040
17 16 6216
16 9 1356
5 13 2536
8 14 9615
16 14 3932
18 11 8003
1 11 7326
8 16 698
5 7 7158
3 5 7012
19 17 8996
16 3 850
9 2 1978
10 14 4745
9 3 1704
19 18 9530
19 13 907
9 14 5792
15 7 903
9 16 1963
1 6 4548
1 17 284
8 15 8742
19 20 1921
9 17 6078
5 7 341
13 16 9306
2 17 7976
11 12 5738
15 19 1519
6 4 8740
9 8 799
14 20 7629
16 11 6776
5 12 282
0 7
0 1
0 5
1 16 17 7541 10
0 4
1 10 1 7867 9
1 17 15 9635 6
0 3
0 8
0 2
*/