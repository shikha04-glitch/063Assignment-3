// campus_planner_user.cpp
// Compile: g++ -std=c++14 campus_planner_user.cpp -O2 -o campus_planner_user

#include <bits/stdc++.h>
using namespace std;

/* ---------------------------
   Building ADT
   --------------------------- */
struct Building {
    int id;
    string name;
    string locationDetails;
    Building(int _id = 0, const string& _name = "", const string& _loc = "") :
        id(_id), name(_name), locationDetails(_loc) {}
};

/* ---------------------------
   Binary Search Tree (BST)
   --------------------------- */
struct BSTNode {
    Building data;
    BSTNode *left, *right;
    BSTNode(const Building& b): data(b), left(nullptr), right(nullptr) {}
};

class BST {
public:
    BST(): root(nullptr) {}
    void insertBuilding(const Building& b) { root = insertRec(root, b); }
    BSTNode* search(int id) { return searchRec(root, id); }
    void inorder() { inorderRec(root); cout << "\n"; }
    void preorder() { preorderRec(root); cout << "\n"; }
    void postorder() { postorderRec(root); cout << "\n"; }
    int height() { return heightRec(root); }
    void clear() { clearRec(root); root = nullptr; }

private:
    BSTNode* root;
    BSTNode* insertRec(BSTNode* node, const Building& b) {
        if (!node) return new BSTNode(b);
        if (b.id < node->data.id) node->left = insertRec(node->left, b);
        else if (b.id > node->data.id) node->right = insertRec(node->right, b);
        else node->data = b; // update
        return node;
    }
    BSTNode* searchRec(BSTNode* node, int id) {
        if (!node) return nullptr;
        if (id == node->data.id) return node;
        if (id < node->data.id) return searchRec(node->left, id);
        return searchRec(node->right, id);
    }
    void inorderRec(BSTNode* node) {
        if (!node) return;
        inorderRec(node->left);
        printNode(node->data);
        inorderRec(node->right);
    }
    void preorderRec(BSTNode* node) {
        if (!node) return;
        printNode(node->data);
        preorderRec(node->left);
        preorderRec(node->right);
    }
    void postorderRec(BSTNode* node) {
        if (!node) return;
        postorderRec(node->left);
        postorderRec(node->right);
        printNode(node->data);
    }
    int heightRec(BSTNode* node) {
        if (!node) return 0;
        return 1 + max(heightRec(node->left), heightRec(node->right));
    }
    void clearRec(BSTNode* node) {
        if (!node) return;
        clearRec(node->left);
        clearRec(node->right);
        delete node;
    }
    static void printNode(const Building& b) {
        cout << "[ID:" << b.id << ", " << b.name << ", " << b.locationDetails << "] ";
    }
};

/* ---------------------------
   AVL Tree
   --------------------------- */
struct AVLNode {
    Building data;
    AVLNode *left, *right;
    int height;
    AVLNode(const Building& b): data(b), left(nullptr), right(nullptr), height(1) {}
};

class AVL {
public:
    AVL(): root(nullptr) {}
    void insertBuilding(const Building& b) { root = insertRec(root, b); }
    int height() { return heightNode(root); }
    void inorder() { inorderRec(root); cout << "\n"; }
    void preorder() { preorderRec(root); cout << "\n"; }
    void postorder() { postorderRec(root); cout << "\n"; }
    void clear() { clearRec(root); root = nullptr; }

private:
    AVLNode* root;

    int heightNode(AVLNode* n) { return n ? n->height : 0; }
    int getBalance(AVLNode* n) { return n ? heightNode(n->left) - heightNode(n->right) : 0; }
    AVLNode* rightRotate(AVLNode* y) {
        AVLNode* x = y->left;
        AVLNode* T2 = x->right;
        x->right = y;
        y->left = T2;
        y->height = 1 + max(heightNode(y->left), heightNode(y->right));
        x->height = 1 + max(heightNode(x->left), heightNode(x->right));
        cout << "Rotation: RightRotate at node " << y->data.id << "\n";
        return x;
    }
    AVLNode* leftRotate(AVLNode* x) {
        AVLNode* y = x->right;
        AVLNode* T2 = y->left;
        y->left = x;
        x->right = T2;
        x->height = 1 + max(heightNode(x->left), heightNode(x->right));
        y->height = 1 + max(heightNode(y->left), heightNode(y->right));
        cout << "Rotation: LeftRotate at node " << x->data.id << "\n";
        return y;
    }

    AVLNode* insertRec(AVLNode* node, const Building& key) {
        if (!node) return new AVLNode(key);
        if (key.id < node->data.id) node->left = insertRec(node->left, key);
        else if (key.id > node->data.id) node->right = insertRec(node->right, key);
        else { node->data = key; return node; } // update

        node->height = 1 + max(heightNode(node->left), heightNode(node->right));
        int balance = getBalance(node);

        // LL
        if (balance > 1 && key.id < node->left->data.id) return rightRotate(node);
        // RR
        if (balance < -1 && key.id > node->right->data.id) return leftRotate(node);
        // LR
        if (balance > 1 && key.id > node->left->data.id) {
            node->left = leftRotate(node->left);
            cout << "Rotation: LeftRotate (LR case) applied\n";
            return rightRotate(node);
        }
        // RL
        if (balance < -1 && key.id < node->right->data.id) {
            node->right = rightRotate(node->right);
            cout << "Rotation: RightRotate (RL case) applied\n";
            return leftRotate(node);
        }
        return node;
    }

    void inorderRec(AVLNode* node) {
        if (!node) return;
        inorderRec(node->left);
        cout << "[ID:" << node->data.id << ", " << node->data.name << "] ";
        inorderRec(node->right);
    }
    void preorderRec(AVLNode* node) {
        if (!node) return;
        cout << "[ID:" << node->data.id << ", " << node->data.name << "] ";
        preorderRec(node->left);
        preorderRec(node->right);
    }
    void postorderRec(AVLNode* node) {
        if (!node) return;
        postorderRec(node->left);
        postorderRec(node->right);
        cout << "[ID:" << node->data.id << ", " << node->data.name << "] ";
    }
    void clearRec(AVLNode* node) {
        if (!node) return;
        clearRec(node->left);
        clearRec(node->right);
        delete node;
    }
};

/* ---------------------------
   Graph (Adjacency list + matrix)
   --------------------------- */
class Graph {
public:
    int n;
    vector<vector<pair<int,int>>> adj; // (neighbor, weight)
    vector<vector<int>> adjMatrix;
    const int INF = 1000000000;

    Graph(): n(0) {}
    void resize(int nodes) {
        n = nodes;
        adj.assign(n, vector<pair<int,int>>());
        adjMatrix.assign(n, vector<int>(n, INF));
        for (int i = 0; i < n; ++i) adjMatrix[i][i] = 0;
    }

    void addEdge(int u, int v, int w=1, bool undirected=true) {
        if (u<0 || u>=n || v<0 || v>=n) return;
        adj[u].push_back(make_pair(v,w));
        adjMatrix[u][v] = min(adjMatrix[u][v], w);
        if (undirected) {
            adj[v].push_back(make_pair(u,w));
            adjMatrix[v][u] = min(adjMatrix[v][u], w);
        }
    }

    vector<int> BFS(int src) {
        vector<int> order;
        if (src<0 || src>=n) return order;
        vector<int> vis(n, 0);
        queue<int> q;
        q.push(src); vis[src]=1;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            order.push_back(u);
            for (size_t i=0;i<adj[u].size();++i) {
                int v = adj[u][i].first;
                if (!vis[v]) { vis[v]=1; q.push(v); }
            }
        }
        return order;
    }

    vector<int> DFS(int src) {
        vector<int> order;
        if (src<0 || src>=n) return order;
        vector<int> vis(n,0);
        stack<int> st;
        st.push(src);
        while (!st.empty()) {
            int u = st.top(); st.pop();
            if (vis[u]) continue;
            vis[u]=1;
            order.push_back(u);
            for (int i = (int)adj[u].size()-1; i>=0; --i) {
                int v = adj[u][i].first;
                if (!vis[v]) st.push(v);
            }
        }
        return order;
    }

    pair<vector<int>, vector<int>> dijkstra(int src) {
        vector<int> dist(n, INF), parent(n, -1);
        if (src<0 || src>=n) return make_pair(dist, parent);
        priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> pq;
        dist[src] = 0;
        pq.push(make_pair(0, src));
        while (!pq.empty()) {
            pair<int,int> top = pq.top(); pq.pop();
            int d = top.first;
            int u = top.second;
            if (d != dist[u]) continue;
            for (size_t i=0;i<adj[u].size();++i) {
                int v = adj[u][i].first;
                int w = adj[u][i].second;
                if (dist[v] > dist[u] + w) {
                    dist[v] = dist[u] + w;
                    parent[v] = u;
                    pq.push(make_pair(dist[v], v));
                }
            }
        }
        return make_pair(dist, parent);
    }

    pair<int, vector<tuple<int,int,int>>> kruskalMST() {
        struct Edge{ int u,v,w; };
        vector<Edge> edges;
        for (int u=0; u<n; ++u) {
            for (size_t i=0;i<adj[u].size();++i) {
                int v = adj[u][i].first;
                int w = adj[u][i].second;
                if (u < v) edges.push_back({u,v,w});
            }
        }
        sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b){ return a.w < b.w; });
        vector<int> parent(n);
        for (int i=0;i<n;++i) parent[i]=i;
        function<int(int)> findp = [&](int x)->int { return parent[x]==x ? x : parent[x]=findp(parent[x]); };
        auto unite = [&](int a,int b){ a=findp(a); b=findp(b); if (a!=b) parent[b]=a; };
        int total = 0;
        vector<tuple<int,int,int>> chosen;
        for (size_t i=0;i<edges.size();++i) {
            int u = edges[i].u, v = edges[i].v, w = edges[i].w;
            if (findp(u) != findp(v)) {
                unite(u, v);
                total += w;
                chosen.push_back(make_tuple(u,v,w));
            }
        }
        return make_pair(total, chosen);
    }

    void printAdjList(const vector<string>& names) {
        cout << "Adjacency List:\n";
        for (int i=0;i<n;++i) {
            cout << i << "(" << (i < (int)names.size() ? names[i] : to_string(i)) << "): ";
            for (size_t j=0;j<adj[i].size();++j) cout << "[" << adj[i][j].first << "," << adj[i][j].second << "] ";
            cout << "\n";
        }
    }
    void printAdjMatrix(const vector<string>& names) {
        cout << "Adjacency Matrix (INF = no-edge):\n   ";
        for (int j=0;j<n;++j) cout << j << "(" << (j < (int)names.size() ? names[j] : to_string(j)) << ") ";
        cout << "\n";
        for (int i=0;i<n;++i) {
            cout << i << "(" << (i < (int)names.size() ? names[i] : to_string(i)) << "): ";
            for (int j=0;j<n;++j) {
                if (adjMatrix[i][j] >= INF/2) cout << "INF ";
                else cout << adjMatrix[i][j] << " ";
            }
            cout << "\n";
        }
    }
};

/* ---------------------------
   Expression Tree (infix->postfix->tree->eval)
   --------------------------- */
struct ExprNode { string val; ExprNode *left, *right; ExprNode(const string& v): val(v), left(nullptr), right(nullptr) {} };

bool isOperator(const string& s) { return s=="+"||s=="-"||s=="*"||s=="/"; }
int precedence(const string& op) { if (op=="+"||op=="-") return 1; if (op=="*"||op=="/") return 2; return 0; }

vector<string> tokenizeExpr(const string& s) {
    vector<string> tokens;
    int n=s.size();
    for (int i=0;i<n;) {
        if (isspace((unsigned char)s[i])) { ++i; continue; }
        if (isdigit((unsigned char)s[i]) || s[i]=='.') {
            string num;
            while (i<n && (isdigit((unsigned char)s[i])||s[i]=='.')) num.push_back(s[i++]);
            tokens.push_back(num);
        } else {
            tokens.push_back(string(1, s[i])); ++i;
        }
    }
    return tokens;
}

vector<string> infixToPostfix(const vector<string>& tokens) {
    vector<string> output;
    stack<string> st;
    for (size_t i=0;i<tokens.size();++i) {
        string tk = tokens[i];
        if (tk.empty()) continue;
        if (isdigit((unsigned char)tk[0]) || (tk.size()>1 && isdigit((unsigned char)tk[1]))) {
            output.push_back(tk);
        } else if (tk == "(") st.push(tk);
        else if (tk == ")") {
            while(!st.empty() && st.top()!="(") { output.push_back(st.top()); st.pop(); }
            if (!st.empty() && st.top()=="(") st.pop();
        } else if (isOperator(tk)) {
            while(!st.empty() && isOperator(st.top()) && precedence(st.top()) >= precedence(tk)) {
                output.push_back(st.top()); st.pop();
            }
            st.push(tk);
        } else {
            st.push(tk);
        }
    }
    while(!st.empty()) { output.push_back(st.top()); st.pop(); }
    return output;
}

ExprNode* buildExprTree(const vector<string>& postfix) {
    stack<ExprNode*> st;
    for (size_t i=0;i<postfix.size();++i) {
        string tk = postfix[i];
        if (isOperator(tk)) {
            ExprNode *r = nullptr, *l = nullptr;
            if (!st.empty()) { r = st.top(); st.pop(); }
            if (!st.empty()) { l = st.top(); st.pop(); }
            ExprNode* node = new ExprNode(tk);
            node->left = l; node->right = r;
            st.push(node);
        } else {
            st.push(new ExprNode(tk));
        }
    }
    if (st.empty()) return nullptr;
    ExprNode* root = st.top(); st.pop();
    return root;
}

double evalExpr(ExprNode* root) {
    if (!root) return 0.0;
    if (!root->left && !root->right) return stod(root->val);
    double L = evalExpr(root->left);
    double R = evalExpr(root->right);
    if (root->val == "+") return L+R;
    if (root->val == "-") return L-R;
    if (root->val == "*") return L*R;
    if (root->val == "/") return R == 0 ? numeric_limits<double>::infinity() : L/R;
    return 0.0;
}

void printExprInorder(ExprNode* node) {
    if (!node) return;
    bool op = isOperator(node->val);
    if (op) cout << "(";
    printExprInorder(node->left);
    cout << node->val;
    printExprInorder(node->right);
    if (op) cout << ")";
}

/* ---------------------------
   Save / Load CSV helpers
   --------------------------- */
void saveBuildingsCSV(const string& filename, const vector<Building>& buildings) {
    ofstream fout(filename);
    if (!fout) { cout << "Cannot open file to write.\n"; return; }
    fout << "id,name,location\n";
    for (size_t i=0;i<buildings.size();++i) {
        fout << buildings[i].id << "," << buildings[i].name << "," << buildings[i].locationDetails << "\n";
    }
    fout.close();
    cout << "Saved buildings to " << filename << "\n";
}

vector<Building> loadBuildingsCSV(const string& filename) {
    vector<Building> res;
    ifstream fin(filename);
    if (!fin) { cout << "Cannot open file to read.\n"; return res; }
    string line;
    getline(fin, line); // header
    while (getline(fin, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        string id_s, name, loc;
        getline(ss, id_s, ',');
        getline(ss, name, ',');
        getline(ss, loc, ',');
        int id = stoi(id_s);
        res.push_back(Building(id, name, loc));
    }
    fin.close();
    cout << "Loaded " << res.size() << " buildings from " << filename << "\n";
    return res;
}

/* ---------------------------
   Interactive Menu & main
   --------------------------- */
void printMainMenu() {
    cout << "\n====== Campus Navigation & Utility Planner (User Menu) ======\n";
    cout << "1. Add building (manually)\n";
    cout << "2. Load buildings from CSV\n";
    cout << "3. Save current buildings to CSV\n";
    cout << "4. List buildings (BST traversals)\n";
    cout << "5. Show BST and AVL heights\n";
    cout << "6. Build/Reset campus graph (use current buildings)\n";
    cout << "7. Add edge to graph\n";
    cout << "8. Show graph (adj list & matrix)\n";
    cout << "9. BFS from a node\n";
    cout << "10. DFS from a node\n";
    cout << "11. Dijkstra shortest path\n";
    cout << "12. Kruskal MST (utility layout)\n";
    cout << "13. Evaluate expression (expression tree)\n";
    cout << "14. Demo dataset (reset to sample)\n";
    cout << "0. Exit\n";
    cout << "============================================================\n";
    cout << "Enter choice: ";
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    vector<Building> buildings;
    BST bst;
    AVL avl;
    Graph graph;
    vector<string> names;

    auto rebuildTreesFromBuildings = [&]() {
        bst.clear();
        avl.clear();
        for (size_t i=0;i<buildings.size();++i) {
            bst.insertBuilding(buildings[i]);
            avl.insertBuilding(buildings[i]);
        }
    };

    auto rebuildGraphNames = [&]() {
        names.clear();
        for (size_t i=0;i<buildings.size();++i) names.push_back(buildings[i].name);
    };

    auto loadDemo = [&]() {
        buildings.clear();
        buildings.push_back(Building(101,"Admin","Administration Block"));
        buildings.push_back(Building(102,"Library","Central Library"));
        buildings.push_back(Building(103,"CSE","Computer Science Dept"));
        buildings.push_back(Building(104,"ECE","Electronics Dept"));
        buildings.push_back(Building(105,"HostelA","Student Hostel A"));
        buildings.push_back(Building(106,"Cafeteria","Main Cafeteria"));
        buildings.push_back(Building(107,"Auditorium","Main Auditorium"));

        rebuildTreesFromBuildings();
        graph.resize((int)buildings.size());
        rebuildGraphNames();

        // demo edges
        graph.addEdge(0,1,50);
        graph.addEdge(1,2,80);
        graph.addEdge(2,3,120);
        graph.addEdge(2,5,60);
        graph.addEdge(5,4,90);
        graph.addEdge(3,6,200);
        graph.addEdge(1,5,140);
        graph.addEdge(0,2,160);

        cout << "Demo dataset loaded.\n";
    };

    // initially load demo
    loadDemo();

    while (true) {
        printMainMenu();
        int choice;
        if (!(cin >> choice)) { cout << "Invalid input, exiting.\n"; break; }
        if (choice == 0) { cout << "Exiting. Bye!\n"; break; }

        if (choice == 1) {
            cout << "Enter Building ID (int): "; int id; cin >> id; cin.ignore();
            cout << "Enter Building name: "; string name; getline(cin, name);
            cout << "Enter Location details: "; string loc; getline(cin, loc);
            buildings.push_back(Building(id, name, loc));
            rebuildTreesFromBuildings();
            rebuildGraphNames();
            cout << "Building added.\n";
        } else if (choice == 2) {
            cout << "Enter CSV filename to load: ";
            string fn; cin >> fn;
            vector<Building> loaded = loadBuildingsCSV(fn);
            if (!loaded.empty()) {
                buildings = loaded;
                rebuildTreesFromBuildings();
                graph.resize((int)buildings.size());
                rebuildGraphNames();
            }
        } else if (choice == 3) {
            cout << "Enter CSV filename to save: ";
            string fn; cin >> fn;
            saveBuildingsCSV(fn, buildings);
        } else if (choice == 4) {
            cout << "BST Inorder: ";
            bst.inorder();
            cout << "BST Preorder: ";
            bst.preorder();
            cout << "BST Postorder: ";
            bst.postorder();
        } else if (choice == 5) {
            cout << "BST height: " << bst.height() << "\n";
            cout << "AVL height: " << avl.height() << "\n";
        } else if (choice == 6) {
            int n = (int)buildings.size();
            graph.resize(n);
            rebuildGraphNames();
            cout << "Graph resized to " << n << " nodes (use option 7 to add edges or 14 to reload demo).\n";
        } else if (choice == 7) {
            if (graph.n == 0) { cout << "Graph is empty. Build graph first (option 6 or 14).\n"; continue; }
            cout << "Enter u v weight (indices 0.." << graph.n-1 << "): ";
            int u,v,w; cin >> u >> v >> w;
            graph.addEdge(u,v,w);
            cout << "Edge added.\n";
        } else if (choice == 8) {
            if (graph.n==0) { cout << "Graph empty.\n"; continue; }
            graph.printAdjList(names);
            cout << "\n";
            graph.printAdjMatrix(names);
        } else if (choice == 9) {
            if (graph.n==0) { cout << "Graph empty.\n"; continue; }
            cout << "Enter source index (0.." << graph.n-1 << "): ";
            int s; cin >> s;
            vector<int> ord = graph.BFS(s);
            cout << "BFS order: ";
            for (size_t i=0;i<ord.size();++i) cout << ord[i] << "("<<names[ord[i]]<<") ";
            cout << "\n";
        } else if (choice == 10) {
            if (graph.n==0) { cout << "Graph empty.\n"; continue; }
            cout << "Enter source index (0.." << graph.n-1 << "): ";
            int s; cin >> s;
            vector<int> ord = graph.DFS(s);
            cout << "DFS order: ";
            for (size_t i=0;i<ord.size();++i) cout << ord[i] << "("<<names[ord[i]]<<") ";
            cout << "\n";
        } else if (choice == 11) {
            if (graph.n==0) { cout << "Graph empty.\n"; continue; }
            cout << "Enter source index and target index (0.." << graph.n-1 << "): ";
            int s,t; cin >> s >> t;
            pair<vector<int>, vector<int>> dp = graph.dijkstra(s);
            vector<int> dist = dp.first;
            vector<int> parent = dp.second;
            if (dist[t] >= graph.INF) cout << "No path from " << s << " to " << t << "\n";
            else {
                cout << "Distance = " << dist[t] << "\nPath: ";
                vector<int> path;
                for (int cur = t; cur != -1; cur = parent[cur]) path.push_back(cur);
                reverse(path.begin(), path.end());
                for (size_t i=0;i<path.size();++i) cout << path[i] << "("<<names[path[i]]<<") ";
                cout << "\n";
            }
        } else if (choice == 12) {
            if (graph.n==0) { cout << "Graph empty.\n"; continue; }
            pair<int, vector<tuple<int,int,int>>> mst = graph.kruskalMST();
            int total = mst.first;
            vector<tuple<int,int,int>> edges = mst.second;
            cout << "Total MST weight = " << total << "\nEdges:\n";
            for (size_t i=0;i<edges.size();++i) {
                int u = get<0>(edges[i]);
                int v = get<1>(edges[i]);
                int w = get<2>(edges[i]);
                cout << u << "("<<names[u]<<") - " << v << "("<<names[v]<<") : " << w << "\n";
            }
        } else if (choice == 13) {
            cout << "Enter arithmetic expression (e.g., 100 + 50*(2+3) - 25/5):\n";
            cin.ignore();
            string line;
            getline(cin, line);
            vector<string> toks = tokenizeExpr(line);
            vector<string> postfix = infixToPostfix(toks);
            cout << "Postfix: ";
            for (size_t i=0;i<postfix.size();++i) cout << postfix[i] << " ";
            cout << "\n";
            ExprNode* root = buildExprTree(postfix);
            cout << "Inorder: ";
            printExprInorder(root);
            cout << "\n";
            cout << "Result = " << evalExpr(root) << "\n";
        } else if (choice == 14) {
            loadDemo();
        } else {
            cout << "Invalid option.\n";
        }
    }

    return 0;
}
