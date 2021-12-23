#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <string>
#include <stack>
#include <algorithm>
using namespace std;
const int inf = 2000000000;
class Graph{

private:
    int m_numberOfNodes;
    int m_numberOfEdges;

    vector <vector<int>> m_costMatrix;

public:
    Graph(int numberOfNodes, int numberOfEdges, vector <vector<int>> &costMatrix){
        m_numberOfNodes = numberOfNodes;
        m_numberOfEdges = numberOfEdges;
        m_costMatrix = costMatrix;
        //costMatrix.clear();
    }

    virtual ~Graph() {
        m_costMatrix.clear();
    }

    int getNumberOfNodes() const {
        return m_numberOfNodes;
    }

    void setNumberOfNodes(int mNumberOfNodes) {
        m_numberOfNodes = mNumberOfNodes;
    }

    int getNumberOfEdges() const {
        return m_numberOfEdges;
    }

    void setNumberOfEdges(int mNumberOfEdges) {
        m_numberOfEdges = mNumberOfEdges;
    }

    const vector<vector<int>> &getCostMatrix() const {
        return m_costMatrix;
    }

    void setCostMatrix(const vector<vector<int>> &mCostMatrix) {
        m_costMatrix = mCostMatrix;
    }
    void setEdgeCost(int i, int j, int value){
        if(i <= m_costMatrix.size()){
            if(j <= m_costMatrix[i].size()){
                m_costMatrix[i][j] = value;
            }
        }
    }
    //problema royfloyd de pe infoarena
    void royFloydWarshall(){
        for(int k = 1; k <= m_numberOfNodes; ++k){
            for(int i = 1; i <= m_numberOfNodes; ++i){
                for(int j = 1; j <= m_numberOfNodes; ++j){
                    m_costMatrix[i][j] = min(m_costMatrix[i][j], m_costMatrix[i][k] + m_costMatrix[k][j]);
                }
            }
        }
    }
    void bfs(int start, vector <bool> &viz){
        queue <int> q;
        viz[start] = true;
        q.push(start);
        while(!q.empty()){
            int current = q.front();
            q.pop();
            unsigned int sz = m_costMatrix[current].size();
            for(unsigned int i = 0; i < sz; ++i){
                int neighbour = m_costMatrix[current][i];
                if(!viz[neighbour]){
                    viz[neighbour] = true;
                    q.push(neighbour);
                }
            }
        }
    }

    void dfs(int start, vector <bool> &viz){
        viz[start] = true;
        unsigned int sz = m_costMatrix[start].size();
        for(unsigned int i = 0; i < sz; ++i){
            int neighbour = m_costMatrix[start][i];
            if(!viz[neighbour]){
                dfs(neighbour, viz);
            }
        }
    }

    void dfsStack(int start, vector <int> &viz, stack <int> &s){
        viz[start] = 1;
        unsigned int sz = m_costMatrix[start].size();
        for(unsigned int i = 0; i < sz; ++i){
            int neighbour = m_costMatrix[start][i];
            if(!viz[neighbour]){
                dfsStack(neighbour, viz, s);
            }
        }
        s.push(start);
    }



    int getNumberOfConnectedComponents(){
        vector <bool> viz (m_numberOfNodes + 1, false);
        if(m_numberOfNodes == 0){
            return 0;
        }
        int cc = 0;
        for(int i = 1; i <= m_numberOfNodes; ++i){
            if(!viz[i]){
                ++cc;
                dfs(i, viz);
            }
        }
        return cc;
    }

    int bfs_dist(int start, vector <bool> &viz, vector <int> &dist){
        viz.resize(m_numberOfNodes + 1);
        dist.resize(m_numberOfNodes + 1);
        for(int i = 0; i <= m_numberOfNodes; ++i){
            viz[i] = false;
            dist[i] = 0;
        }
        queue <int> q;
        q.push(start);
        viz[start] = true;
        int lastVisitedNode = start;
        dist[start] = 1;
        while(!q.empty()){
            int current = q.front();
            unsigned int neighbours = m_costMatrix[current].size();
            q.pop();
            for(int i = 0; i < neighbours; ++i){
                int neighbour = m_costMatrix[current][i];
                if(!viz[neighbour]){
                    viz[neighbour] = true;
                    dist[neighbour] = dist[current] + 1;
                    q.push(neighbour);
                    lastVisitedNode = neighbour;
                }
            }
        }
        return lastVisitedNode;
    }

    void dfsTopSort(int start, vector <bool> &isVisited, vector <int> &order) {
        isVisited[start] = true;
        unsigned int sz = m_costMatrix[start].size();
        for(unsigned int i = 0; i < sz; ++i){
            int neighbour = m_costMatrix[start][i];
            if(!isVisited[neighbour]){
                isVisited[neighbour] = true;
                dfsTopSort(neighbour, isVisited, order);
            }
        }
        order.push_back(start);
    }

    vector <int> topologicalSort(){
        vector <bool> isVisited (m_numberOfNodes + 1, false);
        vector <int> order;
        for(int i = 1; i <= m_numberOfNodes; ++i){
            if(!isVisited[i]){
                dfsTopSort(i, isVisited, order);
            }
        }
        reverse(order.begin(), order.end());
        return order;
    }

    int getDiameter(){
        vector <bool> viz;
        vector <int> dist;
        int lastVisitedNode = bfs_dist(1, viz, dist);
        int k = bfs_dist(lastVisitedNode, viz, dist);
        return dist[k];
    }

    void dijkstra(int start, vector <vector < pair <int, int>>> &adjList, vector <int> &dist){
        vector <bool> viz(m_numberOfNodes + 1, false);
        dist.clear();
        dist.resize(m_numberOfNodes + 1, inf);
        priority_queue <pair <int,int>, vector <pair <int, int>>, greater < pair <int, int>>> q;
        dist[start] = 0;
        q.push({0, start});
        while(!q.empty()){
            int current = q.top().second;
            q.pop();
            if(!viz[current]){
                viz[current] = true;
                unsigned int sz = adjList[current].size();
                for(unsigned int i = 0 ; i < sz; ++i){
                    int neighbour = adjList[current][i].first;
                    int cost = adjList[current][i].second;
                    if(dist[current] + cost < dist[neighbour]){
                        dist[neighbour] = dist[current] + cost;
                        q.push({dist[neighbour], neighbour});
                    }
                }
            }
        }

    }

    bool hasEulerianPath(){
        //strict pt problema de pe infoarena, ciclueuler
        for(int i = 1; i <= m_numberOfNodes; ++i){
            if(m_costMatrix[i].size() % 2){
                return false;
            }
        }
        return true;
    }

};

class disjointSet{
private:
    int m_numberOfElements;
    vector <int> m_parent;
    vector <int> m_ranks;
public:
    explicit disjointSet(int numberOfElements){
        m_numberOfElements = numberOfElements;
        for(int i = 0; i <= numberOfElements; ++i){
            m_parent.push_back(i);
            m_ranks.push_back(0);
        }
    }

    virtual ~disjointSet() = default;

    int getNumberOfElements() const {
        return m_numberOfElements;
    }

    void setNumberOfElements(int mNumberOfElements) {
        m_numberOfElements = mNumberOfElements;
    }

    const vector<int> &getParent() const {
        return m_parent;
    }

    void setParent(const vector<int> &mParent) {
        m_parent = mParent;
    }


    int findParent(int i){
        int root = i;
        while(m_parent[root] != root){
            root = m_parent[root];
        }
        //path compression
        while(i != root){
            int next = m_parent[i];
            m_parent[i] = root;
            i = next;
        }

        return root;
    }

    void unify(int x, int y){
        int px = findParent(x);
        int py = findParent(y);
        if(px != py){
            if(m_ranks[px] > m_ranks[py]){
                m_parent[py] = px;
            }
            else if(m_ranks[px] < m_ranks[py]){
                m_parent[px] = py;
            }
            else{
                m_parent[py] = px;
                m_ranks[px] ++;
            }
        }
    }

};

class minimumSpanningTree{
private:
    int m_numberOfNodes;
    int m_numberOfEdges;
    int m_cost = 0;
    vector <pair <int, pair<int, int>>> m_edges;
    vector <pair <int, int>> m_mst;
    int m_mstEdges = 0;
    disjointSet m_dsj = disjointSet(0);
public:
    minimumSpanningTree(int mNumberOfNodes, int mNumberOfEdges, const vector<pair<int, pair<int, int>>> &mEdges){
        m_numberOfNodes = mNumberOfNodes;
        m_numberOfEdges = mNumberOfEdges;
        m_edges = mEdges;
        m_dsj = disjointSet(mNumberOfNodes);
        createMst();

    }

    virtual ~minimumSpanningTree() {
    }

    int getNumberOfNodes() const {
        return m_numberOfNodes;
    }

    void setNumberOfNodes(int mNumberOfNodes) {
        m_numberOfNodes = mNumberOfNodes;
    }

    int getNumberOfEdges() const {
        return m_numberOfEdges;
    }

    void setNumberOfEdges(int mNumberOfEdges) {
        m_numberOfEdges = mNumberOfEdges;
    }

    int getCost() const {
        return m_cost;
    }
    void createMst(){
        sort(m_edges.begin(), m_edges.end());
        unsigned int sz = m_edges.size();
        for(int i = 0; i < sz; ++i){
            int a, b, c;
            //extragem din m_edges nodul sursa, destinatia si costul
            a = m_edges[i].second.first;
            b = m_edges[i].second.second;
            c = m_edges[i].first;
            //daca nodurile nu se afla in aceeasi multime (componenta conexa)
            if(m_dsj.findParent(a) != m_dsj.findParent(b)){
                m_dsj.unify(a, b);
                m_cost += c;
                m_mstEdges += 1;
                m_mst.push_back(make_pair(a, b));
            }
            if(m_mstEdges == m_numberOfNodes - 1){
                break;
            }

        }
    }
    void setCost(int mCost) {
        m_cost = mCost;
    }

    int getMstEdges() const {
        return m_mstEdges;
    }

    const vector<pair<int, pair<int, int>>> &getEdges() const {
        return m_edges;
    }

    void setEdges(const vector<pair<int, pair<int, int>>> &mEdges) {
        m_edges = mEdges;
    }

    const vector<pair<int, int>> &getMst() const {
        return m_mst;
    }

    void setMst(const vector<pair<int, int>> &mMst) {
        m_mst = mMst;
    }
};


void pbRoyFloid(){
    //problema royfloyd de pe infoarena
    ifstream f("royfloyd.in");
    ofstream g("royfloyd.out");
    int n, m = 0;
    f >> n;

    vector <vector<int>> v;
    v.resize(n + 1, vector <int> (n + 1, 0));

    for(int i = 1; i <= n; ++i){
        for(int j = 1; j <= n; ++j){
            f >> v[i][j];
            ++m;
            if(i != j && v[i][j] == 0){
                v[i][j] = inf;
            }
        }
    }
    Graph graf(n, m, v);
    graf.royFloydWarshall();
    v = graf.getCostMatrix();
    for(int i = 1; i <= n; ++i){
        for(int j = 1; j <= n; ++j){
            if(v[i][j] != inf){
                g << v[i][j] << " ";
            }
            else{
                g << 0 << " ";
            }
        }
        g << "\n";
    }
}

void pbDarb(const string& source, const string& output){
    ifstream f(source);
    ofstream g(output);
    int n, a, b, m;
    f >> n;
    //cout << n;
    vector <vector <int>> v (n + 1);
    for(int i = 1; i < n; ++i){
        f >> a >> b;
        v[a].push_back(b);
        v[b].push_back(a);
        m += 2;
    }
    Graph graph(n, m, v);
    g << graph.getDiameter();
    f.close();
    g.close();
}

void pbSortareTopologica(const string& source, const string& output){
    //problema de pe infoarena
    ifstream f(source);
    ofstream g(output);
    int n, m;
    f >> n >> m;
    vector <vector <int>> v(n + 1);
    for(int i = 0; i < m; ++i){
        int x, y;
        f >> x >> y;
        v[x].push_back(y);
    }
    Graph graph(n, m, v);
    vector <int> order = graph.topologicalSort();

    for(unsigned int i = 0; i < order.size(); ++i){
        g << order[i] << " ";
    }
    f.close();
    g.close();
}
void pbDfs(const string &source, const string &output){
    ifstream f(source);
    ofstream g(output);
    int n, m;
    f >> n >> m;
    vector <vector <int>> v(n + 1);
    for(int i = 0; i < m; ++i){
        int x, y;
        f >> x >> y;
        v[x].push_back(y);
        v[y].push_back(x);
    }
    Graph graph(n, m, v);
    g << graph.getNumberOfConnectedComponents();
    f.close();
    g.close();
}

void pbDisJoint(const string &source, const string &output){
    ifstream f(source);
    ofstream g(output);
    int n, m, question, a, b;
    f >> n >> m;
    disjointSet dsj(n + 1);
    for(int i = 0; i < m; ++i){
        f >> question >> a >> b;
        if(question == 1){
            dsj.unify(a, b);
        }
        else if(question == 2){
            if(dsj.findParent(a) != dsj.findParent(b)){
                g << "NU\n";
            }
            else{
                g << "DA\n";
            }
        }
    }
}

void pbApm(const string& source, const string& output){
    ifstream f (source);
    ofstream g (output);
    int n, m, a, b, c;
    vector <pair <int, pair <int, int>>> edges;
    vector <pair <int, int>> apm;
    f >> n >> m;
    for(int i = 0; i < m; ++i){
        f >> a >> b >> c;
        edges.push_back(make_pair(c, make_pair(a, b)));
    }
    minimumSpanningTree tree1(n, m, edges);
    apm = tree1.getMst();
    g << tree1.getCost() << "\n";
    g << tree1.getMstEdges() << "\n";
    unsigned int sz = apm.size();
    for(unsigned int i = 0; i < sz; ++i){
        g << apm[i].first << " " << apm[i].second << "\n";
    }

    f.close();
    g.close();
}
void pbDijkstra(const string &source, const string &output){
    ifstream f(source);
    ofstream g(output);
    int n, m;
    f >> n >> m;
    vector <vector <pair <int, int>>> v(n + 1);
    vector <vector <int>> k(n + 1);
    vector <int> dist(n + 1, inf);
    for(int i = 0; i < m; ++i){
        int a, b, c;
        f >> a >> b >> c;
        v[a].push_back(make_pair(b, c));
        k[a].push_back(b);
    }
    Graph graph(n, m, k);
    graph.dijkstra(1, v, dist);
    for(int i = 2; i <= n; ++i){
        if(dist[i] == inf){
            g << 0 << " ";
        }
        else{
            g << dist[i] << " ";
        }
    }
    f.close();
    g.close();

}
int main() {
    pbDijkstra("dijkstra.in", "dijkstra.out");

    return 0;
}
