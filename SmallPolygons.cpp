#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cassert>
#include <limits.h>
#include <cfloat>
#include <string>
#include <string.h>
#include <sstream>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stack>
#include <queue>

using namespace std;

typedef long long ll;

const int MAX_NP = 1600;
const int UNKNOWN = -1;
const int COUNTER_CLOCKWISE = 1;
const int CLOCKWISE = -1;
const int ONLINE_BACK = 2;
const int ONLINE_FRONT = -2;
const int ON_SEGMENT = 0;
const double EPS = 1e-10;

int pointCount;
int pointY[MAX_NP];
int pointX[MAX_NP];
double pointsDistance[MAX_NP][MAX_NP];
int pointUsedCount[MAX_NP];
int nodeCount;

string int2string(int number){
  stringstream ss; 
  ss << number;
  return ss.str();
}

class Vector{
  public:
  int id;
  double x;
  double y;

  Vector(double y = 0.0, double x = 0.0){
    this->y = y;
    this->x = x;
  }

  bool operator==(const Vector& v) const{
    return (x == v.x && y == v.y);
  }

  bool operator<(const Vector& v) const{
    return (x != v.x)? x < v.x : y < v.y;
  }

  Vector operator+(Vector p){
    return Vector(y + p.y, x + p.x);
  }

  Vector operator-(Vector p){
    return Vector(y - p.y, x - p.x);
  }

  Vector operator*(double k){
    return Vector(k * y, k * x);
  }

  Vector operator/(double k){
    return Vector(y / k, x / k);
  }

  double norm(){
    return x*x + y*y;
  }

  double abs(){
    return sqrt(norm());
  }
};

double norm(Vector a){
  return a.x*a.x + a.y*a.y;
}

double abs(Vector a){
  return sqrt(norm(a));
}

struct Point {
  Vector *p1;
  Vector *p2;
  Vector *p3;
  double dist;

  Point(Vector *p1, Vector *p2, Vector *p3, double dist){
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
    this->dist = dist;
  }

  bool operator >(const Point &e) const{
    return dist > e.dist;
  }    
};

struct Circle{
  Vector center;    // 中心座標
  double radius;    // 半径
};

class Triangle{
  public:
    const Vector *p1, *p2, *p3;   // 頂点座標

    bool operator==(const Triangle& t) const{
      return (*p1 == *t.p1 && *p2 == *t.p2 && *p3 == *t.p3) ||
        (*p1 == *t.p1 && *p2 == *t.p3 && *p3 == *t.p2) ||
        (*p1 == *t.p2 && *p2 == *t.p3 && *p3 == *t.p1) ||
        (*p1 == *t.p2 && *p2 == *t.p1 && *p3 == *t.p3) ||
        (*p1 == *t.p3 && *p2 == *t.p1 && *p3 == *t.p2) ||
        (*p1 == *t.p3 && *p2 == *t.p2 && *p3 == *t.p1);
    }

    bool operator<(const Triangle& t) const{
      return !(*getMinVertex() == *t.getMinVertex())? 
        *getMinVertex() < *t.getMinVertex() :
        !(*getMidVertex() == *t.getMidVertex())?
        *getMidVertex() < *t.getMidVertex() :
        *getMaxVertex() < *t.getMaxVertex();
    }

    bool hasCommonPoints(const Triangle& t) const{
      return *p1 == *t.p1 || *p1 == *t.p2 || *p1 == *t.p3 ||
        *p2 == *t.p1 || *p2 == *t.p2 || *p2 == *t.p3 ||
        *p3 == *t.p1 || *p3 == *t.p2 || *p3 == *t.p3;
    }

  private:
    inline const Vector* getMinVertex() const{
      return *p1 < *p2 && *p1 < *p3 ? p1 : (*p2 < *p3)? p2 : p3;
    }

    inline const Vector* getMidVertex() const{
      return ((*p1 < *p2 && *p2 < *p3) || (*p3 < *p2 && *p2 < *p1))? p2 :
        ((*p2 < *p3 && *p3 < *p1) || (*p1 < *p3 && *p3 < *p2))? p3 : p1;
    }

    inline const Vector* getMaxVertex() const{
      return (*p2 < *p1 && *p3 < *p1)? p1 : (*p3 < *p2)? p2 : p3;
    }
};

Vector vectorList[10000];

class Delaunay2d{
  public:
    typedef const set<Vector>               ConstVertexSet;
    typedef ConstVertexSet::const_iterator  ConstVertexIter;

    typedef set<Triangle>                   TriangleSet;
    typedef set<Triangle>::iterator         TriangleIter;

    typedef map<Triangle, bool>             TriangleMap;

    static void getDelaunayTriangles(ConstVertexSet &pVertexSet, TriangleSet *triangleSet){
      Triangle hugeTriangle;{
        double maxX, maxY; maxX = maxY = DBL_MIN;
        double minX, minY; minX = minY = DBL_MAX;

        for(ConstVertexIter it = pVertexSet.begin(); it != pVertexSet.end(); it++){
          double y = it->y;
          double x = it->x;

          maxX = max(maxX, x);
          minX = min(minX, x);

          maxY = max(maxY, y);
          minY = min(minY, y);
        }

        double centerX = (maxX - minX) * 0.5;   // 中心x座標
        double centerY = (maxY - minY) * 0.5;   // 中心y座標

        double dx = maxX - centerX;
        double dy = maxY - centerY;
        double radius = sqrt(dx*dx + dy*dy) + 1.0; // 半径

        Vector *p1 = &vectorList[MAX_NP-1];
        p1->x = centerX - sqrt(3.0) * radius;
        p1->y = centerY - radius;

        Vector *p2 = &vectorList[MAX_NP-2];
        p2->x = centerX + sqrt(3.0) * radius;
        p2->y = centerY - radius;

        Vector *p3 = &vectorList[MAX_NP-3];
        p3->x = centerX;
        p3->y = centerY + 2.0 * radius;

        hugeTriangle.p1 = p1;
        hugeTriangle.p2 = p2;
        hugeTriangle.p3 = p3;
      }

      triangleSet->insert(hugeTriangle);

      for(ConstVertexIter vIter = pVertexSet.begin(); vIter != pVertexSet.end(); vIter++){
        const Vector *p = &*vIter;

        TriangleMap rddcMap;

        for(TriangleIter tIter = triangleSet->begin(); tIter != triangleSet->end();){
          Triangle t = *tIter;

          Circle c;{
            double x1 = t.p1->x; double y1 = t.p1->y;
            double x2 = t.p2->x; double y2 = t.p2->y;
            double x3 = t.p3->x; double y3 = t.p3->y;

            double m = 2.0 * ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
            double x = ((y3-y1)*(x2*x2-x1*x1+y2*y2-y1*y1)
                      +(y1-y2)*(x3*x3-x1*x1+y3*y3-y1*y1))/m;
            double y = ((x1-x3)*(x2*x2-x1*x1+y2*y2-y1*y1)
                      +(x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1))/m;

            c.center.x = x;
            c.center.y = y;

            double dx = t.p1->x - x;
            double dy = t.p1->y - y;
            double radius = sqrt(dx*dx + dy*dy);

            c.radius = radius;
          }

          double dx = c.center.x - p->x;
          double dy = c.center.y - p->y;
          double dist = sqrt(dx*dx + dy*dy);

          if(dist < c.radius){
            Triangle t1;
            t1.p1 = p; t1.p2 = t.p1; t1.p3 = t.p2;
            addElementToRedundanciesMap(&rddcMap, t1);

            Triangle t2;
            t2.p1 = p; t2.p2 = t.p2; t2.p3 = t.p3;
            addElementToRedundanciesMap(&rddcMap, t2);

            Triangle t3;
            t3.p1 = p; t3.p2 = t.p3; t3.p3 = t.p1;
            addElementToRedundanciesMap(&rddcMap, t3);

            triangleSet->erase(tIter++);
          }else{
            tIter++;
          }
        }

        for(TriangleMap::iterator iter = rddcMap.begin(); iter != rddcMap.end(); iter++){
          if(iter->second){
            triangleSet->insert(iter->first);
          }
        }
      }

      for(TriangleIter tIter = triangleSet->begin(); tIter != triangleSet->end();){
        if(hugeTriangle.hasCommonPoints(*tIter)){
          triangleSet->erase(tIter++);
        }else{
          tIter++;
        }
      }
    };

  private:
    static inline void addElementToRedundanciesMap(TriangleMap *pRddcMap, Triangle &t){
      TriangleMap::iterator it = pRddcMap->find(t);

      if(it != pRddcMap->end() && it->second){
        pRddcMap->erase(it);
        pRddcMap->insert(TriangleMap::value_type(t, false));
      }else{
        pRddcMap->insert(TriangleMap::value_type(t, true));
      }
    }
};

struct Node{
  int id;                       // ノードのID
  int degree;                   // ノードの次数
  int originDegree;
  int triangleId;               // このノードに対応している三角形
  set<int> neighbors;           // 隣接ノード一覧
  set<int> originNeighbors;     // 元の隣接ノード
  const Vector *p1, *p2, *p3;   // 頂点座標
  double y;                     // ノードのY座標
  double x;                     // ノードのX座標
  double area;                  // 面積
  bool used;                    // 使用済みかどうか
  bool removed;                 // 削除されたかどうか

  Node(int nodeId = UNKNOWN){
    this->id = nodeId;
    this->used = false;
    this->removed = false;
    this->degree = 0;
  }

  /*
   * 隣接ノードの追加、次数が1増える
   */
  void addNeighbor(int nodeId){
    neighbors.insert(nodeId);
    degree = neighbors.size();
  }

  void addOriginNeighbor(int nodeId){
    originNeighbors.insert(nodeId);
    originDegree = originNeighbors.size();
  }

  /*
   * 隣接ノードの削除、次数が1減る
   */
  void removeNeighbor(int nodeId){
    neighbors.erase(nodeId);
    degree = neighbors.size();
  }

	void clearNeighbor(){
		neighbors.clear();
		degree = 0;
	}

  bool isNeighbor(Node *node){
    int sameCnt = 0;
    sameCnt += (p1 == node->p1 || p1 == node->p2 || p1 == node->p3);
    sameCnt += (p2 == node->p1 || p2 == node->p2 || p2 == node->p3);
    sameCnt += (p3 == node->p1 || p3 == node->p2 || p3 == node->p3);

    return sameCnt == 2;
  }

  bool operator >(const Node &e) const{
    return degree * area > e.degree * e.area;
  }    
};

struct Edge{
  int u;
  int v;
  double cost;

  Edge(int u = UNKNOWN, int v = UNKNOWN, double cost = UNKNOWN){
    this->u    = u;
    this->v    = v;
    this->cost = cost;
  }

  bool operator >(const Edge &e) const{
    return cost > e.cost;
  }  
};

bool comp(const Edge &e1, const Edge &e2){
  return e1.cost < e2.cost;
}

vector<Edge> edgeList;

struct Tree{
  int id;
  set<int> nodes;   // ノードの一覧
  int nodeCount;    // ノードの数
};

// ノード一覧
Node nodeList[10000];

int unionPar[MAX_NP];
int unionRank[MAX_NP];

// n要素で初期化
void initUnionFind(int n){
  for(int i = 0; i < n; i++){
    unionPar[i] = i;
    unionRank[i] = 0;
  }
}

// 木の根を求める
int unionFind(int x){
  if(unionPar[x] == x){
    return x;
  }else{
    return unionPar[x] = unionFind(unionPar[x]);
  }
}

// xとyの属する集合を併合
void unionUnite(int y, int x){
  x = unionFind(x);
  y = unionFind(y);
  if(x == y) return;

  if(unionRank[x] < unionRank[y]){
    unionPar[x] = y;
  }else{
    unionPar[y] = x;
    if(unionRank[x] == unionRank[y]) unionRank[x] += 1;
  }
}

bool unionSame(int y, int x){
  return unionFind(y) == unionFind(x);
}

vector<int> lines;

class SmallPolygons{
  public:
    void init(vector<int> &points){

      for(int id = 0; id < pointCount; id++){
        pointY[id] = points[id*2];
        pointX[id] = points[id*2+1];
        Vector p;
        p.id = id;
        p.y = points[id*2];
        p.x = points[id*2+1];

        vectorList[id] = p;
      }

      initializePointDistance();
    };

    inline Node *getNode(int id){
      return &nodeList[id];
    }

    inline Vector *getVector(int id){
      return &vectorList[id];
    }

    void prim(){
      memset(pointUsedCount, 0, sizeof(pointUsedCount));

      set<int> nodeIdList;
      priority_queue< Edge, vector<Edge>, greater<Edge> > pque;

      int rootId = 0;

      Node *root = getNode(rootId);

      while(root->removed){
        rootId += 1;
        root = getNode(rootId);
      }

      set<int>::iterator it = root->originNeighbors.begin();

      while(it != root->originNeighbors.end()){
        Node *to = getNode(*it);
        Edge edge(root->id, to->id, to->area);

        if(!to->removed){
          //fprintf(stderr,"%d <---> %d\n", root->id, to->id);
          pque.push(Edge(root->id, to->id, to->area));
        }else if(to->removed){
          //fprintf(stderr,"%d <-x-> %d\n", root->id, to->id);
        }

        it++;
      }

      int cnt = 0;

      while(!pque.empty()){
        cnt++;
        Edge e = pque.top(); pque.pop();

        Node *from = getNode(e.u);
        Node *to   = getNode(e.v);

        if(pointUsedCount[to->p1->id] > 0 && pointUsedCount[to->p2->id] > 0 && pointUsedCount[to->p3->id] > 0){
          fprintf(stderr,"%d -> %d: No new point\n", from->id, to->id);
          continue;
        }

        fprintf(stderr,"%d: %d <---> %d\n", cnt, e.u, e.v);

        if(from->removed){
          fprintf(stderr,"(prim)Node %d is removed!\n", from->id);
          continue;
        }

        set<int>::iterator it  = nodeIdList.find(from->id);
        set<int>::iterator iti = nodeIdList.find(to->id);

        if(it == nodeIdList.end()){
          pointUsedCount[from->p1->id] += 1;
          pointUsedCount[from->p2->id] += 1;
          pointUsedCount[from->p3->id] += 1;
        }

        if(iti == nodeIdList.end()){
          pointUsedCount[to->p1->id] += 1;
          pointUsedCount[to->p2->id] += 1;
          pointUsedCount[to->p3->id] += 1;
        }

        from->addNeighbor(to->id);
        to->addNeighbor(from->id);
        nodeIdList.insert(from->id);
        nodeIdList.insert(to->id);

        set<int>::iterator that = to->originNeighbors.begin();

        while(that != to->originNeighbors.end()){
          Node *next = getNode(*that);

          if(!next->removed){
            //fprintf(stderr,"Add edge %d -> %d\n", to->id, next->id);
            pque.push(Edge(to->id, next->id, next->area));
          }
          that++;
        }
      }
    }

    void kruskal(){
      fprintf(stderr,"kruskal =>\n");
      sort(edgeList.begin(), edgeList.end(), comp);
      initUnionFind(nodeCount);

      int E = edgeList.size();
      fprintf(stderr,"Edge count = %d\n", E);
      set<int> nodeIdList;

      for(int i = 0; i < E; i++){
        Edge e = edgeList[i];

        Node *from = getNode(e.u);

        if(!unionSame(e.u, e.v)){
          Node *to = getNode(e.v);

          if(from->removed || to->removed) continue;
          if((pointUsedCount[from->p1->id] > 0 && pointUsedCount[from->p2->id] > 0 && pointUsedCount[from->p3->id] > 0) &&
            (pointUsedCount[to->p1->id] > 0 && pointUsedCount[to->p2->id] > 0 && pointUsedCount[to->p3->id] > 0)){
            //continue;
          }

          unionUnite(e.u, e.v);

          from->addNeighbor(to->id);
          to->addNeighbor(from->id);
          nodeIdList.insert(from->id);
          nodeIdList.insert(to->id);
          
          //fprintf(stderr,"node %d <-> node %d\n", e.u, e.v);
        }
      }

      set<int>::iterator it = nodeIdList.begin();

      while(it != nodeIdList.end()){
        Node *node = getNode((*it));

        pointUsedCount[node->p1->id] += 1;
        pointUsedCount[node->p2->id] += 1;
        pointUsedCount[node->p3->id] += 1;
        it++;
      }
    }

    /*
     * 三角形の面積を計算
     */
    double calcArea(const Vector *p1, const Vector *p2, const Vector *p3){
      double dx, dy;

      dy = (p1->y - p2->y);
      dx = (p1->x - p2->x);

      double a = sqrt(dy*dy + dx*dx);

      dy = (p2->y - p3->y);
      dx = (p2->x - p3->x);

      double b = sqrt(dy*dy + dx*dx);

      dy = (p3->y - p1->y);
      dx = (p3->x - p1->x);

      double c = sqrt(dy*dy + dx*dx);

      double s = (a + b + c) / 2.0;

      double S = sqrt(s * (s - a) * (s - b) * (s - c));

      return S;
    }

    void initializePointDistance(){
      for(int i = 0; i < pointCount; i++){
        int p1_y = pointY[i];
        int p1_x = pointX[i];

        for(int j = i+1; j < pointCount; j++){
          int p2_y = pointY[j];
          int p2_x = pointX[j];

          int dy = p2_y - p1_y;
          int dx = p2_x - p1_x;

          double dist = sqrt(dy*dy + dx*dx);

          pointsDistance[i][j] = dist;
          pointsDistance[j][i] = dist;
        }
      }
    }

    void insertVertex(const Vector *p1, const Vector *p2, const Vector *p, vector<int> &vlist){
      vector<int>::iterator first = vlist.begin();
      int listSize = vlist.size();
      int p1index = find(vlist.begin(), vlist.end(), p1->id) - vlist.begin();
      int p2index = find(vlist.begin(), vlist.end(), p2->id) - vlist.begin();

      int right = max(p1index, p2index);
      int left  = min(p1index, p2index);

      if(left == 0 && right == listSize-1){
        vlist.insert(first + right + 1, p->id);
      }else{
        vlist.insert(first + right, p->id);
      }
    }

    string createPolygon(int nodeId){
			fprintf(stderr,"create polygon = %d\n", nodeId);
      vector<int> vlist;
      queue<int> que;
      map<int, bool> checkList;

      que.push(nodeId);
      int p1index, p2index, p3index;
      int listSize;

      while(!que.empty()){
        int id = que.front(); que.pop();

        if(checkList[id]) continue;
        checkList[id] = true;

        Node *node = getNode(id);
        node->used = true;

        listSize = vlist.size();

        if(listSize == 0){
          vlist.push_back(node->p1->id);
          vlist.push_back(node->p2->id);
          vlist.push_back(node->p3->id);
        }else{
          p1index = find(vlist.begin(), vlist.end(), node->p1->id) - vlist.begin();
          p2index = find(vlist.begin(), vlist.end(), node->p2->id) - vlist.begin();
          p3index = find(vlist.begin(), vlist.end(), node->p3->id) - vlist.begin();

          if(p1index < vlist.size() && p2index < vlist.size() && p3index < vlist.size()){
            fprintf(stderr,"p1 = %d, p2 = %d, p3 = %d\n", node->p1->id, node->p2->id, node->p3->id);
          }

          if(p1index >= listSize){
            insertVertex(node->p2, node->p3, node->p1, vlist);
          }else if(p2index >= listSize){
            insertVertex(node->p1, node->p3, node->p2, vlist);
          }else if(p3index >= listSize){
            insertVertex(node->p1, node->p2, node->p3, vlist);
          }
        }

        set<int>::iterator it = node->neighbors.begin();

        while(it != node->neighbors.end()){
          int nid = (*it);

          if(!checkList[nid]){
            que.push(nid);
          }

          it++;
        }
      }

      string result = "";
      listSize = vlist.size();

      int fixCount = edgeListClearSimple(vlist);
      //int fixCount = edgeListClear(vlist);
      //int fixCount = 0;
      int limit = 0;
      int i = 0;

      while(i < limit){
        if(i%2){
          fixCount = edgeListClearSimple(vlist);
        }else{
          fixCount = edgeListClear(vlist);
        }
        i += 1;
      }

      lines = vlist;

      return polygon2string(vlist);
    }

    string polygon2string(vector<int> &lines){
      int listSize = lines.size();
      string result = "";

      for(int i = 0; i < listSize; i++){
        result += int2string(lines[i]);

        if(i != listSize-1){
          result += " ";
        }
      }

      return result;
    }

    bool lineCross(vector<int> &vlist){
      int listSize = vlist.size();

      for(int i = 0; i < listSize; i++){
        Vector *p1 = &vectorList[vlist[i%listSize]];
        Vector *p2 = &vectorList[vlist[(i+1)%listSize]];


        for(int j = i+1; j < i+listSize; j++){
          Vector *p3 = &vectorList[vlist[j%listSize]];
          Vector *p4 = &vectorList[vlist[(j+1)%listSize]];

          if(intersect2(*p1, *p2, *p3, *p4)){
            return true;
          }
        }
      }

      return false;
    }

    int edgeListClearSimple(vector<int> &vlist){
      int listSize = vlist.size();
      int fixCount = 0;

      for(int i = 0; i < listSize; i++){
        Vector *p1 = &vectorList[vlist[i%listSize]];
        Vector *p2 = &vectorList[vlist[(i+1)%listSize]];
        Vector *p3 = &vectorList[vlist[(i+2)%listSize]];
        Vector *p4 = &vectorList[vlist[(i+3)%listSize]];

        if(intersect2(*p1, *p2, *p3, *p4)){
          fprintf(stderr,"intersect! %d <-> %d, %d <-> %d\n", p1->id, p2->id, p3->id, p4->id);
          swap(vlist[(i+1)%listSize],  vlist[(i+2)%listSize]);
          fixCount += 1;
        }
      }

      fprintf(stderr,"fixCount = %d\n", fixCount);

      return fixCount;
    }

    int edgeListClear(vector<int> &vlist){
      int listSize = vlist.size();
      int fixCount = 0;

      for(int i = 0; i < listSize; i++){
        Vector *p1 = &vectorList[vlist[i%listSize]];
        Vector *p2 = &vectorList[vlist[(i+1)%listSize]];

        double minDist = DBL_MAX;
        int swapId = -1;

        for(int j = i+1; j < i+listSize; j++){
          Vector *p3 = &vectorList[vlist[j%listSize]];
          Vector *p4 = &vectorList[vlist[(j+1)%listSize]];

          if(intersect2(*p1, *p2, *p3, *p4)){
          //if(intersect(p1, p2, p3, p4)){
            
            double d1 = pointsDistance[p1->id][p3->id] + pointsDistance[p2->id][p4->id];
            double d2 = pointsDistance[p1->id][p4->id] + pointsDistance[p3->id][p2->id];
            double d3 = pointsDistance[p3->id][p2->id] + pointsDistance[p1->id][p4->id];
            double d4 = pointsDistance[p4->id][p2->id] + pointsDistance[p3->id][p1->id];
          
            fprintf(stderr,"intersect! %d <-> %d, %d <-> %d\n", p1->id, p2->id, p3->id, p4->id);

            if(d4 > d1 && d3 > d1 && d2 > d1){
              swap(vlist[(i+1)%listSize],  vlist[j%listSize]);
            }else if(d4 > d2 && d3 > d2 && d1 > d2){
              swap(vlist[(i+1)%listSize],  vlist[(j+1)%listSize]);
            }else if(d4 > d3 && d2 > d3 && d1 > d3){
              swap(vlist[i%listSize],  vlist[j%listSize]);
            }else{
              swap(vlist[i%listSize],  vlist[(j+1)%listSize]);
            }
            fixCount += 1;

            /*
            if(minDist > dist){
              fprintf(stderr,"update minDist\n");
              minDist = dist;
              swapId = j;
            }
            */
          }
        }

        if(minDist < DBL_MAX){
          Vector *p3 = &vectorList[vlist[swapId%listSize]];
          Vector *p4 = &vectorList[vlist[(swapId+1)%listSize]];

          fprintf(stderr,"intersect! %d <-> %d, %d <-> %d\n", p1->id, p2->id, p3->id, p4->id);
          double d0 = pointsDistance[p1->id][p2->id] + pointsDistance[p3->id][p4->id];
          double d1 = pointsDistance[p1->id][p3->id] + pointsDistance[p2->id][p4->id];
          double d2 = pointsDistance[p1->id][p4->id] + pointsDistance[p3->id][p2->id];
          double d3 = pointsDistance[p3->id][p2->id] + pointsDistance[p1->id][p4->id];
          double d4 = pointsDistance[p4->id][p2->id] + pointsDistance[p3->id][p1->id];
          
          if(d4 > d1 && d3 > d1 && d2 > d1 && d0 > d1){
            swap(vlist[(i+1)%listSize],  vlist[swapId%listSize]);
          }else if(d4 > d2 && d3 > d2 && d1 > d2 && d0 > d2){
            swap(vlist[(i+1)%listSize],  vlist[(swapId+1)%listSize]);
          }else if(d4 > d3 && d2 > d3 && d1 > d3 && d0 > d3){
            swap(vlist[i%listSize],  vlist[swapId%listSize]);
          }else if(d0 > d4){
            swap(vlist[i%listSize],  vlist[(swapId+1)%listSize]);
          }
          fixCount += 1;

          /*
          swap(vlist[i%listSize],  vlist[swapId%listSize]);
          fixCount += 1;
          */

          break;
        }
      }

      fprintf(stderr,"fixCount = %d\n", fixCount);

      return fixCount;
    }

    int direction(Vector *p0, Vector *p1, Vector *p2){
      return (p1->x - p0->x) * (p2->y - p0->y) - (p2->x - p0->x) * (p1->y - p0->y);
    }

    bool onSegment(Vector *pi, Vector *pj, Vector *pk){
      if((min(pi->x, pj->x) <= pk->x && pk->x <= max(pi->x, pj->x)) && (min(pi->y, pj->y) <= pk->y && pk->y <= max(pi->y, pj->y))){
        return true;
      }else{
        return false;
      }
    }

    // ベクトルaとbの内積
    double dot(Vector a, Vector b){
      return a.x * b.x + a.y * b.y;
    }

    // ベクトルaとbの外積
    double cross(Vector a, Vector b){
      return a.x * b.y - a.y * b.x;
    }

    int ccw(Vector p0, Vector p1, Vector p2){
      Vector a = p1 - p0;
      Vector b = p2 - p0;

      if(cross(a, b) > EPS){
        return COUNTER_CLOCKWISE;
      }
      if(cross(a, b) < -EPS){
        return CLOCKWISE;
      }
      if(dot(a, b) < -EPS){
        return ONLINE_BACK;
      }
      if(a.norm() < b.norm()){
        return ONLINE_FRONT;
      }

      return ON_SEGMENT;
    }

    double getCrossPointDistance(Vector p1, Vector p2, Vector p3, Vector p4){
      Vector base = p4 - p1;
      double d1 = abs(cross(base, p1 - p3));
      double d2 = abs(cross(base, p2 - p3));
      double t = d1 / (d1 + d2);
      Vector result = p1 + (p2 - p1) * t;
      double dy = p1.y - result.y;
      double dx = p1.x - result.x;
      double dy2 = p2.y - result.y;
      double dx2 = p2.x - result.x;

      return (dy*dy + dx*dx + dy2*dy2 + dx2*dx2);
    }

    double getDistanceLP(Vector p1, Vector p2, Vector p){
      return abs(cross(p2 - p1, p - p1) / abs(p2 - p1));
    }

    double getDistanceSP(Vector p1, Vector p2, Vector p){
      if(dot(p2 - p1, p - p1) < 0.0) return abs(p - p1);
      if(dot(p1 - p2, p - p2) < 0.0) return abs(p - p2);
      return getDistanceLP(p1, p2, p);
    }

    bool intersect2(Vector p1, Vector p2, Vector p3, Vector p4){
      return ((ccw(p1, p2, p3) * ccw(p1, p2, p4) <= 0) && (ccw(p3, p4, p1) * ccw(p3, p4, p2) < 0));
    }

    bool intersect(Vector *p1, Vector *p2, Vector *p3, Vector *p4){
      int d1 = direction(p3, p4, p1);
      int d2 = direction(p3, p4, p2);
      int d3 = direction(p1, p2, p3);
      int d4 = direction(p1, p2, p4);

      if(((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))){
        return true;
      }else{
        return false;
      }
    }

    /*
     * グラフの辺のみを作成する
     */
    void createEdge(){
			fprintf(stderr,"\ncreate edges =>\n");

      for(int i = 0; i < nodeCount; i++){
        Node *from = getNode(i);

        if(from->removed){
          continue;
        }

        for(int j = i+1; j < nodeCount; j++){
          Node *to = getNode(j);

          if(to->removed) continue;

          if(from->isNeighbor(to)){
						// 辺の追加
            Edge edgeA(from->id, to->id, to->area);
            Edge edgeB(to->id, from->id, from->area);

            edgeList.push_back(edgeA);
            edgeList.push_back(edgeB);

            from->addOriginNeighbor(to->id);
            to->addOriginNeighbor(from->id);
          }
        }
      }
    }

		/*
		 * グラフを作成する。削除済みのノードは使用しないように。
		 */
		void createGraph(){
			fprintf(stderr,"\ncreate graph =>\n");

      memset(pointUsedCount, 0, sizeof(pointUsedCount));

      createEdge();
      prim();
      //kruskal();
      /*
      for(int i = 0; i < nodeCount; i++){
        Node *nodeA = getNode(i);

				if(nodeA->removed) continue;

        for(int j = i+1; j < nodeCount; j++){
          Node *nodeB = getNode(j);

					if(nodeB->removed) continue;

          if(nodeA->isNeighbor(nodeB)){
            //fprintf(stderr,"node %d <-> node %d\n", nodeA->id, nodeB->id);
            nodeA->addNeighbor(nodeB->id);
            nodeB->addNeighbor(nodeA->id);

						// 辺の追加
            Edge edgeA(nodeA->id, nodeB->id, nodeB->area);
            Edge edgeB(nodeB->id, nodeA->id, nodeA->area);

            edgeList.push_back(edgeA);
            edgeList.push_back(edgeB);
          }
        }
      }
      */
		}

		/*
		 * グラフの掃除を行う
		 */
		int cleanGraph(){
			fprintf(stderr,"clean graph =>\n");
			int removeNodeCount = 0;
      priority_queue< Node, vector<Node>, greater<Node>  > pque;

      for(int i = 0; i < nodeCount; i++){
				Node node = nodeList[i];
				if(node.removed) continue;
        pque.push(node);
      }

      while(!pque.empty()){
        Node n = pque.top(); pque.pop();

        if(n.degree == 0){
          //fprintf(stderr,"Node %d removed!: degree = %d\n", n.id, n.degree);
					removeNodeCount += 1;
          removeNode(n.id);
        }else if(canNodeRemove(n.id)){
          fprintf(stderr,"Node %d remove!\n", n.id);
					removeNodeCount += 1;
          //fprintf(stderr,"p1 = %d, p2 = %d, p3 = %d\n", n.p1->id, n.p2->id, n.p3->id);
          removeNode(n.id);
        }
      }

			return removeNodeCount;
		}

		/*
		 * グラフのクリアを行う
		 */
		void resetGraph(){
			fprintf(stderr,"reset graph =>\n");
			edgeList.clear();

			for(int id = 0; id < nodeCount; id++){
				Node *node = getNode(id);

				node->clearNeighbor();
			}
		}

		/*
		 * uとvが接続しているかどうかを調べる
		 */
		bool isConnect(int u, int v){
			queue<int> que;
			map<int, bool> checkList;

			que.push(u);

			while(!que.empty()){
				int id = que.front(); que.pop();

				if(checkList[id]) continue;
				checkList[id] = true;

				if(id == v) return true;

				Node *node = getNode(id);

      	set<int>::iterator it = node->neighbors.begin();

      	while(it != node->neighbors.end()){
        	int nid = (*it);
        	Node *neighbor = getNode(nid);

					if(!neighbor->removed){
						que.push(nid);
					}

        	it++;
      	}
			}

			return false;
		}

    Node createNode(int nodeId, Triangle t){
      Node node;
      node.id   = nodeId;
      node.p1   = t.p1;
      node.p2   = t.p2;
      node.p3   = t.p3;
      node.y    = (t.p1->y + t.p2->y + t.p3->y) / 3.0;
      node.x    = (t.p1->x + t.p2->x + t.p3->x) / 3.0;
      node.removed = false;
      node.used = false;
      node.area = calcArea(t.p1, t.p2, t.p3);

      return node;
    }

    void removeNode(int nodeId){
      Node *node = getNode(nodeId);

      set<int>::iterator it = node->neighbors.begin();
      node->used = true;
      node->removed = true;

      while(it != node->neighbors.end()){
        int nid = (*it);
        Node *neighbor = getNode(nid);
        // 隣接ノードから自分を削除
        neighbor->removeNeighbor(nodeId);
        it++;
      }

      /*
      pointUsedCount[node->p1->id] = max(0, pointUsedCount[node->p1->id] - 1);
      pointUsedCount[node->p2->id] = max(0, pointUsedCount[node->p2->id] - 1);
      pointUsedCount[node->p3->id] = max(0, pointUsedCount[node->p3->id] - 1);
      */
    }

    /*
     * ノードが削除出来るかどうかを確認
     */
    bool canNodeRemove(int nodeId){
      Node *node = getNode(nodeId);

			vector<int> neighbors;

      if(node->degree <= 1){
        /*
        fprintf(stderr,"Node = %d, p%d -> %d, p%d -> %d, p%d -> %d\n", nodeId,
            node->p1->id, pointUsedCount[node->p1->id],
            node->p2->id, pointUsedCount[node->p2->id],
            node->p3->id, pointUsedCount[node->p3->id]);
            */
        if(pointUsedCount[node->p1->id] > 1 && pointUsedCount[node->p2->id] > 1 && pointUsedCount[node->p3->id] > 1){
          return true;
        }
      }else{
        return false;
      }

			// 次数が1のノードは消せない
			// if(node->degree == 1) return false;
      //
      set<int>::iterator it = node->originNeighbors.begin();

      // 隣接ノードに次数が1のノードがある場合は削除しない
      while(it != node->originNeighbors.end()){
        int nid = (*it);
        Node *neighbor = getNode(nid);

        if(neighbor->degree == 1) return false;
				neighbors.push_back(nid);

        it++;
      }

			if(node->degree == 2){
				node->removed = true;

				if(!isConnect(neighbors[0], neighbors[1])){
					node->removed = false;
					return false;
				}

				node->removed = false;
				return true;
			}else if(false && node->degree == 3){
				node->removed = true;

				if(!isConnect(neighbors[0], neighbors[1])){
					node->removed = false;
					return false;
				}

				if(!isConnect(neighbors[1], neighbors[2])){
					node->removed = false;
					return false;
				}

				return true;
			}

			return false;
    }

    vector<string> choosePolygons(vector<int> points, int n){
      vector<string> result;
      pointCount = points.size()/2;
      fprintf(stderr,"N = %d\n", n);

      init(points);
      set<Vector> vertices;
      set<Triangle> triangles;

      for(int id = 0; id < pointCount; id++){
        Vector v;
        v.id = id;
        v.y = points[id*2];
        v.x = points[id*2+1];

        vertices.insert(v);
      }

      Delaunay2d::getDelaunayTriangles(vertices, &triangles);

      fprintf(stderr,"triangle num = %lu\n", triangles.size());

      set<Triangle>::iterator it = triangles.begin();

      int nodeId = 0;
      nodeCount = 0;

      /*
       * 各三角形をノード化
       */
      while(it != triangles.end()){
        Triangle t = (*it);

        Node node = createNode(nodeId, t);
        nodeList[nodeId] = node;

        //fprintf(stderr,"Node %d - centerY = %4.2f, centerX = %4.2f, area = %4.2f\n", node.id, node.y, node.x, node.area);
        it++;
        nodeId += 1;
        nodeCount += 1;
      }

      /*
      createEdge();
      // クラスカルで最小全域木を作成
      kruskal();
      */

			int limit = 1;

      /*
       * グラフを整備する
       */
			for(int i = 0; i < limit; i++){
				createGraph();

				int removedNodeCount = cleanGraph();
				if(removedNodeCount == 0) break;

				if(i != limit-1){
					resetGraph();
				}
			}

      set<int> notUsePoints;
      for(int i = 0; i < pointCount; i++){
        if(pointUsedCount[i] == 0){
          fprintf(stderr,"Point %d is not used\n", i);
          notUsePoints.insert(i);
        }
      }

      for(int id = 0; id < nodeCount; id++){
        Node *node = getNode(id);

				if(node->removed || node->used) continue;

        if(!node->used && node->degree <= 1){
          result.push_back(createPolygon(id));
        }

        /*
        string str = "";
        str += int2string(node->p1->id);
        str += " ";
        str += int2string(node->p2->id);
        str += " ";
        str += int2string(node->p3->id);

        result.push_back(str);
        */
      }

      int cnt = notUsePoints.size();

      for(int i = 0; i < cnt; i++){
        addNewPoint(lines, notUsePoints);
        //edgeListClearSimple(lines);
      }

      result.clear();
      result.push_back(polygon2string(lines));

      return result;
    }

    void addNewPoint(vector<int> &lines, set<int> &notUsePoints){
      fprintf(stderr,"add new point =>\n");
      int lineCount = lines.size();
      priority_queue< Point, vector<Point>, greater<Point>  > pque;

      for(int i = 0; i < lineCount; i++){
        Vector *p1 = getVector(lines[i]);
        Vector *p2 = getVector(lines[(i+1)%lineCount]);

        set<int>::iterator it = notUsePoints.begin();

        while(it != notUsePoints.end()){
          Vector *p = getVector(*it);

          double dist = getDistanceSP(*p1, *p2, *p);

          if(dist < 150.0){
            pque.push(Point(p1, p2, p, dist));
          }

          it++;
        }
      }

      while(!pque.empty()){
        Point pt = pque.top(); pque.pop();
        vector<int> tempLines = lines;

        insertVertex(pt.p1, pt.p2, pt.p3, lines);

        if(lineCross(lines)){
          lines = tempLines;
        }else{
          notUsePoints.erase(pt.p3->id);
          break;
        }
      }
    }
};

int main(){
  int n, np, point;
  SmallPolygons sp;
  vector<int> ps;

  cin >> np;

  for(int i = 0; i < np; i++){
    cin >> point;
    ps.push_back(point);
  }

  cin >> n;

  vector<string> ret = sp.choosePolygons(ps, n);

  int size = ret.size();

  cout << size << endl;

  for(int i = 0; i < size; i++){
    string str = ret[i];
    cout << str << endl;
  }

  return 0;
}
