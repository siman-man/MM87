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
typedef pair<int, int> nextRootId;

const int MAX_NP            =  1600;  // 頂点の最大数
const int UNKNOWN           =    -1;  // 未定義
const int COUNTER_CLOCKWISE =     1;
const int CLOCKWISE         =    -1;
const int ONLINE_BACK       =     2;
const int ONLINE_FRONT      =    -2;
const int ON_SEGMENT        =     0;
const double EPS            = 1e-10;  // 誤差

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
  double x, y;

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

struct Neighbors{
  vector<int> nlist;

  int first(){ return nlist[0]; }
  int second(){ return nlist[1]; }
  int third(){ return nlist[2]; }
};

struct Node{
  int id;                       // ノードのID
  int degree;                   // ノードの次数
  int originDegree;             // 元の次元数
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
    this->id      = nodeId;
    this->used    = false;
    this->removed = false;
    this->degree  = 0;
  }

  /*
   * 隣接ノードの追加、次数が1増える
   */
  void addNeighbor(int nodeId){
    neighbors.insert(nodeId);
    degree = neighbors.size();
  }

  /*
   * 本来接続していたノード
   */
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

  /*
   * 隣接ノードを完全に削除する
   */
  void deleteNeighbor(int nodeId){
    originNeighbors.erase(nodeId);
    originDegree = originNeighbors.size();
  }

  /*
   * 隣接しているノードをリセット
   */
	void clearNeighbor(){
		neighbors.clear();
		degree = 0;
	}

  /*
   * 共有している点が2点あれば隣接している
   */
  bool isNeighbor(Node *node){
    int sameCnt = 0;
    sameCnt += (p1->id == node->p1->id || p1->id == node->p2->id || p1->id == node->p3->id);
    sameCnt += (p2->id == node->p1->id || p2->id == node->p2->id || p2->id == node->p3->id);
    sameCnt += (p3->id == node->p1->id || p3->id == node->p2->id || p3->id == node->p3->id);

    return sameCnt == 2;
  }

  /*
   * 共有している点が1点あれば接している
   */
  bool isAttach(Node *node){
    int sameCnt = 0;
    sameCnt += (p1->id == node->p1->id || p1->id == node->p2->id || p1->id == node->p3->id);
    sameCnt += (p2->id == node->p1->id || p2->id == node->p2->id || p2->id == node->p3->id);
    sameCnt += (p3->id == node->p1->id || p3->id == node->p2->id || p3->id == node->p3->id);

    return sameCnt > 0;
  }

  bool operator >(const Node &n) const{
    return area < n.area;
  }    
};

/*
 * 線分の情報
 */
struct Line{
  int u;  // 始点
  int v;  // 終点

  Line(int u = UNKNOWN, int v = UNKNOWN){
    this->u = u;
    this->v = v;
  }
};

/*
 * 各ノードを連結する辺情報
 */
struct Edge{
  int u;        // 始点
  int v;        // 終点
  double cost;  // コスト(終点側ノードの面積)

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

// 多角形を表す構造体
struct Polygon{
  int id;
  set<int> nodes;     // ノードの一覧
  double totalArea;   // 多角形の総面積

  Polygon(){
    totalArea = 0.0;
  }

  bool operator >(const Polygon &p) const{
    return totalArea < p.totalArea;
  }    
};

// ノード一覧
Node nodeList[10000];

typedef pair< bool, pair<Polygon, Polygon> > Polygons;

class SmallPolygons{
  public:
    /*
     * 初期化処理
     */
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

    /*
     * ノード情報の取得
     */
    inline Node *getNode(int id){
      return &nodeList[id];
    }

    /*
     * 頂点情報の取得
     */
    inline Vector *getVector(int id){
      return &vectorList[id];
    }
    
    /*
     * 隣接ノード一覧の更新
     */
    void updateNeighbor(Node *from){
      from->clearNeighbor();

      for(int id = 0; id < nodeCount; id++){
        Node *to = getNode(id);

        if(to->removed || from->id == id) continue;

        if(from->isNeighbor(to)){
          from->addOriginNeighbor(to->id);
          to->addOriginNeighbor(from->id);
          from->addNeighbor(to->id);
          to->addNeighbor(from->id);
        }
      }
    }

    /*
     * 頂点の使用されている数の更新を行う
     */
    void refreshPointUsedCount(){
      memset(pointUsedCount, 0, sizeof(pointUsedCount));

      for(int id = 0; id < nodeCount; id++){
        Node *node = getNode(id);

        if(node->removed) continue;

        pointUsedCount[node->p1->id] += 1;
        pointUsedCount[node->p2->id] += 1;
        pointUsedCount[node->p3->id] += 1;
      }
    }

    /*
     * プリム法を使って最小木を作成（出来るところまで）
     */
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
          //fprintf(stderr,"%d -> %d: No new point\n", from->id, to->id);
          continue;
        }

        //fprintf(stderr,"%d: %d <---> %d\n", cnt, e.u, e.v);

        if(from->removed) continue;

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

    /*
     * 三角形の面積を計算
     */
    double calcTriangleArea(const Vector *p1, const Vector *p2, const Vector *p3){
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

    /*
     * 頂点間の距離を事前に計算しておく
     */
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

    /*
     * 頂点の挿入を行う
     */
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

    /*
     * 三角形を交換
     */
    bool swapTriangle(int nodeID_A, int nodeID_B){
      fprintf(stderr,"Swap Triangle! %d <-> %d\n", nodeID_A, nodeID_B);
      Node *nodeA = getNode(nodeID_A);  
      Node *nodeB = getNode(nodeID_B);

      Line shareLine = getShareLine(nodeID_A, nodeID_B);
      Line crossLine = getCrossLine(nodeID_A, nodeID_B);

      Vector *p1 = getVector(shareLine.u);
      Vector *p2 = getVector(shareLine.v);
      Vector *p3 = getVector(crossLine.u);
      Vector *p4 = getVector(crossLine.v);

      double areaA = calcTriangleArea(p3, p4, p1);
      double areaB = calcTriangleArea(p3, p4, p2);

      if(max(areaA, areaB) > max(nodeA->area, nodeB->area)){
        cleanMe(nodeA);
        cleanMe(nodeB);

        nodeA->p1 = p3;
        nodeA->p2 = p4;
        nodeA->p3 = p1;
        nodeA->area = areaA;

        nodeB->p1 = p3;
        nodeB->p2 = p4;
        nodeB->p3 = p2;
        nodeB->area = areaB;

        updateNeighbor(nodeA);
        updateNeighbor(nodeB);
      }else{
        return false;
      }

      return true;
    }

    void cleanMe(Node *node){
      set<int>::iterator it = node->originNeighbors.begin();

      while(it != node->originNeighbors.end()){
        Node *neighbor = getNode(*it);

        neighbor->deleteNeighbor(node->id);

        it++;
      }

      node->originNeighbors.clear();
    }

    /*
     * 多角形の作成を行う
     */
    Polygon createPolygon(int nodeId){
			//fprintf(stderr,"create polygon = %d\n", nodeId);
      vector<int> vlist;
      queue<int> que;
      map<int, bool> checkList;

      Polygon polygon;

      que.push(nodeId);

      while(!que.empty()){
        int id = que.front(); que.pop();

        if(checkList[id]) continue;
        checkList[id] = true;
        // ノードのリストのに追加
        polygon.nodes.insert(id);

        Node *node = getNode(id);
        node->used = true;

        // 合計面積に追加
        polygon.totalArea += node->area;

        set<int>::iterator nid = node->neighbors.begin();

        while(nid != node->neighbors.end()){
          que.push(*nid);
          nid++;
        }
      }

      return polygon;
    }

    /*
     * 与えられた２つのノードが削除出来るかどうかを確認
     */
    bool checkDivide(int leftId, int rightId){
      Node *one, *two;
      Line line = getShareLine(rightId, leftId);

      // 点が2点しか共有していない場合は、削除すると独立した点になってしまうので分割しない
      if(pointUsedCount[line.u] <= 2 || pointUsedCount[line.v] <= 2) return true;

      Neighbors rightNeighbors = createNeighbors(rightId);
      Neighbors leftNeighbors = createNeighbors(leftId);

      if(rightNeighbors.first() == leftId){
        one = getNode(rightNeighbors.second());

        if(leftNeighbors.first() == rightId){
          two = getNode(leftNeighbors.second());
        }else{
          two = getNode(leftNeighbors.first());
        }
      }else{
        one = getNode(rightNeighbors.first());

        if(leftNeighbors.first() == rightId){
          two = getNode(leftNeighbors.second());
        }else{
          two = getNode(leftNeighbors.first());
        }
      }

      return one->isAttach(two);
    }

    /*
     * 分割出来るノードかどうかを調べる
     */
    bool canDivide(Node *node){
      if(node->degree != 2) return false;

      set<int>::iterator id = node->neighbors.begin();
      int cnt = 0;

      while(id != node->neighbors.end()){
        Node *neighbor = getNode(*id);

        cnt += (neighbor->degree == 2 && !checkDivide(node->id, neighbor->id));

        id++;
      }

      return (cnt > 0);
    }

    /*
     * ノードを削除して次のポリゴンの要素の一部を返す
     *   removeID_A: 削除するノードのID
     *   removeID_B: 削除するノードのID
     */
    Polygon divideNode(int removeID_A, int removeID_B){
      Neighbors neighbors = createNeighbors(removeID_A);
      Polygon polygon;

      removeNode(removeID_A);

      if(neighbors.first() != removeID_B){
        polygon = createPolygon(neighbors.first());
      }else{
        polygon = createPolygon(neighbors.second());
      }

      return polygon;
    }

    /*
     * 多角形を二等分する
     */
    Polygons dividePolygon(Polygon polygon){
      //fprintf(stderr,"dividePolygon =>\n");
      priority_queue< Node, vector<Node>, greater<Node> > pque;
      Polygons polygons;

      set<int>::iterator id = polygon.nodes.begin();

      // 自身と隣接しているノードの次数が2のノードだけを入れる
      while(id != polygon.nodes.end()){
        Node *node = getNode(*id);

        if(canDivide(node)){
          pque.push(*node);     
        }
        id++;
      }

      if(pque.size() == 0){
        //fprintf(stderr,"Failed divide polygon...\n");
        polygons.first = false;
        return polygons;
      }else{
        //fprintf(stderr,"Success divide polygon!\n");
        int removeID_A = pque.top().id;
        polygons.first = true;

        Neighbors neighbors = createNeighbors(removeID_A);
        Node *right = getNode(neighbors.first());
        Node *left  = getNode(neighbors.second());

        if(right->degree == 2 && !checkDivide(removeID_A, right->id)){
          polygons.second.first  = divideNode(removeID_A, right->id);
          polygons.second.second = divideNode(right->id, removeID_A);
        }else if(left->degree == 2 && !checkDivide(removeID_A, left->id)){
          polygons.second.first  = divideNode(removeID_A, left->id);
          polygons.second.second = divideNode(left->id, removeID_A);
        }
      }

      return polygons;
    }

    vector<int> polygon2vlist(Polygon polygon){
			//fprintf(stderr,"polygon2vlist =>\n");
      vector<int> vlist;
      queue<int> que;
      map<int, bool> checkList;
      int p1index, p2index, p3index;
      int listSize;

      que.push(*polygon.nodes.begin());

      while(!que.empty()){
        int id = que.front(); que.pop();

        if(checkList[id]) continue;
        checkList[id] = true;

        Node *node = getNode(id);

        listSize = vlist.size();

        if(listSize == 0){
          vlist.push_back(node->p1->id);
          vlist.push_back(node->p2->id);
          vlist.push_back(node->p3->id);
        }else{
          p1index = find(vlist.begin(), vlist.end(), node->p1->id) - vlist.begin();
          p2index = find(vlist.begin(), vlist.end(), node->p2->id) - vlist.begin();
          p3index = find(vlist.begin(), vlist.end(), node->p3->id) - vlist.begin();

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
          if(!checkList[nid]) que.push(nid);
          it++;
        }
      }

      return vlist;
    }

    /*
     * 解答用に頂点のリストを文字列に変換
     */
    string vlist2string(vector<int> &lines){
      int listSize = lines.size();
      string result = "";

      for(int i = 0; i < listSize; i++){
        result += int2string(lines[i]);

        if(i != listSize-1) result += " ";
      }

      return result;
    }

    /*
     * 多角形を構成している、線分が交差していないかどうかを調べる
     */
    bool lineCross(vector<int> &vlist){
      int listSize = vlist.size();

      for(int i = 0; i < listSize; i++){
        Vector *p1 = &vectorList[vlist[i%listSize]];
        Vector *p2 = &vectorList[vlist[(i+1)%listSize]];

        for(int j = i+1; j < i+listSize; j++){
          Vector *p3 = &vectorList[vlist[j%listSize]];
          Vector *p4 = &vectorList[vlist[(j+1)%listSize]];

          if(intersect(*p1, *p2, *p3, *p4)) return true;
        }
      }

      return false;
    }

    int direction(Vector *p0, Vector *p1, Vector *p2){
      return (p1->x - p0->x) * (p2->y - p0->y) - (p2->x - p0->x) * (p1->y - p0->y);
    }

    // 線分上に存在しているかどうかを判定
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

      if(cross(a, b) > EPS) return COUNTER_CLOCKWISE;
      if(cross(a, b) < -EPS) return CLOCKWISE;
      if(dot(a, b) < -EPS) return ONLINE_BACK;
      if(a.norm() < b.norm()) return ONLINE_FRONT;

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

    bool intersect(Vector p1, Vector p2, Vector p3, Vector p4){
      return ((ccw(p1, p2, p3) * ccw(p1, p2, p4) <= 0) && (ccw(p3, p4, p1) * ccw(p3, p4, p2) < 0));
    }

    /*
     * グラフの辺のみを作成する
     */
    void createEdge(){
			fprintf(stderr,"\ncreate edges =>\n");

      for(int i = 0; i < nodeCount; i++){
        Node *from = getNode(i);
        for(int j = i+1; j < nodeCount; j++){
          Node *to = getNode(j);

          if(from->isNeighbor(to)){
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
      //cleanTriangles();
      prim();
		}

    void cleanTriangles(){
      for(int id = 0; id < nodeCount; id++){
        Node *node = getNode(id);

        if(node->originDegree == 2){
          swapTriangle(node->id, *node->originNeighbors.begin());
        }
      }
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

        // 次数0のノードは削除
        if(n.degree == 0){
					removeNodeCount += 1;
          removeNode(n.id);
        }else if(canNodeRemove(n.id)){
          fprintf(stderr,"Node %d remove!\n", n.id);
					removeNodeCount += 1;
          removeNode(n.id);
        }
      }

			return removeNodeCount;
		}

		/*
		 * グラフのリセットを行う
		 */
		void resetGraph(){
			fprintf(stderr,"reset graph =>\n");

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

    /*
     * ノードAとノードBが共有している線分を取り出す
     */
    Line getShareLine(int nodeID_A, int nodeID_B){
      Node *nodeA = getNode(nodeID_A);
      Node *nodeB = getNode(nodeID_B);

      int p1A = nodeA->p1->id;
      int p2A = nodeA->p2->id;
      int p3A = nodeA->p3->id;

      int p1B = nodeB->p1->id;
      int p2B = nodeB->p2->id;
      int p3B = nodeB->p3->id;

      if(p1A != p1B && p1A != p2B && p1A != p3B){
        return Line(p2A, p3A);
      }
      if(p2A != p1B && p2A != p2B && p2A != p3B){
        return Line(p1A, p3A);
      }

      return Line(p1A, p2B);
    }

    /*
     * 分割する線分を取得する
     */
    Line getCrossLine(int nodeID_A, int nodeID_B){
      Node *nodeA = getNode(nodeID_A);
      Node *nodeB = getNode(nodeID_B);

      map<int, int> pList;

      pList[nodeA->p1->id] += 1;
      pList[nodeA->p2->id] += 1;
      pList[nodeA->p3->id] += 1;

      pList[nodeB->p1->id] += 1;
      pList[nodeB->p2->id] += 1;
      pList[nodeB->p3->id] += 1;

      map<int, int>::iterator it = pList.begin();
      Line line;

      while(it != pList.end()){
        int cnt = (*it).second;

        if(cnt == 1){
          if(line.u == UNKNOWN){
            line.u = (*it).first;
          }else{
            line.v = (*it).first;
          }
        }

        it++;
      }

      return line;
    }


    Neighbors createNeighbors(int nodeId){
      Node *node = getNode(nodeId);
      Neighbors result;

      set<int>::iterator it = node->neighbors.begin();

      while(it != node->neighbors.end()){
        result.nlist.push_back(*it);
        it++; 
      }

      return result;
    }

    /*
     * ノードを作成する
     */
    Node createNode(int nodeId, Triangle t){
      Node node;
      node.id      = nodeId;
      node.p1      = t.p1;
      node.p2      = t.p2;
      node.p3      = t.p3;
      node.y       = (t.p1->y + t.p2->y + t.p3->y) / 3.0;
      node.x       = (t.p1->x + t.p2->x + t.p3->x) / 3.0;
      node.removed = false;
      node.used    = false;
      node.area    = calcTriangleArea(t.p1, t.p2, t.p3);

      return node;
    }

    /*
     * ノードの追加を行う
     */
    void addNode(Node node){
      nodeList[nodeCount] = node;

      nodeCount += 1;
    }

    /*
     * ノードの削除を行う
     *   nodeId: ノードID
     */
    void removeNode(int nodeId){
      Node *node = getNode(nodeId);

      set<int>::iterator it = node->neighbors.begin();
      node->removed = true;

      // 隣接ノードから自分を削除
      while(it != node->neighbors.end()){
        int nid = (*it);
        Node *neighbor = getNode(nid);
        neighbor->removeNeighbor(nodeId);
        it++;
      }

      // 自分の隣接ノードを削除
      node->neighbors.clear();
      node->degree = 0;

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

      if(node->degree <= 1){
        if(pointUsedCount[node->p1->id] > 1 && pointUsedCount[node->p2->id] > 1 && pointUsedCount[node->p3->id] > 1){
          return true;
        }
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

      // ドロネー三角分割
      Delaunay2d::getDelaunayTriangles(vertices, &triangles);

      fprintf(stderr,"triangle num = %lu\n", triangles.size());

      set<Triangle>::iterator it = triangles.begin();

      nodeCount = 0;

      /*
       * 各三角形をノード化
       */
      while(it != triangles.end()){
        Triangle t = (*it);

        Node node = createNode(nodeCount, t);
        addNode(node);

        it++;
      }

			createGraph();
			cleanGraph();

      set<int> notUsePoints;
      for(int i = 0; i < pointCount; i++){
        if(pointUsedCount[i] == 0){
          //fprintf(stderr,"Point %d is not used\n", i);
          notUsePoints.insert(i);
        }
      }

      Polygon rootPolygon;

      /*
       * 多角形を構成するノードの一覧を取得
       */
      for(int id = 0; id < nodeCount; id++){
        Node *node = getNode(id);

				if(node->removed || node->used) continue;
        rootPolygon = createPolygon(id);
      }


      vector<int> lines = polygon2vlist(rootPolygon);
      int cnt = notUsePoints.size();

      for(int i = 0; i < cnt; i++){
        addNewPoint(lines, notUsePoints, rootPolygon);
      }

      result.clear();

      /*
      lines = polygon2vlist(rootPolygon);
      string str = vlist2string(lines);
      result.push_back(str);
      */

      priority_queue< Polygon, vector<Polygon>, greater<Polygon> > pque;
      pque.push(rootPolygon);

      queue<Polygon> finalPolygons;
      int dcnt = 0;
      int dlimit = n-1;

      while(!pque.empty()){
        Polygon poly = pque.top(); pque.pop();

        if(dcnt < dlimit){
          Polygons polygons = dividePolygon(poly);

          if(polygons.first){
            pque.push(polygons.second.first);
            pque.push(polygons.second.second);

            dcnt++;
            refreshPointUsedCount();
          }else{
            finalPolygons.push(poly);
          }
        }else{
          finalPolygons.push(poly);
        }
      }

      fprintf(stderr,"finalPolygons size = %lu\n", finalPolygons.size());

      while(!finalPolygons.empty()){
        Polygon poly = finalPolygons.front(); finalPolygons.pop();

        vector<int> lines = polygon2vlist(poly);
        string str = vlist2string(lines);
        result.push_back(str);
        
        /*
        set<int>::iterator that = poly.nodes.begin();

        while(that != poly.nodes.end()){
          Node *node = getNode(*that);
          string str = "";
          str += int2string(node->p1->id);
          str += " ";
          str += int2string(node->p2->id);
          str += " ";
          str += int2string(node->p3->id);

          result.push_back(str);
          that++;
        }

        */
      }

      return result;
    }

    void addNewPoint(vector<int> &lines, set<int> &notUsePoints, Polygon &rootPolygon){
      //fprintf(stderr,"add new point =>\n");
      int lineCount = lines.size();
      priority_queue< Point, vector<Point>, greater<Point>  > pque;

      for(int i = 0; i < lineCount; i++){
        Vector *p1 = getVector(lines[i]);
        Vector *p2 = getVector(lines[(i+1)%lineCount]);

        set<int>::iterator it = notUsePoints.begin();

        while(it != notUsePoints.end()){
          Vector *p = getVector(*it);

          double dist = getDistanceSP(*p1, *p2, *p);

          if(dist < 200.0){
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
          Triangle t;
          t.p1 = pt.p1;
          t.p2 = pt.p2;
          t.p3 = pt.p3;
          Node node = createNode(nodeCount, t);
          addNode(node);
          updateNeighbor(&node);
          rootPolygon.nodes.insert(node.id);

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
  for(int i = 0; i < np; i++){ cin >> point; ps.push_back(point);}
  cin >> n;
  vector<string> ret = sp.choosePolygons(ps, n);
  int size = ret.size();
  cout << size << endl;
  for(int i = 0; i < size; i++){string str = ret[i]; cout << str << endl;}
  return 0;
}
