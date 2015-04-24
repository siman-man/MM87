#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
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

const int MAX_NP = 5000;
const int UNKNOWN = -1;

int np;
int pointCount;
int pointY[MAX_NP];
int pointX[MAX_NP];
int pointsDistance[MAX_NP][MAX_NP];
int nodeCount;

string int2string(int number){
  stringstream ss; 
  ss << number;
  return ss.str();
}

struct Vector{
  int id;
  double x;
  double y;

  bool operator==(const Vector& v) const{
    return (x == v.x && y == v.y);
  }

  bool operator<(const Vector& v) const{
    return (x != v.x)? x < v.x : y < v.y;
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

        Vector *p1 = new Vector;
        p1->x = centerX - sqrt(3.0) * radius;
        p1->y = centerY - radius;

        Vector *p2 = new Vector;
        p2->x = centerX + sqrt(3.0) * radius;
        p2->y = centerY - radius;

        Vector *p3 = new Vector;
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
  int triangleId;               // このノードに対応している三角形
  set<int> neighbors;           // 隣接ノード一覧
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

class SmallPolygons{
  public:
    void init(vector<int> &points){

      for(int id = 0; id < pointCount; id++){
        pointY[id] = points[id*2];
        pointX[id] = points[id*2+1];
      }

      initializePointDistance();
    };

    Node *getNode(int id){
      return &nodeList[id];
    }

    /*
     * 現在各ノードがどのグループに属しているかを調べる
     */
    void divideCheck(){
    }

    void kruskal(){
      fprintf(stderr,"kruskal =>\n");
      sort(edgeList.begin(), edgeList.end(), comp);
      initUnionFind(np);

      int E = edgeList.size();

      for(int i = 0; i < E; i++){
        Edge e = edgeList[i];

        Node *nodeB = getNode(e.v);

        set<int>::iterator it = nodeB->neighbors.begin();
        bool flag = true;

        while(it != nodeB->neighbors.end()){
          int nid = (*it);
          
          if(nid != e.u && unionSame(e.u, nid)){
            flag = false;
          }

          it++;
        }

        if(flag && !unionSame(e.u, e.v)){
          Node *nodeA = getNode(e.u);

          if(nodeA->removed || nodeB->removed) continue;

          unionUnite(e.u, e.v);

          nodeA->addNeighbor(nodeB->id);
          nodeB->addNeighbor(nodeA->id);
          
          //fprintf(stderr,"node %d <-> node %d\n", e.u, e.v);
        }
      }
    }

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

          int dist = round(sqrt(dy*dy + dx*dx));

          pointsDistance[i][j] = dist;
          pointsDistance[j][i] = dist;
        }
      }
    }

    string createPolygon(int nodeId){
			fprintf(stderr,"create polygon = %d\n", nodeId);
      vector<int> list;
      queue<int> que;
      map<int, bool> checkList;

      que.push(nodeId);
      int p1index;
      int p2index;
      int p3index;
      int listSize;
      int right, left;

      while(!que.empty()){
        int id = que.front(); que.pop();

        if(checkList[id]) continue;
        checkList[id] = true;

				fprintf(stderr,"%d -> ", id);

        Node *node = getNode(id);
        node->used = true;

        listSize = list.size();
        vector<int>::iterator first = list.begin();
        fprintf(stderr,"nodeID = %d, vertex -> %d - %d - %d\n", node->id, node->p1->id, node->p2->id, node->p3->id);

        if(listSize == 0){
          if(node->p1->id < node->p2->id){
            fprintf(stderr,"add %d - %d - %d\n", node->p1->id, node->p2->id, node->p3->id);
            list.push_back(node->p1->id);
            list.push_back(node->p2->id);
            list.push_back(node->p3->id);
          }else{
            fprintf(stderr,"add %d - %d - %d\n", node->p2->id, node->p1->id, node->p3->id);
            list.push_back(node->p2->id);
            list.push_back(node->p1->id);
            list.push_back(node->p3->id);
          }
        }else{
          p1index = find(list.begin(), list.end(), node->p1->id) - list.begin();
          p2index = find(list.begin(), list.end(), node->p2->id) - list.begin();
          p3index = find(list.begin(), list.end(), node->p3->id) - list.begin();

          fprintf(stderr,"p1index = %d, p2index = %d, p3index = %d\n", p1index, p2index, p3index);

          if(p1index >= listSize){
            right = max(p2index, p3index);
            left  = min(p2index, p3index);
            int idx = first + right - first;

            if(left == 0 && right == listSize-1){
              fprintf(stderr,"add %d to %d\n", node->p1->id, idx+1);
              list.insert(first + right + 1, node->p1->id);
            }else{
              fprintf(stderr,"add %d to %d\n", node->p1->id, idx);
              list.insert(first + right, node->p1->id);
            }
          }else if(p2index >= listSize){
            right = max(p1index, p3index);
            left  = min(p1index, p3index);
            int idx = first + right - first;

            if(left == 0 && right == listSize-1){
              fprintf(stderr,"add %d to %d\n", node->p2->id, idx+1);
              list.insert(first + right + 1, node->p2->id);
            }else{
              fprintf(stderr,"add %d to %d\n", node->p2->id, idx);
              list.insert(first + right, node->p2->id);
            }
          }else if(p3index >= listSize){
            right = max(p1index, p2index);
            left  = min(p1index, p2index);
            int idx = first + right - first;

            if(left == 0 && right == listSize-1){
              fprintf(stderr,"add %d to %d\n", node->p3->id, idx);
              list.insert(first + right + 1, node->p3->id);
            }else{
              fprintf(stderr,"add %d to %d\n", node->p3->id, idx);
              list.insert(first + right, node->p3->id);
            }
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

			fprintf(stderr,"\n");

      string result = "";
      listSize = list.size();

      for(int i = 0; i < list.size(); i++){
        result += int2string(list[i]);

        if(i != listSize-1){
          result += " ";
        }
      }

      return result;
    }

		/*
		 * グラフを作成する。削除済みのノードは使用しないように。
		 */
		void createGraph(){
			fprintf(stderr,"\ncreate graph =>\n");

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

        if(canNodeRemove(n.id)){
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
    }

    /*
     * ノードが削除出来るかどうかを確認
     */
    bool canNodeRemove(int nodeId){
      Node *node = getNode(nodeId);

      set<int>::iterator it = node->neighbors.begin();
			vector<int> neighbors;

			// 次数が1のノードは消せない
			if(node->degree == 1) return false;

      // 隣接ノードに次数が1のノードがある場合は削除しない
      while(it != node->neighbors.end()){
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
			}else if(node->degree == 3){
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

      while(it != triangles.end()){
        Triangle t = (*it);

        Node node = createNode(nodeId, t);
        nodeList[nodeId] = node;

        //fprintf(stderr,"Node %d - centerY = %4.2f, centerX = %4.2f, area = %4.2f\n", node.id, node.y, node.x, node.area);
        it++;
        nodeId += 1;
        nodeCount += 1;
      }

			int limit = 10;

			for(int i = 0; i < limit; i++){
				createGraph();
				int removedNodeCount = cleanGraph();

				if(removedNodeCount == 0) break;

				if(i != limit-1){
					resetGraph();
				}
			}

      fprintf(stderr,"\nConnect!\n\n");

      for(int id = 0; id < nodeCount; id++){
        Node *node = getNode(id);

        if(node->removed){
          continue;
        }

        set<int>::iterator it = node->neighbors.begin();

        //fprintf(stderr,"Node %d - degree = %d -> ", node->id, node->degree);
        while(it != node->neighbors.end()){
          it++;
        }
      }

      kruskal();

      for(int id = 0; id < nodeCount; id++){
        Node *node = getNode(id);

				if(node->removed) continue;

        if(!node->used && node->degree <= 1){
          result.push_back(createPolygon(id));
        }
      }

      return result;
    }
};

int main(){
  int n, point;
  SmallPolygons sp;
  vector<int> points;

  cin >> np;

  for(int i = 0; i < np; i++){
    cin >> point;
    points.push_back(point);
  }

  cin >> n;

  vector<string> ret = sp.choosePolygons(points, n);

  int size = ret.size();

  cout << size << endl;

  for(int i = 0; i < size; i++){
    string str = ret[i];
    fprintf(stderr,"%s\n",str.c_str());
    cout << str << endl;
  }

  return 0;
}
