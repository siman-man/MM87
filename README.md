# MM87
SmallPolygons

# Problem Statement
You are given a set of points with integer coordinates on a two-dimensional plane and an integer N. 
あなたは二次元座標上の点と数値Nが与えられます

You have to construct at most N polygons which have these points as vertices. More specifically,
あなたは多くてもN個の多角形を構築して下さい。より正確には

- Each vertex of each polygon must be a point from the given set.
多角形の各頂点は与えられた点の集合である必要があります。

- Each point from the given set must belong to exactly one of the constructed polygons.
与えられた各点は多角形の構成要素の1つである必要があります(全部の点を最低1回は使ってください)

- Each polygon must be simple.
各多角形はシンプルである必要があります。

- Edges of different polygons can't intersect (but one polygon can lie completely within another).
各多角形の辺は他の多角形の辺と交差してはいけません

Your task is to minimize the sum of areas of polygons that you constructed.
あなたの目的は、構成する多角形の合計面積を最小化することです。

# Implementation Details
実装詳細

Your code should implement one method choosePolygons(vector <int> points, int N). 
あなたはchoosePolygons(vector<int> points, int N)を仮引数としてメソッドを作成します。

points gives you the set of points: point i (0-based) has coordinates (points[2*i], points[2*i+1]). 
各座標は(points[2*i], points[2*i+1])で与えられます

N is the number of polygons. You must return a set of polygons you've constructed from these points as a vector<string>
Nは多角形の数です。あなたは構築した多角形を返す必要があります

Each element of your return must describe one polygon and should be formatted as a space-separated list of 
あなたの返す各要素は、それぞれ多角形の構成要素を表しており、点が時計回りか、反時計周りとなっています。

its vertices in clockwise or counterclockwise order. Each vertex is given as its 0-based index in points.
拡張点は0スタートです

# Scoring
スコア計算

Your score for an individual test case will be the sum of areas of polygons you've constructed. 
書くテストにおけるあなたのスコアは構成した多角形の面積の合計値となります。

If your return has invalid format or specifies any invalid polygons, your score for the test case will be 0.
もしあなたが不正な多角形を構成した場合は、あなたのスコアは0点となります。

Your overall score will be calculated in the following way
あなたのスコアは次のように計算されます

for each test case where your score is not 0, you get 1 point for each competitor you beat on this test case 
全ての0点ではないテストケースにおいて、あなたのスコアが一番良い場合はあなたは1点を得ます

 (i.e., your score on a test case is smaller than this competitor's score)
 言い方を変えると、あなたのスコアは他の参加者よりも低いスコアです
 and 0.5 points for each competitor you tie with (a tie with yourself is not counted); 

 また、0.5ポイントは
 finally, the sum of points is divided by (the number of competitors - 1).
 そして最後に、ポイントの合計を参加者の人数で割ります。

 Test Case Generation
 テストケース作成方法

Each test case is generated as follows:
各テストケースは以下の手順で作成されます。

The number of points in the set is chosen either between 20 and 99, between 100 and 499, or between 500 and 1500, 
点の数はそれぞれ20 - 99, 100 - 499, 500 - 1500の中から選ばれます

all inclusive (either interval is chosen with equal probability).
(これらは等しい確率で選択されます)

The points coordinates are chosen between 0 and 699, inclusive, so that all points are distinct.
店の頂点は0 - 699の間から選択されます。 拡張点はそれぞれ独立しています。

The number of polygons allowed is chosen between 2 and 20, inclusive.
多角形の数は2 - 20の間で選択可能です。

All values are chosen uniformly and independently, at random.
それぞれの値は、一定になるようにランダムで選択されます。

Note that example 0 (seed 1) does not conform to these constraints. 
seed1に関してはこの制約にとらわれません。


# 戦略（次に実装するもの)

- [ ] 線分を表す構造体を作る
- [ ] 与えられた頂点に対して凸包を作る
  - [ ] 凸包は辺のリストで構成される(連結リスト)
- [ ] 凸包の中に含まれている1点を選び出して、凸包を構成している辺との三角形を作って、他の点を含まないように大きな三角形を作れる辺と接続する
  - [ ] 点を選ぶ(重心から遠いほうが良い？)
  - [ ] 選んだ点と多角形を構成している辺の1を取り出し三角形を構成する
    - [ ] 新たに出来た線分が既存の線分とかぶっていないか(多角形を構成している線分と比較すれば良いO(n))
    - [ ] 作成した三角形が、現在多角形に含まれていない点を内包しないかどうかを調べる
  - [ ] 新たに辺を作成し追加する
- [ ] 全部つなぎ終わったら1つの多角形が完成する
  - [ ] 構造としては辺のリストを持つ
- [ ] 構成した多角形から２つの辺を選び出しそれぞれ切断してつなぐ
  - [ ] A<->B C<->Dの２つを切断 -> A<->D B<->Cと繋ぎ直す(A<->D or B<->Cが他の線分と交差する場合は失敗)


# 戦略実装に必要な関数(既存の関数を書き換える)

list<Line> createConvexHull();
bool pointInclude(Vector *p1, Vector *p2, Vector *p3, Vector *p);
void addNewPoint(Vector *p);