# MM87
SmallPolygons

# Problem Statement
You are given a set of points with integer coordinates on a two-dimensional plane and an integer N. 
あなたはNの数の二次元座標上の点が与えられます
You have to construct at most N polygons which have these points as vertices. More specifically,
あなたはNの多角形を構築して下さい。より正確には

- Each vertex of each polygon must be a point from the given set.
多角形の各頂点は与えられた点の集合である必要があります。

- Each point from the given set must belong to exactly one of the constructed polygons.
与えられた各点は多角形の構成要素の1つである必要があります

- Each polygon must be simple.
各多角形はシンプルである必要があります。

- Edges of different polygons can't intersect (but one polygon can lie completely within another).
各多角形の辺は他の多角形の辺と交差してはいけません

Your task is to minimize the sum of areas of polygons that you constructed.
あなたの目的は、構成する多角形の合計面積を最小化することです。

Implementation Details
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

Scoring
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