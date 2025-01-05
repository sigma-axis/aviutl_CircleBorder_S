# AviUtl 円形縁取り拡張編集フィルタプラグイン

円形で縁取りするフィルタ効果などを追加する拡張編集フィルタプラグインです．

[ダウンロードはこちら．](https://github.com/sigma-axis/aviutl_CircleBorder_S/releases) [紹介動画．](https://www.nicovideo.jp/watch/sm43972213)

![使用例](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/0bc4f398-3566-413f-ae41-f17b1f3aa123)


次のフィルタ効果が追加されます．

1.  **縁取りσ**

    円形で縁取りをします．拡張編集標準のものは正方形なので角度によって線の幅にムラができやすいですが，こちらは角度によらず線の幅がほぼ均一になります．内側縁取りもできます．[[詳細](#縁取りσ)]

    ![縁取りσ](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/097b2def-a3dc-45ee-9103-0e1f9ca1888b)

    - [ティム様のスクリプト](https://tim3.web.fc2.com/sidx.htm)の「縁取りT」より高速に動作します．

      - 設定によっては 20 倍程度の差が出ます．
      - 「距離グラデーション」に相当する機能はありません．


1.  **角丸めσ**

    図形の角を丸めます．丸まった角は円弧の形に切り取られます．またオブジェクトの縁部分を透明ピクセルに削ることもできます．[[詳細](#角丸めσ)]

    ![角丸めσ](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/20c148d9-49e6-4d6b-9a0b-282cbe7e5649)


1.  **アウトラインσ**

    図形の境界に沿ったラインになります．ラインの幅や凹凸部分の曲率もコントロールできます．[[詳細](#アウトラインσ)]

    ![アウトラインσ](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/ae140f91-d4a3-482b-a97b-8d15d1036a0b)


## 動作要件

- AviUtl 1.10 + 拡張編集 0.92

  http://spring-fragrance.mints.ne.jp/aviutl
  - 拡張編集 0.93rc1 等の他バージョンでは動作しません．

- Visual C++ 再頒布可能パッケージ（\[2015/2017/2019/2022\] の x86 対応版が必要）

  https://learn.microsoft.com/ja-jp/cpp/windows/latest-supported-vc-redist

- patch.aul の `r43 謎さうなフォーク版58` (`r43_ss_58`) 以降

  https://github.com/nazonoSAUNA/patch.aul/releases/latest


## 導入方法

`aviutl.exe` と同階層にある `plugins` フォルダ内に `CircleBorder_S.eef` ファイルをコピーしてください．


## 使い方

各フィルタ効果には[方式](#方式-について)というパラメタが共通してついています．これは円形縁取りの計算アルゴリズムを指定するもので，場面場面によって得手不得手があったり計算速度に違いがあるなど特徴があります．初期値は[2値化倍精度](#2値化倍精度)になっていますが，場面に応じて適したものを選ぶと綺麗な縁取りが実現できるようになります．詳細は[こちら](#方式-について)．

### 縁取りσ

拡張編集の `縁取り` に似ていますが，円をベースにしているため角度によらず線の幅にムラがなく，回転しても違和感のない描画になります．内側への縁取りもできます．

![縁取りσ](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/097b2def-a3dc-45ee-9103-0e1f9ca1888b)

#### 各種パラメタ

![縁取りσの GUI](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/39729641-0c52-416f-b612-945dd974e354)


- サイズ

  縁取りのサイズをピクセル単位で指定します．正の値だと外側に縁取り（通常の縁取りと同等），負の値だと内側に縁取りします．`0.0` だと縁取り効果を無効化します．

  最小値は `-500.0`, 最大値は `500.0`, 初期値は `5.0`.

- 凹半径

  図形の凹んだ部分に丸みを持たせます．

  ![凹半径の適用例](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/c4552043-82a8-49d1-a4ad-d1bce896fbcc)

  - `サイズ` が負の場合は，図形の凸部分の内側縁取りに丸みを持たせます．

  - 手順的には例えば `サイズ` が `5.0` ピクセル，`凹半径` が `10.0` ピクセルだった場合，`15.0` ピクセル外側に縁取りをしてから `10.0` ピクセル内側に縁取りをすることで実現しています．

    - 2回分の計算が必要なため処理速度もその分遅くなります．

    - 「穴」のあるオブジェクトに対して適用した場合，一定の `凹半径` の大さ以上で急に「穴」が塞がるような不連続な変化が起きることがあるので注意してください．

  最小値は `0.0`, 最大値は `500.0`, 初期値は `0.0`.

- 透明度

  縁部分の透明度を % 単位で指定します．

  最小値は `0.0`, 最大値は `100.0`, 初期値は `0.0`.

- 内透明度

  縁取りの元となった図形の透明度を % 単位で指定します．`100.0` を指定すると，縁部分だけが表示されて「穴抜き」のような形状になります．

  最小値は `0.0`, 最大値は `100.0`, 初期値は `0.0`.

- ぼかし

  縁部分のぼかし幅を，`サイズ` との比で % 単位で指定します．

  最小値は `0.0`, 最大値は `100.0`, 初期値は `0.0`.

- αしきい値 / 基準α和

  `方式` で必要なパラメタを指定します．詳しくは[こちら](#方式-について)．

  最小値は `0.0`, 最大値は `100.0`, 初期値は `50.0`.

- 画像X / 画像Y

  `パターン画像ファイル` を設定している場合のみ有効，パターン画像ファイルの位置を調節します．

  画像ファイルの左上座標が，元のオブジェクトの左上座標と一致する配置が $(0,0)$ の基準点になります．

  最小値は `-4000`, 最大値は `4000`, 初期値は `0`.

- 方式

  円形縁取りのアルゴリズムを指定します．詳しくは[こちら](#方式-について)．

  初期値は [`2値化倍精度`](#2値化倍精度).

- 縁色の設定

  縁取りの色を指定します．パターン画像ファイルが設定されている場合は無視されます．

  初期値は `RGB( 0 , 0 , 0 )` （黒）です．

- パターン画像ファイル

  縁取り部分に適用するパターン画像ファイルを指定します．画像ファイルをドラッグ&ドロップでも設定できます．

  画像ファイルのパスは可能な限り相対パスで保存・管理されます．詳しくは[こちら](#パターン画像のファイルパスについて)．


### 角丸めσ

オブジェクトの出っ張った部分を円弧の形に丸めます．またオブジェクトの縁部分を透明ピクセルに削ることもできます．

![角丸めσ](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/20c148d9-49e6-4d6b-9a0b-282cbe7e5649)

- 手順的には，一度内側縁取りの計算をして今度はそれを打ち消すように外側に縁取りをすることで角の丸みを実現しています．縁取り計算が2回分必要なため，その分処理速度も遅くなります．

> [!TIP]
> 指定した半径の円盤がどこにも入りきらないような小さい / 狭いオブジェクトに対して適用すると，全て透明ピクセルになり見えなくなってしまいます．
>
> このフィルタ効果の初期値は `半径` が `32.0` になっているため，フィルタ効果を追加した直後にオブジェクトが表示されなくなってしまうことがあります．`半径` を小さくすれば見えるようになるので，適切な値に調整してください．


#### 各種パラメタ

![角丸めσの GUI](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/dda5b0a4-ff15-4240-893d-4fe1adc7145c)

- 半径

  丸めた角の円弧の半径をピクセル単位で指定します．大きすぎるとオブジェクト全体が見えなくなってしまうので注意してください．

  最小値は `0.0`, 最大値は `500.0`, 初期値は `32.0` です．

- 縁の縮小

  オブジェクトの縁部分を透明ピクセルに削ります．削る幅をピクセル単位で指定します．大きすぎるとオブジェクト全体が見えなくなってしまうので注意してください．

  `サイズも縮小` が ON の場合は，ここで指定したピクセル数だけオブジェクトの縦横サイズが縮みます．

  最小値は `0.0`, 最大値は `500.0`, 初期値は `0.0` です．

- 透明度

  `半径` や `縁の縮小`，`ぼかし` の効果で透明になるピクセルの透明度を % 単位で指定します．`0.0` だと何もしません．

  `サイズも縮小` が ON の場合はこのパラメタは無効です．

  最小値は `0.0`, 最大値は `100.0`, 初期値は `100.0` です．

- ぼかし

  オブジェクトの縁をぼかします．ぼかし幅の大きさをピクセル単位で指定します．大きすぎるとオブジェクト全体が見えなくなってしまうので注意してください．

  最小値は `0.0`, 最大値は `500.0`, 初期値は `0.0` です．

- αしきい値 / 基準α和

  `方式` で必要なパラメタを指定します．詳しくは[こちら](#方式-について)．

  最小値は `0.0`, 最大値は `100.0`, 初期値は `50.0`.

- サイズも縮小

  `縁の縮小` で指定したピクセル数だけ，オブジェクトの縦横サイズを縮めます．また，ON の場合 `透明度` の指定が無効になります．

  初期値は OFF.

- 方式

  利用する円形縁取りのアルゴリズムを指定します．詳しくは[こちら](#方式-について)．

  初期値は [`2値化倍精度`](#2値化倍精度).


### アウトラインσ

オブジェクトの縁部分を別オブジェクトとして抜き出して描画します．少し隙間の空いた縁取りにしたり，他のフィルタ効果を組み合わせるなどして自由に加工できます．`凸半径` や `凹半径` で曲率の調整もできます．

![アウトラインσ](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/ae140f91-d4a3-482b-a97b-8d15d1036a0b)

- 手順的には，縁取り計算を最大4回組み合わせて実現しています．その分処理速度も遅くなる傾向があります．

#### 各種パラメタ

![アウトラインσの GUI](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/974ec1b7-d881-437e-8c13-93e7c547918e)

- 距離

  アウトラインの位置を，オブジェクト境界からの距離としてピクセル単位で指定します．正の値だとオブジェクトの外側に，負の値だと内側に配置します．

  ここで指定した位置が「基準線」となって，アウトラインの片側になります．`凸半径` や `凹半径` で曲率の曲率指定もこの基準線に対して適用されます．アウトラインのもう片側はここを基準にして `ライン幅` 分だけ縁取り計算で決定します．

  最小値は `-500.0`, 最大値は `500.0`, 初期値は `10.0`.

- ライン幅

  アウトラインの幅をピクセル単位で指定します．正の値だと `距離` で指定した「基準線」から外側方向に，負の値だと内側方向に幅を広げます．`0.0` だとオブジェクトサイズの計算だけして非表示になります．

  ラインの形状ではなく，穴のない塗りつぶしをしたい場合は最小値の `-4000.0` を指定してください．

  最小値は `-4000.0`, 最大値は `500.0`, 初期値は `10.0`.

- 凸半径 / 凹半径

  `距離` で指定した基準線の曲率半径の最低保証値をピクセル単位で指定します．凸側（元オブジェクトの出っ張った部分）と凹側（元オブジェクトの凹んだ部分）とで独立に指定します．

  - 手順的には縁取り計算を複数回適用することで実現しています．これらに `0.0` 以外の値を指定すると最大で2回縁取り計算の回数が増えることになり，その分だけ遅くなります．

  最小値は `0.0`, 最大値は `500.0`, 初期値は `0.0`.

- ぼかし

  アウトラインにぼかしを適用します．ぼかし幅をピクセル単位で指定，ただし `ライン幅` （の絶対値）の約 0.5 倍が指定できる最大値です．

  最小値は `0.0`, 最大値は `500.0`, 初期値は `0.0`.

- αしきい値 / 基準α和

  `方式` で必要なパラメタを指定します．詳しくは[こちら](#方式-について)．

  最小値は `0.0`, 最大値は `100.0`, 初期値は `50.0`.

- 画像X / 画像Y

  `パターン画像ファイル` を設定している場合のみ有効，パターン画像ファイルの位置を調節します．

  画像ファイルの左上座標が，元のオブジェクトの左上座標と一致する配置が $(0,0)$ の基準点になります．

  最小値は `-4000`, 最大値は `4000`, 初期値は `0`.

- 方式

  利用する円形縁取りのアルゴリズムを指定します．詳しくは[こちら](#方式-について)．

  初期値は [`2値化倍精度`](#2値化倍精度).

- 手順

  `凸半径` と `凹半径` を指定している場合，縁取り操作を外方向や内方向に複数回かけることで曲率を操作しています．その際の縁取り操作の適用順序を指定します．凹凸の激しい図形の場合，この適用順序で結果が変わってくることがあります．

  選択肢は以下の通り:

  |手順|特徴|
  |---|:---|
  |`速い方`|`距離` の正負に応じて切り替えます．<br>`距離` を時間経過で変化させると `0.0` の前後で不連続になることがあります．|
  |`縮小→拡大→縮小`|アウトラインが全体的に小さめになります．|
  |`拡大→縮小→拡大`|アウトラインが全体的に大きめになります．|

  - `凸半径` と `凹半径` のどちらか一方が `0.0` の場合，この指定による影響はありません．
  - `距離` が正の場合は `縮小→拡大→縮小` が，負の場合は `拡大→縮小→拡大` が速い傾向があります．
  - 特にこだわりがなく，`距離` を時間変化させない / 変化させても `0.0` の前後をまたがない場合は `速い方` の指定で十分です．

- 縁色の設定

  アウトラインの色を指定します．パターン画像ファイルが設定されている場合は無視されます．

  初期値は `RGB( 255 , 255 , 255 )` （白）です．

- パターン画像ファイル

  アウトラインに適用するパターン画像ファイルを指定します．画像ファイルをドラッグ&ドロップでも設定できます．

  画像ファイルのパスは可能な限り相対パスで保存・管理されます．詳しくは[こちら](#パターン画像のファイルパスについて)．


## `方式` について

円形縁取りの計算アルゴリズムを指定します．5種類実装していて，境界付近の形状にそれぞれ特徴があり，得手不得手があります．計算速度にも差があります．場面に応じて最適なものを選んでください．

またこのアルゴリズムの選択によっては，追加のパラメタ `αしきい値` や `基準α和` を指定して調整できます．

初期値は[`2値化倍精度`](#2値化倍精度)になっています．

- おおよその計算時間の目安は以下の表の通りです．

  |方式|サイズ 300<br>縁取り幅 20|サイズ 1000<br>縁取り幅 500|サイズ 300 + ぼかし 20<br>縁取り幅 20|
  |:---:|:---:|:---:|:---:|
  |[2値化](#2値化)|470--500|6,600--8,000|650--750|
  |[2値化倍精度](#2値化倍精度)|510--540|7,600--8,600|730--840|
  |[総和](#総和)|740--840|150,000--160,000|1,100--1,300|
  |[最大値(安定)](#最大値安定)|1,120--1,300|395,000--410,000|1,790--1,900|
  |[最大値(高速)](#最大値高速)|800--1,000|84,000--90,000|2,700--3,000|
  |標準の縁取り|280--300|3,300--3,700|330--410|
  |縁取りT|10,900--11,400|518,000--524,000|20,000--20,600|

  - `図形` オブジェクトの `星型` (あるいはそれにフィルタ効果の `ぼかし` をかけたもの) に対して，各フィルタ効果の計算にかかった時間を計測．単位は $\mu \mathrm{s}$ (マイクロ秒).

  - 最後の2つは比較参考です．

    - 「標準の縁取り」は patch.aul による高速化を適用した状態，「縁取りT」は「高精度」のチェックを外した状態での計測です．

  - 自環境による大雑把な計測で，正確性に関しては保証はありません．あくまでも参考程度にとらえてください．


### 2値化

![2値化の適用例](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/43050195-27d3-4045-88da-abd20444c87c)

元画像をα値の大小によって透明部分と不透明部分の2種類に分類し，それに基づいて高速な手法で円形縁取りします．

追加のパラメタの `αしきい値` は，透明部分と不透明部分を分ける境界のα値を % 単位で指定します．初期値は `50.0`.

- 以下のような特徴があります:

  1.  実装されている中では最も高速なアルゴリズムです．平均・最悪計算時間ともに [Landau の記号](https://ja.wikipedia.org/wiki/%E3%83%A9%E3%83%B3%E3%83%80%E3%82%A6%E3%81%AE%E8%A8%98%E5%8F%B7)では理論上最速の $O(WH)$ です ($W$: 縁取り後の画像幅, $H$: 縁取り後の画像高さ).

  1.  出来上がりの画像は非常にジャギーになります．また一度2値化する関係上，伸ばし棒のようにα値のグラデーションで表現された緩やかな斜め線も階段状の見え方になります．

- 使える場面の例:

  1.  後でぼかしなどをかけるため，ジャギーであることがあまり問題にならないような場面．

  1.  ピクセル風のゲーム画面など，ジャギーな画像が馴染むような場面．


### 2値化倍精度

![2値化倍精度の適用例](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/a59deb4a-2b26-4de9-acbb-455fa110e07c)

[`2値化`](#2値化)と同様の計算手法ですが，縦横それぞれ2倍のサイズの画像に対して処理をした後半分に縮小，という操作に相当する計算をします．`2値化` よりも計算速度は落ちますが遜色ない速度で，簡易的なアンチエイリアスが実現できます．

追加のパラメタの `αしきい値` は，透明部分と不透明部分を分ける境界のα値を % 単位で指定します．初期値は `50.0`.

- 以下のような特徴があります:

  1.  実装されている中では2番目に高速なアルゴリズムです．平均・最悪計算時間も Landau の記号では `2値化` と同じく理論上最速の $O(WH)$ です．

  1.  出来上がりの画像は `2値化` に比べて軽くアンチエイリアスがかかったものになります．一方で `2値化` と同様に一度2値化する関係上，伸ばし棒のようにα値のグラデーションで表現された緩やかな斜め線も階段状の見え方になります．

- 使える・使えない場面の例:

  1.  多くの図形やイラストなどの枠線として遜色なく利用できると思います．

  1.  テキストの縁取りとしては，特に小さい文字だとα値の繊細なグラデーションによる斜め線などがうまく拾えず，ジャギーに見えてしまうこともあります．


### 総和

![総和の適用例](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/28a71c9a-afe4-4300-a7ca-26d0289cf8d4)

各点に対して，その点周辺の円形範囲にあるピクセルのα値の総和をとって，その和からα値を算出します．テキストオブジェクトの `縁取り文字` にも使われている手法で，α値の繊細な変化も拾うことができアンチエイリアスも綺麗です．

追加のパラメタの `基準α和` は，不透明度が 100% になるために必要なα値の総和の最小値をコントロールできます．`0.0` に近付くほどジャギーさが目立ちやすくなりますが境界がはっきりします．`100.0` に近付くほどジャギーさは目立ちにくくなりますが境界ぼやけます．初期値は `50.0`.

- 以下のような特徴があります:

  1.  実装されている中では遅めの分類のアルゴリズムです．平均・最悪計算時間はともに Landau の記号で $O(rWH)$ です ($r$: 円の半径, $W$: 縁取り後の画像幅, $H$: 縁取り後の画像高さ).

  1.  他のアルゴリズムと比べてとがった部分の縁取りサイズが小さくなります．一方で `2値化` とは違い，伸ばし棒のようにα値のグラデーションで表現された緩やかな斜め線も縁取りの境界に綺麗に反映されやすくなります．

- 使える・使えない場面の例:

  1.  テキストオブジェクトの `縁取り文字` の代わりに利用できます．同じ手法なのでほとんど見た目が変わらない上に，複数文字まとめて処理できるので軽くなる場面が多いです ([TIPS](#tips) も参照).

  1.  実装されている中では遅めのアルゴリズムなので，[角丸めσ](#角丸めσ)や[アウトラインσ](#アウトラインσ)など縁取りを複数回繰り返すフィルタ効果は輪をかけて遅くなりがちです．


### 最大値(安定)

![最大値の適用例](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/18a5157f-0af9-473a-a93b-5932a8962dee)

各点に対して，その点周辺の円形範囲にあるピクセルのα値の最大値を探して，その値を縁取り後のα値として採用します．元画像のα値がそのまま反映されるため，元画像のアンチエイリアスの具合がそのまま反映されます．

追加のパラメタ (`αしきい値` など) はありません．

- 以下のような特徴があります:

  1.  実装されている中では一番遅いアルゴリズムです ([`総和`](#総和) の3倍くらい遅い). 平均・最悪計算時間は `総和` と同じく Landau の記号で $O(rWH)$ です．

  1.  元画像のアンチエイリアスの具合がそのまま反映されるため，もともとα値のないオブジェクト（図形オブジェクトの四角形や，透過情報のない画像ファイルなど）だと [`2値化`](#2値化) と同様ジャギーな画像になります．

      図形オブジェクトの円などは，境界付近のグラデーションがそのまま反映されるため綺麗に見えやすいですが，星型のとがった部分などは階段状のグラデーションになってしまいます．

  1.  [`最大値(高速)`](#最大値高速) とは異なり，オブジェクトの形状によって極端に遅くなることはありません．

- 使える・使えない場面の例:

  1.  円オブジェクトや，円オブジェクトを描画に利用したカスタムオブジェクトなどは綺麗に見えることが多いです．

  1.  ぼかしのかかったオブジェクトも，ぼかし具合がそのまま反映されます．

  1.  ただ満遍なく綺麗なグラデーションが施されたオブジェクトはあまり多くないため，使える場面は非常に限られます．

### 最大値(高速)

[`最大値(安定)`](#最大値安定) と同じ計算結果ですが，多くの画像に対して高速に動作しますが，境界のぼやけた画像など一部極端に遅くなるものもあります．

追加のパラメタ (`αしきい値` など) はありません．

- 以下のような特徴があります:

  1.  実装されている中では遅めの分類のアルゴリズムです (ほとんどの場合 [`総和`](#総和) よりやや速い). 最悪計算時間は最も遅く，Landau の記号で $O(r^2WH)$ です．

  1.  結果の画像の質に関しては [`最大値(安定)`](#最大値安定) と全く同じです．

  1.  [`最大値(安定)`](#最大値安定) とは異なり，境界がぼやけた画像に対しては極端に重くなることがあります．

- 使える・使えない場面の例:

  1.  円オブジェクトや，円オブジェクトを描画に利用したカスタムオブジェクトなどは綺麗に見えることが多いですし，[`最大値(安定)`](#最大値安定) よりも高速です．

  1.  ぼかしをかけた画像など一部画像に対しては極端に重くなるため注意．


## パターン画像のファイルパスについて

パターン画像のファイルパスは可能な限りプロジェクトファイルか AviUtl.exe のあるフォルダからの相対パスとして記録管理するようにしています．

これで動画編集ファイルを整理する際プロジェクトファイルと素材ファイルを同じフォルダにまとめておけば，フォルダごと移動するだけでそのまま移動先のフォルダで編集や閲覧ができるようになります．

スクリプトなどで使う素材ファイルを AviUtl.exe のあるフォルダに置いている場合でも，AviUtl.exe をフォルダごと移動 / コピーでそのまま使うことができます．

スクリプトで指定した場合や `.exo`, `.exa`, `.exc` ファイルなどでも相対パスの形で指定できます．仕様は以下の通り:

1.  `<exe>` から始まる文字列は AviUtl.exe のあるフォルダからの相対パスとみなされます．

    例:

    - AviUtl.exe が `C:\hoge\aviutl.exe` に配置されていて，ファイルパスが `<exe>fuga\image.png` だった場合，`C:\hoge\fuga\image.png` のファイルをパターン画像として読み込みます．

1.  `<aup>` から始まる文字列はプロジェクトファイルのあるフォルダからの相対パスとみなされます．ただしプロジェクトファイルが未保存でファイル名が未指定の場合，AviUtl.exe のあるフォルダからの相対パスとして取り扱われます．

    例:

    - プロジェクトファイルが `C:\foo\bar.aup` に保存されていて，ファイルパスが `<aup>xyz\image.png` だった場合，`C:\foo\xyz\image.png` のファイルをパターン画像として読み込みます．

    - プロジェクトファイルが未保存の場合，AviUtl.exe が `C:\hoge\aviutl.exe` に配置されていて，ファイルパスが `<aup>xyz\image.png` だった場合，`C:\hoge\xyz\image.png` のファイルをパターン画像として読み込みます．

1.  その他の場合，普通にフルパス（あるいは現在実行中の AviUtl.exe の作業フォルダからの相対パス）として扱います．


## スクリプトでの利用について

スクリプト制御やアニメーション効果，カスタムオブジェクトなどでも利用できます．

- 例:

  ```lua
  obj.effect("縁取りσ","サイズ",10, "凹半径",20, "kind",1, "file",[[<aup>images\pattern1.png]])
  ```

  `サイズ` が `10`, `凹半径` が `20`, `方式` が `2値化倍精度`,
  `パターン画像ファイル` がプロジェクトファイルのあるフォルダからの相対パスで，
  `(.aup のあるフォルダ)\images\pattern1.png` の指定で `縁取りσ` を適用．

ほとんどのパラメタ名はトラックバーのボタンにあるテキストの通りですが，一部違うものがあります．

1.  `αしきい値` `基準α和` は `"param_a"` で指定します．

1.  [`方式`](#方式-について) は `"kind"` で指定します．

    対応は次の表の通りです:

    |アルゴリズム|数値|備考|
    |:---|---:|---|
    |`2値化`|`0`||
    |`2値化倍精度`|`1`|デフォルト|
    |`総和`|`2`||
    |`最大値(安定)`|`3`||
    |`最大値(高速)`|`4`||

1.  `縁色の設定` は `"color"` で指定します．

    `number` 型で `0xRRGGBB` の形式です．

1.  `パターン画像ファイル` は `"file"` で指定します．

    `<aup>` や `<exe>` を利用して相対パスでファイルを指定することもできます．[[詳細](#パターン画像のファイルパスについて)]

1.  `手順`（[アウトラインσ](#アウトラインσ)のみ）は `"order"` で指定します．

    対応は次の表の通りです:

    |手順|数値|備考|
    |:---|---:|---|
    |`速い方`|`0`|デフォルト|
    |`縮小→拡大→縮小`|`1`||
    |`拡大→縮小→拡大`|`2`||


## TIPS

1.  テキストオブジェクトの `縁取り文字` `縁取り文字(細)` `縁のみ` `縁のみ(細)` はアルゴリズムとしては[総和](#総和)に相当する方式です．これらの代わりに[縁取りσ](#縁取りσ)を使う場合，次のような違いがあります．

    1.  テキストオブジェクトは1文字ずつ縁取り処理します．そのため，隣接する文字の縁で文字が一部隠れることがあります．

        ![縁取り文字で被ってしまう例](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/144b26b8-4a03-450f-ae72-5ab7792d9f8f)

        一方，縁取りσは複数文字をまとめて処理します．そのため縁部分が文字に被ることはありません．

        ![縁取りσだと被らない](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/4b22805c-ed9a-4212-8622-6d1d25cdd8c3)

        - `文字毎に個別オブジェクト` の場合は縁取りσでも1文字ずつなので当てはまりません．

    1.  2--3文字以上のテキストなら多くの場面で縁取りσのほうが高速です．

        ただし制御文字 `<p>` などで，間に大きな空白ある場合などは遅くなることがあります．フォントの種類やサイズなどにもよります．また `文字毎に個別オブジェクト` の場合でも遅くなる傾向があります．

    1.  テキストオブジェクトの場合は制御文字 `<#>` で文字毎に縁取りの色を変えられます．また，`<s>` で文字毎にサイズを変えることで，縁取りのサイズも変化させられます．縁取りσの場合はどの文字にも一律で同じ色・サイズの縁取りしかかけられません．

    1.  テキストオブジェクトの場合はライン幅を2種類からしか選べませんが，縁取りσは任意サイズの指定ができます．

    1.  テキストオブジェクトの縁取り文字は，実は真円より左右に平均 0.5 ピクセルずつ横長になっています．一方で縁取りσは正確に真円なので，厳密に同じ表示にはできません．

        - 参考: https://github.com/nazonoSAUNA/patch.aul/blob/fbadd93fb505687c1258c13ed39cc0dc916d135a/patch/patch_fast_text_border.cpp#L498

          実数値から整数値への近似で `0.5` を加算してから `(int)` キャストしていますが（いわゆる四捨五入），真円のための正しい式は `0.5` を加算しないもの（いわゆる小数点以下切り捨て）です．その結果，左右に平均 0.5 ピクセルずつ「太った」縁取りになっています．

          このコードは patch.aul による高速化版ですが，patch.aul を適用しなくても縁取り文字でも全く同じ計算結果になったため，拡張編集に元々からある仕様のようです．


1.  角丸めσの `ぼかし` は拡張編集標準の境界ぼかしと類似の効果がありますが，凹んだ部分のある図形などでは境界ぼかしより綺麗に見えることがあります．

    ![境界ぼかしと角丸めσのぼかしの例](https://github.com/sigma-axis/aviutl_CircleBorder_S/assets/132639613/676c7d74-f626-4090-a356-0b292e7a447b)

1.  小さい穴がある図形で `凹半径` を大きくしていくと，一定のタイミングで穴が塞がれてしまいます．`凹半径` や[角丸めσ](#角丸めσ)など一部の機能は不連続的な変化になるので，時間変化などをさせる場合は不自然に見えないよう注意してください．


## 既知の問題

1.  [角丸めσ](#角丸めσ)で `半径` を `0.1` ずつ増やしていくと，角が単調に丸まっていくのではなく振動するように出っ張ったり引っ込んだりすることがあります．

    これは小数点以下のピクセルサイズに対して小数点以下切り捨て・切り上げによる誤差蓄積が原因です．現状見通しのいい解決案がないため，`半径` を時間変化させる際には整数単位で動かすなどの工夫が必要です．


## 改版履歴

- **v1.14** (2025-01-06)

  - 「アウトラインσ」で最大画像サイズを超えると例外が発生していたのを修正．

- **v1.13** (2024-11-23)

  - 「縁取りσ」で，`サイズ` が `0.1` 以上 `0.4` 以下で，画像が位置ずれしたり，ゴミ画像が残ることがあったのを修正．

- **v1.12** (2024-08-17)

  - 「縁取りσ」で，`透明度` が 100 など一部条件下で `内透明度` が反映されなかったのを修正．

- **v1.11** (2024-08-16)

  - 「縁取りσ」で `サイズ` が負でかつパターン画像を適用していた場合の次の挙動を修正:
    - パターン画像の位置が正しくなかった．
    - 内側縁取りが元オブジェクト全体を覆っているとき正しくパターンが適用されなかった．

  - 「アウトラインσ」で一部処理が並列実行されていなかったのを修正．

- **v1.10** (2024-08-15)

  - `方式` に `最大値(高速)` を追加．
    - 多くの場合で高速ですが，境界がぼやけた画像など場合によっては極端に遅くなります．
    - 既存の `最大値` は `最大値(安定)` に名前変更．
      - この名前変更による既存のプロジェクトファイルやスクリプトなどの互換性への影響はありません．

  - `方式` の `2値化` と `2値化倍精度` を大幅高速化 (1.3 -- 1.5 倍).

  - `方式` の `最大値(安定)` による縮小処理を高速化．

  - 「縁取りσ」で `サイズ` が正のときのα値合成計算を少し変更．

  - 全体が透明な画像に対して「縁取りσ」をかけると別の画像が現れたのを修正．

  - その他最適化 / コード整理 / バグ修正．

- **v1.02** (2024-07-08)

  - `ぼかし` のパラメタ指定時のぼかし処理を高速化．

  - `方式` に `総和`，`最大値` を選んだ時の処理を高速化．

  - マルチスレッドの並列化処理で環境によってはエラーになることがあったのを修正，無駄処理の原因になっていたコードも修正．

- **v1.01** (2024-07-07)

  - 画像サイズが大きい場合の `総和`，`最大値` のアルゴリズムで極端に重くなっていたのを改善．

- **v1.00** (2024-07-05)

  - 初版．


## ライセンス

このプログラムの利用・改変・再頒布等に関しては MIT ライセンスに従うものとします．

---

The MIT License (MIT)

Copyright (C) 2024-2025 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

https://mit-license.org/


#  Credits

##  aviutl_exedit_sdk

https://github.com/ePi5131/aviutl_exedit_sdk （利用したブランチは[こちら](https://github.com/sigma-axis/aviutl_exedit_sdk/tree/self-use)です．）

---

1条項BSD

Copyright (c) 2022
ePi All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
THIS SOFTWARE IS PROVIDED BY ePi “AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ePi BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#  連絡・バグ報告

- GitHub: https://github.com/sigma-axis
- Twitter: https://x.com/sigma_axis
- nicovideo: https://www.nicovideo.jp/user/51492481
- Misskey.io: https://misskey.io/@sigma_axis
- Bluesky: https://bsky.app/profile/sigma-axis.bsky.social

