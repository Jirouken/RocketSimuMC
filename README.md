# RocketSimuMC
ロケットが弾道飛行した際の落下範囲をモンテカルロシミュレーションにより推定する．

## description
#### モンテカルロシミュレーション
モンテカルロ法とはシミュレーションや数値計算を乱数を用いて行う手法([Wikipedia](https://ja.wikipedia.org/wiki/%E3%83%A2%E3%83%B3%E3%83%86%E3%82%AB%E3%83%AB%E3%83%AD%E6%B3%95))．
[ロケットシミュレーション](https://github.com/Jirouken/ModelRocketSimulator)において確率的に生じる風向・風速の揺らぎを乱数によって再現し，多数シミュレーションを行うことで落下範囲の分布を得る．

#### 風速モデル
シミュレーションで用いる風速モデルは[べき法則](https://www.rikanenpyo.jp/kaisetsu/kisyo/kisyo_011.html)で表される．
これにガウス分布から生成されるノイズを加える(multiplicative noise)．
この時，ガウス分布の平均は0，標準偏差は高度に比例する．
![風速ノイズモデル](https://github.com/Jirouken/RocketSimuMC/blob/master/wind_noise01001.png)
さらに，風向にもガウス分布から生成されるノイズを加える(additive noise)．

## usage
~~~
from RocketSimuMC import RocketSimuMC
import pandas as pd

# 設計データ読み込み
design = pd.read_csv('design.csv')
# 打ち上げ条件読み込み
condition = pd.read_csv('condition.csv')

# モンテカルロシミュレーション
rs = RocketSimuMC()
rs.initialize(design.loc[0], condition.loc[2])
fr, ma = rs.falling_range(n=1000)

print(ma.mean(), ma.std())

plt.figure(figsize=(8, 8))
plt.scatter(fr[:, 0], fr[:, 1])
plt.xlabel('x', fontsize=20)
plt.ylabel('y', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()
~~~
