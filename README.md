# RocketSimuMC
ロケットが弾道飛行した際の落下範囲をモンテカルロシミュレーションにより推定する．

### description
#### モンテカルロシミュレーション
モンテカルロ法とはシミュレーションや数値計算を乱数を用いて行う手法([Wikipedia](https://ja.wikipedia.org/wiki/%E3%83%A2%E3%83%B3%E3%83%86%E3%82%AB%E3%83%AB%E3%83%AD%E6%B3%95))．
[ロケットシミュレーション](https://github.com/Jirouken/ModelRocketSimulator)において確率的に生じる風向・風速の揺らぎを乱数によって再現し，多数シミュレーションを行うことで落下範囲の分布を得る．
