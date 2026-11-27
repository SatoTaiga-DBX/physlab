# 知見ログ

やってみて分かったこと・ハマったことを、題材ごとではなく**テーマごと**に貯める。
同じテーマの発見は新しいファイルを作らず、既存のファイルに追記する。

## 書くときのルール

- 冒頭に出典（どの設計書・どの実行で分かったか）と実測日を書く。
- **実測した数値**と**推測・未検証**を分けて書く（「未検証の見通し」の節を設ける）。
- 結論だけでなく、どのパラメータで観測したかを残す（再現できるように）。
- 一覧は [../README.md](../README.md) の表に 1 行足す。

## 一覧

- [numerical-integration.md](numerical-integration.md) — 数値積分：手法ごとの精度とエネルギー保存
- [visualization.md](visualization.md) — 可視化：matplotlib で動画を作るときの落とし穴
- [architecture.md](architecture.md) — アーキテクチャ：層構造・import-linter・テスト配置を運用して分かったこと
- [analytical-mechanics.md](analytical-mechanics.md) — 解析力学：最小作用の原理・一般座標・ガリレイ変換の実測値
- [latex-textbook.md](latex-textbook.md) — TeX 教科書：LuaLaTeX（MiKTeX）の設定・警告・ハマった点
- [copyright.md](copyright.md) — 著作権：参考書をもとに教科書を書くときの調査（条文・判例）と方針
- [conservation-laws.md](conservation-laws.md) — 保存法則：対称性を壊す実験、積分法と保存量、力学的相似の実測値
- [collisions.md](collisions.md) — 衝突と散乱：崩壊・弾性衝突・断面積・ラザフォード・微小角の実測値とハマった点
