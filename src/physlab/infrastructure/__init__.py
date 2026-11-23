"""外界との接点（描画・ファイル出力）。結果データ（application）を受け取って外へ出す。

matplotlib や ffmpeg を使ってよいのはこの層より外側だけ（依存ルールは
docs/design/003-layered-architecture.md）。
"""
