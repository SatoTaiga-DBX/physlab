"""シミュレーションごとのコマンド。

規約（これを満たすモジュールは ``physlab list`` / ``physlab run`` に自動で現れる）:

* ``NAME``  : 実行時に指定する短い名前（例: ``"pendulum"``）
* ``TITLE`` : 一覧に表示する日本語タイトル
* ``main(argv: list[str] | None = None) -> int | None`` : 引数を解釈して実行する。
  整数を返すと終了コードになる（``None`` は 0）。

任意（教科書の図）:

* ``figures() -> list[Path]`` : 教科書に貼る静止画を ``outputs/figures/`` に書き出し、
  パスの一覧を返す。動画は含めず、数秒で終わるようにする。**引数を取らず**、常に教科書の条件で作る
  （``main`` の引数を変えても、教科書の図が上書きされないように）。
  ``physlab run figures`` と ``physlab run textbook`` が全コマンド分をまとめて呼ぶ。

トピック別にサブパッケージを切る: ``mechanics``、将来 ``waves``、``fluids``。
トピックに属さない道具（``figures``・``textbook``）は直下に置く。
"""
