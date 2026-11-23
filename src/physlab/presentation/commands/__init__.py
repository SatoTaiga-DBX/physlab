"""シミュレーションごとのコマンド。

規約（これを満たすモジュールは ``physlab list`` / ``physlab run`` に自動で現れる）:

* ``NAME``  : 実行時に指定する短い名前（例: ``"pendulum"``）
* ``TITLE`` : 一覧に表示する日本語タイトル
* ``main(argv: list[str] | None = None) -> None`` : 引数を解釈して実行する

トピック別にサブパッケージを切る: ``mechanics``、将来 ``waves``、``fluids``。
"""
